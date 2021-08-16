#!/usr/bin/env python

import os
import glob
import platform
import re
import shutil
import subprocess
import sys
import sysconfig
from setuptools.command.build_ext import build_ext
from setuptools.command.install import install
from distutils.command.clean import clean
from setuptools import setup, Extension
from distutils.command.build import build
from distutils.command.sdist import sdist
from distutils.errors import DistutilsExecError
from distutils.version import StrictVersion
from distutils.util import convert_path
from distutils.sysconfig import customize_compiler
from distutils.ccompiler import show_compilers


# Path to the directory that contains this setup.py file.
base_dir = os.path.abspath(os.path.dirname(__file__))


def find_package_version(package_name):
    try:
        return __import__(package_name).__version__
    except ImportError:
        return None


def is_package_installed(package_name):
    try:
        __import__(package_name)
        return True
    except ImportError:
        return False


def in_conda():
    return os.path.exists(os.path.join(sys.prefix, 'conda-meta'))


def find_version():
    """Extract the current version of AutoDock Vina.
    
    The version will be obtained from (in priority order):
    1. version.py (file created only when using GitHub Actions)
    2. git describe
    3. __init__.py (as failback)

    """
    version_file = os.path.join(base_dir, 'vina', 'version.py')
    if os.path.isfile(version_file):
        with open(version_file) as f:
            version = f.read().strip()
        
        print('Version found: %s (from version.py)' % version)
        return version

    try:
        git_output = subprocess.check_output(['git', 'describe', '--abbrev=7', '--dirty', '--always', '--tags'])
        git_output = git_output.strip().decode()

        if git_output.startswith('v'):
            git_output = git_output[1:]
        version = git_output.replace('dirty', 'mod').replace('-', '+', 1).replace('-', '.')

        print('Version found %s (from git describe)' % version)
        return version
    except:
        pass
    
    init_file = os.path.join(base_dir, 'vina', '__init__.py') 
    with open(init_file) as f:
        for line in f:
            version_match = re.match(r'^__version__ = "(.+?)"$', line)

            if version_match:
                version = version_match.group(1)

                print('Version found %s (from __init__.py)' % version)
                return version
    
    raise RuntimeError('Could not find version string for AutoDock Vina.')


def execute_command(cmd_line):
    """Simple function to execute bash command."""
    args = cmd_line.split()
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    output, errors = p.communicate()
    return output, errors


def locate_ob():
    """Try use pkgconfig to locate Open Babel, otherwise guess default location."""
    # Warn if the (major, minor) version of the installed OB doesn't match these python bindings
    if not is_package_installed("openbabel"):
        raise RuntimeError("Error: Openbabel is not installed.")

    py_ver = StrictVersion(find_package_version('openbabel'))
    py_major_ver, py_minor_ver = py_ver.version[:2]

    if in_conda:
        # It means that Openbabel was installed in an Anaconda env
        data_pathname = sysconfig.get_path('data')
        include_dirs = data_pathname + os.path.sep + 'include' + os.path.sep + 'openbabel{}'.format(py_major_ver)
        library_dirs = data_pathname + os.path.sep + 'lib'

        if os.path.isdir(include_dirs):
            print('Open Babel location automatically determined in Anaconda.')
            return include_dirs, library_dirs
        else:
            print("Warning: We are in an Anaconda env, but Openbabel is not installed here.")

    pcfile = 'openbabel-{}'.format(py_major_ver)
    output, errors = execute_command("pkg-config --modversion %s" % pcfile)

    if output:
        # It means that Openbabel was install with apt-get
        ob_ver = StrictVersion(output.strip())

        if not ob_ver.version[:2] == py_ver.version[:2]:
            print('Warning: Open Babel {}.{}.x is required. Your version ({}) may not be compatible.'
                    .format(py_major_ver, py_minor_ver, ob_ver))
        include_dirs = execute_command("pkg-config --variable=pkgincludedir %s" % pcfile)[0].strip()
        library_dirs = execute_command("pkg-config --variable=libdir %s" % pcfile)[0].strip()

        print('Open Babel location automatically determined by pkg-config.')
        return include_dirs, library_dirs
    else:
        pathnames = ['/usr', '/usr/local']
        for pathname in pathnames:
            include_dirs = pathname + os.path.sep + 'include' + os.path.sep + 'openbabel{}'.format(py_major_ver)
            library_dirs = pathname + os.path.sep + 'lib'

            if os.path.isdir(include_dirs):
                print('Open Babel location was automatically guessed.')
                return include_dirs, library_dirs

        print('Open Babel location was set to default location.')
        return include_dirs, library_dirs


def locate_boost():
    """Try to locate boost."""
    if in_conda:
        data_pathname = sysconfig.get_path('data')
        include_dirs = data_pathname + os.path.sep + 'include'
        library_dirs = data_pathname + os.path.sep + 'lib'
        
        if os.path.isdir(include_dirs + os.path.sep + 'boost'):
            print('Boost library location automatically determined in this conda environment.')
            return include_dirs, library_dirs
        else:
            print('Boost library is not installed in this conda environment.')

    include_dirs = '/usr/local/include'

    if os.path.isdir(include_dirs + os.path.sep + 'boost'):
        if glob.glob('/usr/local/lib/x86_64-linux-gnu/libboost*'):
            return include_dirs, '/usr/local/lib/x86_64-linux-gnu'
        elif glob.glob('/usr/local/lib64/libboost*'):
            return include_dirs, '/usr/local/lib64'
        elif glob.glob('/usr/local/lib/libboost*'):
            return include_dirs, '/usr/local/lib'

    include_dirs = '/usr/include'

    if os.path.isdir(include_dirs + os.path.sep + 'boost'):
        if glob.glob('/usr/lib/x86_64-linux-gnu/libboost*'):
            return include_dirs, '/usr/lib/x86_64-linux-gnu'
        elif glob.glob('/usr/lib64/libboost*'):
            return include_dirs, '/usr/lib64'
        elif glob.glob('/usr/lib/libboost*'):
            return include_dirs, '/usr/lib'

    return None, None


class CustomBuild(build):
    """Ensure build_ext runs first in build command."""
    def run(self):
        # Fix to make it compatible with wheel
        # We copy src directory only when it is not present already.
        # If it is present, it means that we are creating linux wheels using the manylinux docker image
        # The src copy is done outside setup.py before we start creating the wheels
        # Source: https://github.com/pypa/pip/issues/3500
        if not os.path.exists('src'):
            shutil.copytree('../../src', 'src')

        self.run_command('build_ext')
        build.run(self)


class CustomInstall(install):
    """Ensure build_ext runs first in install command."""
    def run(self):
        # This is not called when creating wheels for linux in the docker image
        if not os.path.exists('src'):
            shutil.copytree('../../src', 'src')

        self.run_command('build_ext')
        install.run(self)

        # It means that we are in build/python, so we can safely remove src
        if not os.path.exists('build/python'):
            shutil.rmtree('src')


class CustomSdist(sdist):
    """Add swig interface files into distribution from parent directory."""
    def make_release_tree(self, base_dir, files):
        sdist.make_release_tree(self, base_dir, files)
        link = 'hard' if hasattr(os, 'link') else None
        pkg_dir = os.path.join(base_dir, 'vina')
        self.copy_file(os.path.join('vina', 'vina.i'), pkg_dir, link=link)

    def run(self):
        if not os.path.exists('src'):
            shutil.copytree('../../src', 'src')

        sdist.run(self)

        # It means that we are in build/python, so we can safely remove src
        if not os.path.exists('build/python'):
            shutil.rmtree('src')


class CustomBuildExt(build_ext):
    """Custom build_ext to set SWIG options and print a better error message."""
    def finalize_options(self):
        # Setting include_dirs, library_dirs, swig_opts here instead of in Extension constructor allows them to be
        # overridden using -I and -L command line options to python setup.py build_ext.
        build_ext.finalize_options(self)

        # Boost
        self.boost_include_dir, self.boost_library_dir = locate_boost()

        if self.boost_include_dir is None and self.boost_library_dir is None:
            error_msg = 'Boost library location was not found!\n'
            error_msg += 'Directories searched: conda env, /usr/local/include and /usr/include.'
            raise ValueError(error_msg)
        else:
            print('Boost library location was automatically guessed at %s.' % self.boost_include_dir)

        self.include_dirs.append(self.boost_include_dir)
        self.library_dirs.append(self.boost_library_dir)

        # Openbabel
        #self.ob_include_dir, self.ob_library_dir = locate_ob()
        #self.include_dirs.append(self.ob_include_dir)
        #self.library_dirs.append(self.ob_library_dir)

        # Vina
        self.include_dirs.append('src/lib')
        # SWIG
        # shadow, creates a pythonic wrapper around vina
        # castmode
        self.swig_opts = ['-c++', '-small', '-naturalvar', '-fastdispatch', '-shadow', '-py3']
        self.swig_opts += ['-I%s' % i for i in self.include_dirs]

        print('- include_dirs: %s\n- library_dirs: %s' % (self.include_dirs, self.library_dirs))
        print('- swig options: %s' % self.swig_opts)
        print('- libraries: %s' % self.libraries)

    def swig_sources(self, sources, extension):
        try:
            return build_ext.swig_sources(self, sources, extension)
        except DistutilsExecError:
            print('\nError: SWIG failed.',
                  'You may need to manually specify the location of Open Babel include and library directories. '
                  'For example:',
                  '  python setup.py build_ext -I{} -L{}'.format(self.include_dirs, self.library_dir),
                  '  python setup.py install',
                  sep='\n')
            sys.exit(1)

    def build_extensions(self):
        customize_compiler(self.compiler)

        # Patch for macOS (libboost_thread)
        # Check if we have an include "system"
        include_system = set(self.include_dirs).intersection(['/usr/local/include', '/usr/include'])
        if platform.system() == 'Darwin':
            if include_system:
                # In a conda env on macOS there is no -mt suffix at -lboost_thread
                idx = self.extensions[0].extra_link_args.index('-lboost_thread')
                self.extensions[0].extra_link_args[idx] = '-lboost_thread-mt'
            # To get the right @rpath on macos for libraries
            self.extensions[0].extra_link_args.append('-Wl,-rpath,' + self.library_dirs[0])
            self.extensions[0].extra_link_args.append('-Wl,-rpath,' + '/usr/lib')
        
        print('- extra link args: %s' % self.extensions[0].extra_link_args)

        # Replace current compiler to g++
        self.compiler.compiler_so[0] = "g++"
        self.compiler.compiler_so.insert(2, "-shared")

        remove_flags = ["-Wstrict-prototypes", "-Wall"]
        for remove_flag in remove_flags:
            try:
                self.compiler.compiler_so.remove(remove_flag)
            except ValueError:
                print('Warning: compiler flag %s is not present, cannot remove it.' % remove_flag)
                pass

        self.compiler.compiler_so.append("-std=c++11")
        self.compiler.compiler_so.append("-Wno-long-long")
        self.compiler.compiler_so.append("-pedantic")
        # Source: https://stackoverflow.com/questions/9723793/undefined-reference-to-boostsystemsystem-category-when-compiling
        self.compiler.compiler_so.append('-DBOOST_ERROR_CODE_HEADER_ONLY')

        print('- compiler options: %s' % self.compiler.compiler_so)
        build_ext.build_extensions(self)


# Dirty fix when the compilation is not done in build/python
# Fix for readthedocs
autodock_swig_interface = 'vina/autodock_vina.i'
package_dir = {}
if os.path.exists('build/python'):
    autodock_swig_interface = 'build/python/' + autodock_swig_interface
    package_dir = {'vina': 'build/python/vina'}


obextension = Extension(
    'vina._vina_wrapper',
    sources=['src/lib/random.cpp', 'src/lib/utils.cpp', 'src/lib/vina.cpp',
             'src/lib/quaternion.cpp', 'src/lib/monte_carlo.cpp', 'src/lib/non_cache.cpp',
             'src/lib/mutate.cpp', 'src/lib/szv_grid.cpp', 'src/lib/quasi_newton.cpp',
             'src/lib/parallel_progress.cpp', 'src/lib/model.cpp', 'src/lib/coords.cpp',
             'src/lib/ad4cache.cpp', 'src/lib/grid.cpp', 'src/lib/parallel_mc.cpp', 
             'src/lib/conf_independent.cpp', 'src/lib/parse_pdbqt.cpp',
             'src/lib/cache.cpp', autodock_swig_interface],
    extra_link_args=['-lboost_thread', '-lboost_serialization',
                     '-lboost_filesystem', '-lboost_program_options'],
    #libraries=['openbabel'],
)


setup(
    name='vina',
    version=find_version(),
    author='Diogo Santos Martins, Jerome Eberhardt, Andreas F. Tillack, Stefano Forli',
    author_email='forli@scripps.edu',
    license='Apache-2.0',
    url='https://ccsb.scripps.edu/',
    description='Python interface to AutoDock Vina',
    long_description=open(os.path.join(base_dir, 'README.md')).read(),
    long_description_content_type="text/markdown",
    zip_safe=False,
    cmdclass={'build': CustomBuild, 'build_ext': CustomBuildExt, 'install': CustomInstall, 'sdist': CustomSdist},
    packages=['vina'],
    package_dir=package_dir,
    install_requires=['numpy>=1.18'],
    python_requires='>=3.5.*',
    ext_modules=[obextension],
    #entry_points={"console_scripts": ["vina = vina.vina_cli:main"]},
    classifiers=[
        'Environment :: Console',
        'Environment :: Other Environment',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        #'Operating System :: Microsoft :: Windows',
        'Operating System :: OS Independent',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: C++',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Software Development :: Libraries'
    ]
)

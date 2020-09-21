#!/usr/bin/env python

import os
import glob
import re
import subprocess
import sys
import sysconfig
from setuptools.command.build_ext import build_ext
from setuptools.command.install import install
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
    """Extract the current version of these python bindings from the __init__.py file."""
    try:
        with open(os.path.join(base_dir, 'vina', '__init__.py')) as fp:
            for line in fp:
                version_match = re.match(r'^__version__ = "(.+?)"$', line)
                if version_match:
                    return version_match.group(1)
            raise Exception('Could not find version string in vina/__init__.py.')
    except IOError:
        #return '1.2'
        raise Exception('Could not find vina/__init__.py.')


def execute_command(cmd_line):
    """Simple function to execute bash command."""
    args = cmd_line.split()
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    output, errors = p.communicate()
    return output, errors


def locate_ob():
    """Try use pkgconfig to locate Open Babel, otherwise guess default location."""
    # Warn if the (major, minor) version of the installed OB doesn't match these python bindings
    assert is_package_installed("openbabel"), print("Error: Openbabel is not installed.")

    py_ver = StrictVersion(find_package_version('openbabel'))
    py_major_ver, py_minor_ver = py_ver.version[:2]

    if in_conda:
        # It means that Openbabel was installed in an Anaconda env
        data_pathname = sysconfig.get_path('data')
        include_dirs = data_pathname + os.path.sep + 'include' + os.path.sep + 'openbabel{}'.format(py_major_ver)
        library_dirs = data_pathname + os.path.sep + 'lib'

        if os.path.isdir(include_dirs):
            print('Open Babel location automatically determined in Anaconda:')
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

        print('Open Babel location automatically determined by pkg-config:')
        return include_dirs, library_dirs
    else:
        pathnames = ['/usr', '/usr/local']
        for pathname in pathnames:
            include_dirs = pathname + os.path.sep + 'include' + os.path.sep + 'openbabel{}'.format(py_major_ver)
            library_dirs = pathname + os.path.sep + 'lib'

            if os.path.isdir(include_dirs):
                print('Open Babel location was automatically guessed:')
                return include_dirs, library_dirs

        print('Open Babel location was set to default location:')
        return include_dirs, library_dirs


def locate_boost():
    return "/usr/local/include", "/usr/local/lib"


class CustomBuild(build):
    """Ensure build_ext runs first in build command."""
    def run(self):
        self.run_command('build_ext')
        build.run(self)


class CustomInstall(install):
    """Ensure build_ext runs first in install command."""
    def run(self):
        self.run_command('build_ext')
        install.run(self)


class CustomSdist(sdist):
    """Add swig interface files into distribution from parent directory."""
    def make_release_tree(self, base_dir, files):
        sdist.make_release_tree(self, base_dir, files)
        link = 'hard' if hasattr(os, 'link') else None
        pkg_dir = os.path.join(base_dir, 'vina')
        self.copy_file(os.path.join('.', 'vina-python.i'), pkg_dir, link=link)


class CustomBuildExt(build_ext):
    """Custom build_ext to set SWIG options and print a better error message."""
    def finalize_options(self):
        # Setting include_dirs, library_dirs, swig_opts here instead of in Extension constructor allows them to be
        # overridden using -I and -L command line options to python setup.py build_ext.
        build_ext.finalize_options(self)
        
        # Boost
        self.boost_include_dir, self.boost_library_dir = locate_boost()
        self.include_dirs.append(self.boost_include_dir)
        self.library_dirs.append(self.boost_library_dir)
        # Openbabel
        self.ob_include_dir, self.ob_library_dir = locate_ob()
        self.include_dirs.append(self.ob_include_dir)
        self.library_dirs.append(self.ob_library_dir)
        # Vina
        self.include_dirs.append('../../src/lib')
        # SWIG
        # shadow, creates a pythonic wrapper around vina
        # castmode
        self.swig_opts = ['-c++', '-small', '-naturalvar', '-fastdispatch', '-shadow', '-py3', '-threads']
        self.swig_opts += ['-I%s' % i for i in self.include_dirs]
        
        print('- include_dirs: %s\n- library_dirs: %s' % (self.include_dirs, self.library_dirs))
        print('- swig options: %s' % self.swig_opts)
        print('- libraries: %s' % self.libraries)

    def swig_sources(self, sources, extension):
        try:
            return build_ext.swig_sources(self, sources, extension)
        except DistutilsExecError:
            print('\nError: SWIG failed. Is Open Babel installed?',
                  'You may need to manually specify the location of Open Babel include and library directories. '
                  'For example:',
                  '  python setup.py build_ext -I{} -L{}'.format(self.ob_include_dir, self.ob_library_dir),
                  '  python setup.py install',
                  sep='\n')
            sys.exit(1)

    def build_extensions(self):
        customize_compiler(self.compiler)

        try:
            self.compiler.compiler_so[0] = "g++"
            self.compiler.compiler_so.insert(2, "-shared")
            self.compiler.compiler_so.remove("-Wstrict-prototypes")
            self.compiler.compiler_so.remove("-Wall")
            self.compiler.compiler_so.append("-std=c++11")
            self.compiler.compiler_so.append("-Wno-long-long")
            self.compiler.compiler_so.append("-pedantic")
        except (AttributeError, ValueError):
            pass

        print('- compiler options: %s' % self.compiler.compiler_so)
        build_ext.build_extensions(self)


obextension = Extension(
    'vina._vina',
    sources=glob.glob('../../src/lib/*.cpp') + ['vina/autodock_vina.i'],
    extra_link_args=['-lboost_system', '-lboost_thread', 
                     '-lboost_serialization', '-lboost_filesystem', 
                     '-lboost_program_options'],
    libraries=['openbabel'],
    headers=glob.glob('../../src/lib/*.h'),
)


setup(
    name='vina',
    version=find_version(),
    author='The beautiful people',
    author_email='people@scripps.edu',
    license='Apache-2.0',
    url='https://ccsb.scripps.edu/',
    description='Python interface to AutoDock Vina',
    long_description=open(os.path.join(base_dir, 'README.rst')).read(),
    zip_safe=False,
    cmdclass={'build': CustomBuild, 'build_ext': CustomBuildExt, 'install': CustomInstall, 'sdist': CustomSdist},
    packages=['vina'],
    ext_modules=[obextension],
    classifiers=[
        'Environment :: Console',
        'Environment :: Other Environment',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: OS Independent',
        'Operating System :: POSIX',
        'Operating System :: POSIX :: Linux',
        'Operating System :: Unix',
        'Programming Language :: C++',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Software Development :: Libraries'
    ]
)

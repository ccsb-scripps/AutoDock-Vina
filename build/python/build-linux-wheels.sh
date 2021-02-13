#!/bin/bash

# cd <Autodock-Vina_directory>
# DOCKER_IMAGE=quay.io/pypa/manylinux2014_x86_64
# PLAT="manylinux2014_x86_64"
# sudo docker run --rm -it -e PLAT=$PLAT -v "$(pwd)":/io "$DOCKER_IMAGE" /io/build/python/build-linux-wheels.sh

set -e -u -x

function repair_wheel {
    wheel="$1"
    if ! auditwheel show "$wheel"; then
        echo "Skipping non-platform wheel $wheel"
    else
        auditwheel repair "$wheel" --plat "$PLAT" -w /io/build/python/wheelhouse/
    fi
}


# Install a system package required by our library
yum install -y cmake wget python-devel

# Install SWIG
# Source: http://www.linuxfromscratch.org/blfs/view/svn/general/swig.html
wget https://downloads.sourceforge.net/swig/swig-4.0.2.tar.gz
tar -xvf swig-4.0.2.tar.gz 
cd swig-4.0.2
./configure --prefix=/usr --without-maximum-compile-warnings && make
make install && install -v -m755 -d /usr/share/doc/swig-4.0.2 && cp -v -R Doc/* /usr/share/doc/swig-4.0.2
cd /
rm -rf swig-4.0.2*

# Install Boost
# Source: http://www.linuxfromscratch.org/blfs/view/svn/general/boost.html
wget https://dl.bintray.com/boostorg/release/1.75.0/source/boost_1_75_0.tar.bz2
tar -xvf boost_1_75_0.tar.bz2 
cd boost_1_75_0
./bootstrap.sh --prefix=/usr --with-python=python
./b2 install threading=multi link=shared
cd /
rm -rf boost_1_75_0*

cd /io/build/python

# We need to copy the src file before setuptools create a tmp dir of it
# Source: https://github.com/pypa/pip/issues/3500
cp -r ../../src src

# Compile wheels
for PYBIN in /opt/python/*/bin; do
    "${PYBIN}/pip" install -r requirements.txt
    "${PYBIN}/pip" wheel . --no-deps -w wheelhouse/
done

rm -rf src

# Bundle external shared libraries into the wheels
for whl in wheelhouse/vina*.whl; do
    repair_wheel "$whl"
done

# Install packages and test
for PYBIN in /opt/python/*/bin/; do
    "${PYBIN}/pip" install vina --no-index -f wheelhouse
#    (cd "$HOME"; "${PYBIN}/nosetests" pymanylinuxdemo)
done

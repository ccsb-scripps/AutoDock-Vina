#!/bin/bash

# brew install boost
# Py35: compatibility issue during conda installation with py35, so version ignored
# Py36: issue during compilation with py36, so version ignored
# conda create --name cp37 python=3.7 pip wheel swig numpy delocate twine
# conda create --name cp38 python=3.8 pip wheel swig numpy delocate twine
# conda create --name cp39 python=3.9 pip wheel swig numpy delocate twine
# cd <AutoDock-Vina>/build/python

eval "$(conda shell.bash hook)"
conda activate base

mkdir wheelhouse

PYENVS=('cp37' 'cp38' 'cp39')

for PYENV in "${PYENVS[@]}"; do
	echo "${PYENV}"
	conda activate $PYENV
	python setup.py clean --all bdist_wheel
	delocate-listdeps --all dist/vina*"${PYENV}"*.whl
	delocate-wheel -w wheelhouse -v dist/vina*"${PYENV}"*.whl
	delocate-listdeps --all wheelhouse/vina*"${PYENV}"*.whl
done

rm -rf build dist vina.egg-info

conda activate base

%define DOCSTRING
"#################################################################\n\
# If you used AutoDock Vina in your work, please cite:          #\n\
#                                                               #\n\
# O. Trott, A. J. Olson,                                        #\n\
# AutoDock Vina: improving the speed and accuracy of docking    #\n\
# with a new scoring function, efficient optimization and       #\n\
# multithreading, Journal of Computational Chemistry 31 (2010)  #\n\
# 455-461                                                       #\n\
#                                                               #\n\
# DOI 10.1002/jcc.21334                                         #\n\
#                                                               #\n\
# Please see http://vina.scripps.edu for more information.      #\n\
#################################################################\n"
%enddef

%module(docstring=DOCSTRING, package="vina") vina_wrapper

%begin %{
#define SWIG_PYTHON_2_UNICODE
//#define SWIG_FILE_WITH_INIT
%}

%{
#include "ad4cache.h"
#include "array3d.h"
#include "atom.h"
#include "atom_base.h"
#include "atom_constants.h"
#include "atom_type.h"
#include "bfgs.h"
#include "brick.h"
#include "cache.h"
#include "common.h"
#include "conf.h"
#include "convert_substring.h"
#include "conf_independent.h"
#include "coords.h"
#include "curl.h"
#include "file.h"
#include "grid.h"
#include "grid_dim.h"
#include "igrid.h"
#include "incrementable.h"
#include "int_pow.h"
#include "macros.h"
#include "matrix.h"
#include "model.h"
#include "monte_carlo.h"
#include "mutate.h"
#include "non_cache.h"
#include "parallel.h"
#include "parallel_mc.h"
#include "parallel_progress.h"
#include "parse_error.h"
#include "parse_pdbqt.h"
#include "potentials.h"
#include "precalculate.h"
#include "quasi_newton.h"
#include "quaternion.h"
#include "random.h"
#include "scoring_function.h"
#include "szv_grid.h"
#include "tree.h"
#include "triangular_matrix_index.h"
#include "utils.h"
%}

// Set and reset dlopenflags so that plugin loading works fine for "import _openbabel"
%pythonbegin %{
import sys
if sys.platform.find("linux") != -1:
    dlflags = sys.getdlopenflags()
    import ctypes
    sys.setdlopenflags(dlflags | ctypes.RTLD_GLOBAL)
%}
%pythoncode %{
if sys.platform.find("linux") != -1:
    sys.setdlopenflags(dlflags)
%}

// Add standard C++ library
%include "std_array.i"
%include "std_list.i"
%include "std_map.i"
%include "std_vector.i"
%include "std_string.i"

// Help SWIG to understand some special types, like list of strings
namespace std {
    %template(IntVector) vector<int>;
    %template(DoubleVector) vector<double>;
    %template(DoubleVectorVector) vector<vector<double>>;
    %template(StringVector) vector<string>;
    %template(ConstCharVector) vector<const char*>;
    //%template(OBMolVector) vector<OpenBabel::OBMol*>;
}

// Add numpy
//%include "numpy.i"
//%init %{
//import_array();
//%}

%include "vina.i"

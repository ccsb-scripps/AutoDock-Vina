%{
#include "vina.h"
%}

%include "vina.h"

%extend Vina {
    // Add docstring to python functions
    %feature("autodoc", "3");
    int test() {
        return 1;
    }
};
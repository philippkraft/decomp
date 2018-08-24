%feature("compactdefaultargs");
%feature("autodoc","1") ;


// Include typemaps for STL
%include <std_string.i>
%include <std_vector.i>

// enable exception support
%include "exception.i"
%exception {
	try {
		$action
	} 
	SWIG_CATCH_STDEXCEPT
}
%include "attribute.i"

%module decomp

%{
#include "SOMcomponent.h"
#include "SOM.h"
%}

%include "SOMcomponent.h"

%extend DECOMP::SOMcomponent {
std::string __repr__() {return $self->Name;}
}

%echo "SOMcomponent ok, now SOM..."

%template(component_set) std::vector<DECOMP::SOMcomponent>;
%attribute(DECOMP::SOM, double, C, get_C_pool);
%attribute(DECOMP::SOM, double, CN, get_CN);

%include "SOM.h"


%extend DECOMP::SOM {
	double __getitem__(const SOMcomponent& comp)
	{
		return $self->get_C_pool(comp.Id);
	}
	void __setitem__(const SOMcomponent& comp, double pool_size)
	{
		(*$self)[comp] = pool_size;
	}
	SOM __rmul__(double right)
	{
	    return (*$self) * right;
	}
	std::string __repr__()
	{ return $self->to_string();}
	%pythoncode
	{
    def __iter__(self):
        pools=SOM.get_pool_types()
        for pool in pools:
            yield pool, self[pool]
	}
}


%pythoncode {
	EDC, CELL, LIGN, RC, DOC, CO2 = SOM.get_pool_types()
}


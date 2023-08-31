#ifndef DEF_BOOMAES
#define DEF_BOOMAES

#include <vector>
#include <map>
#include "/opt/gurobi/linux64/include/gurobi_c++.h"
#include "Matrix.hpp"


class mycallback: public GRBCallback
{
public:
	Matrix mat;
	unsigned numvars;
	GRBVar * varsD;
	int version;

	mycallback(int xnumvars, GRBVar* xvars, Matrix const & eqs, int xversion) {
		numvars = xnumvars;
		varsD = xvars;
		mat = eqs;
		version = xversion;
	}
protected:
	void callback();
};

#endif

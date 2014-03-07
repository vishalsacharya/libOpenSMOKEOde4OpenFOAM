
#include "ode/OpenSMOKE_OdeSystemObject.h"
#include "math/stiff-ode-solvers/OdeSystemObject.h"
#include "math/stiff-ode-solvers/StiffOdeSolverObject_Dense.h"
#include "ODESystem.H"

#ifndef OpenSMOKEOde_H
#define OpenSMOKEOde_H

class OpenSMOKEOde
{
public:

	OpenSMOKEOde(const Foam::ODESystem& openFOAM_odeSystem);

	unsigned int NumberOfEquations() const { return number_of_equations_; }

	int Equations(const double t, const Eigen::VectorXd& y, Eigen::VectorXd& dy);
	int Print(const double t, const Eigen::VectorXd& y);

private:

	unsigned int number_of_equations_;
	Foam::scalarField yFoam_;
	Foam::scalarField dyFoam_;


	const Foam::ODESystem& openFOAM_odeSystem_;
};

#include "OpenSMOKEOde.C"

#endif // OpenSMOKEOde_H


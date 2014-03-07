#include <Eigen/Dense>
#include "ode/OpenSMOKE_OdeSystemObject.h"
#include "ode/OpenSMOKE_OdeSystemObject.h"
#include "math/stiff-ode-solvers/OdeSystemObject.h"
#include "math/stiff-ode-solvers/StiffOdeSolverObject_Dense.h"
#include "ODESystem.H"
#include "OpenSMOKEOde.H"

#ifndef OdeSystem_OpenSMOKE_H
#define OdeSystem_OpenSMOKE_H

class OdeSystem_OpenSMOKE : public OpenSMOKE::UserOdeSystemObject 
{
public:

	OdeSystem_OpenSMOKE( OpenSMOKEOde& cell ) :
	cell_(cell)
	{
	}
	
	virtual void GetSystemFunctions(const Eigen::VectorXd &y, const double t, Eigen::VectorXd &dy)
	{
		cell_.Equations(t, y, dy);
	}

	void MyPrint(const Eigen::VectorXd &y, const double t)
	{
		cell_.Print(t, y);
	}

private:

	OpenSMOKEOde& cell_;;
};

#endif // OdeSystem_OpenSMOKE_H


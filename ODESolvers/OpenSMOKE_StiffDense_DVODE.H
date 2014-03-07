/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::OpenSMOKE_DVODE

Description
    Semi-implicit Euler ODE solver of order (0)1.

    The method calculates the new state from:
    \f[
        y_{n+1} = y_n
          + \delta_x\left[I - \delta_x\frac{\partial f}{\partial y}\right]^{-1}
            \cdot \left[f(y_n) + \delta_x\frac{\partial f}{\partial x}\right]
    \f]
    The error is estimated directly from the change in the solution,
    i.e. the difference between the 0th and 1st order solutions:
    \f[
        err_{n+1} = y_{n+1} - y_n
    \f]

SourceFiles
    OpenSMOKE_DVODE.C

\*---------------------------------------------------------------------------*/

#ifndef OpenSMOKE_StiffDense_DVODE_H
#define OpenSMOKE_StiffDense_DVODE_H

#include "ODESolver.H"
#include "ODESystem.H"

// OpenSMOKE
#include <Eigen/Dense>
//#include "OpenSMOKEOde.H"

//#if DVODE_SUPPORT == 1
#include "ode/OpenSMOKE_DVODE_Interface.h"
#include "ode/OpenSMOKE_DVODE.h"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class ODESystem_DVODE : public OpenSMOKE::OpenSMOKE_OdeSystemObject
{
	DEFINE_ODESOLVERINTERFACE_DVODE(ODESystem_DVODE)

	const Foam::ODESystem* openFOAM_odeSystem_;
	
	Foam::scalarField yFoam_;
	Foam::scalarField dyFoam_;

public:

	void Set(const Foam::ODESystem* openFOAM_odeSystem)
	{
		openFOAM_odeSystem_ = openFOAM_odeSystem;
		yFoam_.resize(openFOAM_odeSystem_->nEqns());
		dyFoam_.resize(openFOAM_odeSystem_->nEqns());
	}
 
	int GetSystemFunctions(const double t, double* y,  double* dy)
	{
		for(unsigned int i=0;i<openFOAM_odeSystem_->nEqns();i++)
			yFoam_[i] = y[i];
		openFOAM_odeSystem_->derivatives(t,yFoam_,dyFoam_);
		for(unsigned int i=0;i<openFOAM_odeSystem_->nEqns();i++)
			dy[i] = dyFoam_[i];

		return(0);
	}

	int GetAnalyticalJacobian(const double t,  double* y,  double* J)
	{
		return(0);
	}
 
	int GetWriteFunction(const double t, double *y)
	{
		return 0;
	}
}; 
COMPLETE_ODESOLVERINTERFACE_DVODE(ODESystem_DVODE)

namespace Foam
{


/*---------------------------------------------------------------------------*\
                           Class OpenSMOKE_StiffDense_DVODE Declaration
\*---------------------------------------------------------------------------*/

class OpenSMOKE_StiffDense_DVODE
:
    public ODESolver
{
    // Private data

	Eigen::VectorXd yMin_;
	mutable Eigen::VectorXd yf_;
	mutable Eigen::VectorXd y0_;

	ODESystem_DVODE* ode_system_;

public:

    //- Runtime type information
    TypeName("OpenSMOKE_StiffDense_DVODE");


    // Constructors

        //- Construct from ODE
        OpenSMOKE_StiffDense_DVODE(const ODESystem& ode, const dictionary& dict);

    // Member Functions

	void solve
	(
	    const scalar xStart,
	    const scalar xEnd,
	    scalarField& y,
	    scalar& dxTry
	) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "OpenSMOKE_StiffDense_DVODE.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

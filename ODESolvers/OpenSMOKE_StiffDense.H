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
    Foam::OpenSMOKE_StiffDense

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
    OpenSMOKE_StiffDense.C

\*---------------------------------------------------------------------------*/

#ifndef OpenSMOKE_StiffDense_H
#define OpenSMOKE_StiffDense_H

// ODESolver (standard OpenFOAM)
#include "ODESolver.H"

// OpenSMOKE
#include <Eigen/Dense>
#include "OpenSMOKEOde.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class OpenSMOKE_StiffDense Declaration
\*---------------------------------------------------------------------------*/

class OpenSMOKE_StiffDense
:
    public ODESolver
{
    // Private data

	Eigen::VectorXd yMin_;
	mutable Eigen::VectorXd yf_;
	mutable Eigen::VectorXd y0_;

	OpenSMOKEOde* ode_system_;

public:

    //- Runtime type information
    TypeName("OpenSMOKE_StiffDense");


    // Constructors

        //- Construct from ODE
        OpenSMOKE_StiffDense(const ODESystem& ode, const dictionary& dict);

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

#include "OpenSMOKE_StiffDense.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

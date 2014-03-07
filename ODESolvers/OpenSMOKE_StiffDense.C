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

\*---------------------------------------------------------------------------*/

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(OpenSMOKE_StiffDense, 0);
    addToRunTimeSelectionTable(ODESolver, OpenSMOKE_StiffDense, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::OpenSMOKE_StiffDense::OpenSMOKE_StiffDense(const ODESystem& ode, const dictionary& dict)
:
    ODESolver(ode, dict)
{
	yMin_.resize(n_);
	yf_.resize(n_);
	y0_.resize(n_);
	
	for(label i=0;i<n_;i++)
		yMin_(i) = 0.;

	ode_system_ = new OpenSMOKEOde(odes_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::OpenSMOKE_StiffDense::solve
(
    const scalar xStart,
    const scalar xEnd,
    scalarField& y,
    scalar& dxTry
) const
{
	for(register label i=0;i<n_;i++)
		y0_(i) = y[i];

	OdeSystem_OpenSMOKE ode_object(*ode_system_);
	OpenSMOKE::StiffOdeSolverObject_Dense o( y0_, xStart, &ode_object, "Eigen" );
	
	o.SetMinimumConstraints(yMin_);
	o.SetTollAbs(absTol_[0]);
	o.SetTollRel(relTol_[0]);

	yf_ = o(xEnd,xEnd);

	for(register label i=0;i<n_;i++)
		y[i] = yf_(i);	

	return;
}


// ************************************************************************* //

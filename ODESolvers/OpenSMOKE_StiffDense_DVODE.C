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

#include "ode/OpenSMOKE_DVODE_Interface.h"
#include "ode/OpenSMOKE_DVODE.h"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(OpenSMOKE_StiffDense_DVODE, 0);
    addToRunTimeSelectionTable(ODESolver, OpenSMOKE_StiffDense_DVODE, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::OpenSMOKE_StiffDense_DVODE::OpenSMOKE_StiffDense_DVODE(const ODESystem& ode, const dictionary& dict)
:
    ODESolver(ode, dict)
{
	yMin_.resize(n_);
	yf_.resize(n_);
	y0_.resize(n_);
	
	for(label i=0;i<n_;i++)
		yMin_(i) = 0.;

	ode_system_ = ODESystem_DVODE::GetInstance();
	ode_system_->Set(&odes_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::OpenSMOKE_StiffDense_DVODE::solve
(
    const scalar xStart,
    const scalar xEnd,
    scalarField& y,
    scalar& dxTry
) const
{
	for(register label i=0;i<n_;i++)
		y0_(i) = y[i];

	OpenSMOKE::OpenSMOKE_DVODE<ODESystem_DVODE>   o(ode_system_);
	o.SetDimensions(n_);
	o.SetAbsoluteTolerance(absTol_[0]);
	o.SetRelativeTolerance(relTol_[0]);
	o.SetMaximumNumberOfSteps(10000);
	o.SetAnalyticalJacobian(false);
	o.SetInitialValues(xStart, y0_.data());
	o.Solve(xEnd);	o.Solution(yf_.data());

	for(register label i=0;i<n_;i++)
		y[i] = yf_(i);	

	return;
}


// ************************************************************************* //

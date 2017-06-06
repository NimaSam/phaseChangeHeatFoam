/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "Lee.H"
#include "fvc.H"//add
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseChangeTwoPhaseMixtures
{
    defineTypeNameAndDebug(Lee, 0);
    addToRunTimeSelectionTable
    (
        phaseChangeTwoPhaseMixture,
        Lee,
        components
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeTwoPhaseMixtures::Lee::Lee
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    const word& alpha1Name
)
:
    phaseChangeTwoPhaseMixture(typeName, U, phi, alpha1Name),
    rv_(phaseChangeTwoPhaseMixtureCoeffs_.lookup("rv")),
	rc_(phaseChangeTwoPhaseMixtureCoeffs_.lookup("rc")),
	T0_("T0", TSatLocal().dimensions(), 0.0)
{
    correct();
}
// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::phaseChangeTwoPhaseMixtures::Lee::mDotAlphal() const
{
    const volScalarField& T = alpha1_.db().lookupObject<volScalarField>("T");

    return Pair<tmp<volScalarField> >
    (
    	-rc_*rho2_*min(T-TSatLocal(),T0_)/TSatLocal()
    	,
    	-rv_*rho1_*max(T-TSatLocal(),T0_)/TSatLocal()
    );

}

Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::phaseChangeTwoPhaseMixtures::Lee::mDotP() const
{
    const volScalarField& T = alpha1_.db().lookupObject<volScalarField>("T");
    const volScalarField& p = alpha1_.db().lookupObject<volScalarField>("p");
    volScalarField limitedAlpha1(min(max(alpha1_, scalar(0)), scalar(1)));

    return Pair<tmp<volScalarField> >
    (
         -rc_*rho2_*min(T-TSatLocal(),T0_)/TSatLocal()
         *pos(p-pSat_)/max(p-pSat_,1E-06*pSat_)
         *(1.0 - limitedAlpha1)
    	,
       	 -rv_*rho1_*max(T-TSatLocal(),T0_)/TSatLocal()
    	*neg(p-pSat_)/max(pSat_-p,1E-06*pSat_)
    	*limitedAlpha1
    );
}

Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::phaseChangeTwoPhaseMixtures::Lee::mDotT() const
{
    const volScalarField& T = alpha1_.db().lookupObject<volScalarField>("T");
    volScalarField limitedAlpha1(min(max(alpha1_, scalar(0)), scalar(1)));


    return Pair<tmp<volScalarField> >
    (
    	-rc_*rho2_*neg(T-TSatLocal())/TSatLocal()
    	*(1.0 - limitedAlpha1)
    	,
    	rv_*rho1_*pos(T - TSatLocal())/TSatLocal()
    	*limitedAlpha1
    );
}


void Foam::phaseChangeTwoPhaseMixtures::Lee::correct()
{}



bool Foam::phaseChangeTwoPhaseMixtures::Lee::read()
{
	if (phaseChangeTwoPhaseMixture::read())

    {
        phaseChangeTwoPhaseMixtureCoeffs_ = subDict(type() + "Coeffs");
        phaseChangeTwoPhaseMixtureCoeffs_.lookup("rv") >> rv_;
        phaseChangeTwoPhaseMixtureCoeffs_.lookup("rc") >> rc_;

        return true;
    }
    else
    {
        return false;
    }

}

// ************************************************************************* //

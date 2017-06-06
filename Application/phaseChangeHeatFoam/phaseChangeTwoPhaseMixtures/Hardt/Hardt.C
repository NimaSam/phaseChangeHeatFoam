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

#include "Hardt.H"
#include "fvc.H"//add
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseChangeTwoPhaseMixtures
{
    defineTypeNameAndDebug(Hardt, 0);
    addToRunTimeSelectionTable
    (
        phaseChangeTwoPhaseMixture,
        Hardt,
        components
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeTwoPhaseMixtures::Hardt::Hardt
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    const word& alpha1Name
)
:
    phaseChangeTwoPhaseMixture(typeName, U, phi, alpha1Name),

    rv_(phaseChangeTwoPhaseMixtureCoeffs_.lookup("rv")),
    rc_(phaseChangeTwoPhaseMixtureCoeffs_.lookup("rc")),
    Cv_(phaseChangeTwoPhaseMixtureCoeffs_.lookup("Cv")),
    T0_("T0", TSat().dimensions(),0.0),// Temperature criteria
   	Cm1_("Cm1_", 2.0*Cv_*Hfg_*rho2_/((2.0-Cv_)*pow(2.0*M_PI*R_,0.5)))
{
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
Foam ::volScalarField Foam::phaseChangeTwoPhaseMixtures::Hardt::AbyV() const
{
	volScalarField limitedAlpha1 = min(max(alpha1_, scalar(0)), scalar(1));
	return(
			mag(fvc::grad(limitedAlpha1))
			);
}
Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::phaseChangeTwoPhaseMixtures::Hardt::mDotAlphal() const
{
    const volScalarField& T = alpha1_.db().lookupObject<volScalarField>("T");

    return Pair<tmp<volScalarField> >
    (
    	-rc_*Cm1_*min(T - TSatLocal() ,T0_)*AbyV()/sqrt(pow(TSatLocal(),3.0))
    	//*neg(T - TSatLocal() + 1E-05*TSat())
    	,
        -rv_*Cm1_*max(T - TSatLocal() ,T0_)*AbyV()/sqrt(pow(TSatLocal(),3.0))
        //*pos(T - TSatLocal() - 1E-05*TSat())
    );

}
Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::phaseChangeTwoPhaseMixtures::Hardt::mDotP() const
{
    const volScalarField& T = alpha1_.db().lookupObject<volScalarField>("T");
    const volScalarField& p = alpha1_.db().lookupObject<volScalarField>("p");
    volScalarField limitedAlpha1 = min(max(alpha1_, scalar(0)), scalar(1));

    return Pair<tmp<volScalarField> >
    (
    		-rc_*Cm1_*min(T - TSatLocal(),T0_)*AbyV()/sqrt(pow(TSatLocal(),3.0))
    		*pos(p-pSat_)/max(p-pSat_,1E-6*pSat_)
    		*(1.0-limitedAlpha1)
    	//	*neg(T - TSatLocal() + 1E-06*TSat())
    		,
    		-rv_*Cm1_*max(T - TSatLocal(),T0_)*AbyV()/sqrt(pow(TSatLocal(),3.0))
    		*neg(p-pSat_)/max(pSat_-p,1E-05*pSat_)
    		*limitedAlpha1
    	//	*pos(T - TSatLocal() - 1E-05*TSat())
    );
}

Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::phaseChangeTwoPhaseMixtures::Hardt::mDotT() const
{
    const volScalarField& T = alpha1_.db().lookupObject<volScalarField>("T");
    volScalarField limitedAlpha1 = min(max(alpha1_, scalar(0)), scalar(1));

    return Pair<tmp<volScalarField> >
    (
           -rc_*Cm1_*AbyV()/sqrt(pow(TSatLocal(),3.0))
            *neg(T - TSatLocal())
            *(1.0-limitedAlpha1)
         //   *neg(T - TSatLocal() + 1E-05*TSat())
            ,
            rv_*Cm1_*AbyV()/sqrt(pow(TSatLocal(),3.0))
           *limitedAlpha1
           *pos(T - TSatLocal())
       //    *pos(T - TSatLocal() - 1E-05*TSat())

    );
}

void Foam::phaseChangeTwoPhaseMixtures::Hardt::correct()
{}


bool Foam::phaseChangeTwoPhaseMixtures::Hardt::read()
{
    if (phaseChangeTwoPhaseMixture::read())
    {
        phaseChangeTwoPhaseMixtureCoeffs_ = subDict(type() + "Coeffs");
        phaseChangeTwoPhaseMixtureCoeffs_.lookup("Cv") >> Cv_;
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

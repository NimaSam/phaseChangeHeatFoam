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

#include "phaseChangeTwoPhaseMixture.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseChangeTwoPhaseMixture, 0);
    defineRunTimeSelectionTable(phaseChangeTwoPhaseMixture, components);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeTwoPhaseMixture::phaseChangeTwoPhaseMixture
(
    const word& type,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const word& alpha1Name
)
:
    twoPhaseMixtureI(U, phi, alpha1Name),//add
    phaseChangeTwoPhaseMixtureCoeffs_(subDict(type + "Coeffs")),
    pSat_(lookup("pSat")),
    TSat_(lookup("TSat")),
    TSatLocal_(readBool(lookup("TSatLocal"))),
    Hfg_(lookup("Hfg")),
    R_(lookup("R"))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::phaseChangeTwoPhaseMixture::vDotAlphal() const
{

	volScalarField alphalCoeff(1.0/rho1() - alpha1_*(1.0/rho1() - 1.0/rho2()));
	    Pair<tmp<volScalarField> > mDotAlphal = this->mDotAlphal();

	    return Pair<tmp<volScalarField> >
	    (
	        alphalCoeff*mDotAlphal[0],
	        alphalCoeff*mDotAlphal[1]
	    );
}
Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::phaseChangeTwoPhaseMixture::vDotP() const
{
	 dimensionedScalar pCoeff(1.0/rho1() - 1.0/rho2());
	 Pair<tmp<volScalarField> > mDotP = this->mDotP();

	 return Pair<tmp<volScalarField> >(pCoeff*mDotP[0], pCoeff*mDotP[1]);

}

Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::phaseChangeTwoPhaseMixture::vDotT() const
{
	   volScalarField rhoCp =  rho1()*C1()*alpha1_ + rho2()*C2()*(1.0-alpha1_);
	   //volScalarField TCoeff = (Hfg()+(C1()-C2())*TSat_)/rhoCp;
	   volScalarField TCoeff = Hfg()/rhoCp;
	   //volScalarField TCoeff = Hfg()*(alpha1_/C1() + (1.0 - alpha1_)/C2());
	   Pair<tmp<volScalarField> > mDotT = this->mDotT();

	    return Pair<tmp<volScalarField> >
	    (
	    		TCoeff*mDotT[0],
	    		TCoeff*mDotT[1]
	    );
}

Foam::volScalarField
Foam::phaseChangeTwoPhaseMixture::TSatLocal() const
{
	if (TSatLocal_)
	{
		//Info <<"TSatlocal" << endl;
	    const volScalarField& p = alpha1_.db().lookupObject<volScalarField>("p");
	    volScalarField oneByT = 1.0/TSat() - R_/Hfg_*log(max(p/pSat_,1E-08));

	    return volScalarField
	    (
	           1.0/oneByT
	    );
	}
	else
	{
		//Info <<"TSat" << endl;
		volScalarField one
		(
			IOobject
			(
				"one",
				alpha1_.mesh(),
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			alpha1_.mesh(),
			dimensionedScalar ("one",dimless, 1.0)
	    );

		return volScalarField
		(
			one*TSat_
	    );
	}
}


bool Foam::phaseChangeTwoPhaseMixture::read()
{
    if (twoPhaseMixtureI::read())
    {
        phaseChangeTwoPhaseMixtureCoeffs_ = subDict(type() + "Coeffs");
        lookup("pSat") >> pSat_;
        lookup("TSat") >> TSat_;
        lookup("TSatLocal") >> TSatLocal_;
        lookup("Hfg") >> Hfg_;
        lookup("R") >> R_;
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //

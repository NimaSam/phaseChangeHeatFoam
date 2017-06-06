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

Class
    Foam::phaseChangeTwoPhaseMixture

Description

SourceFiles
    phaseChangeTwoPhaseMixture.C
    newPhaseChangeModel.C

\*---------------------------------------------------------------------------*/

#ifndef phaseChangeTwoPhaseMixture_H
#define phaseChangeTwoPhaseMixture_H

#include "twoPhaseMixtureI/twoPhaseMixtureI.H"//add change all twoPhaseMixture into twoPhaseMixtureI
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "Pair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class phaseChangeTwoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class phaseChangeTwoPhaseMixture
:
    public twoPhaseMixtureI //add
{

protected:

    // Protected data

        dictionary phaseChangeTwoPhaseMixtureCoeffs_;

        //- Saturation vapour Pressure
        dimensionedScalar pSat_; //add
        //- Saturation vapour Temperature
        dimensionedScalar TSat_; //add
        //- true if Saturation vapour Temperature based on pressure
        bool TSatLocal_; //add
        //- latent Heat
        dimensionedScalar Hfg_;//add
        //- gas specific constant
        dimensionedScalar R_;//add

    // Private Member Functions
        //- Disallow copy construct
        phaseChangeTwoPhaseMixture(const phaseChangeTwoPhaseMixture&);

        //- Disallow default bitwise assignment
        void operator=(const phaseChangeTwoPhaseMixture&);


public:

    //- Runtime type information
    TypeName("phaseChangeTwoPhaseMixture");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            phaseChangeTwoPhaseMixture,
            components,
            (
                const volVectorField& U,
                const surfaceScalarField& phi,
                const word& alpha1Name
            ),
            (U, phi, alpha1Name)
        );


    // Selectors

        //- Return a reference to the selected phaseChange model
        static autoPtr<phaseChangeTwoPhaseMixture> New
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            const word& alpha1Name = "alpha1"
        );


    // Constructors

        //- Construct from components
        phaseChangeTwoPhaseMixture
        (
            const word& type,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const word& alpha1Name = "alpha1"
        );


    // Destructor

        virtual ~phaseChangeTwoPhaseMixture()
        {}


    // Member Functions
        //- Return const-access to the saturation vapour Temperature
        const dimensionedScalar& TSat() const //add
        {
        	return	TSat_;
        }
        //- Return const-access to the saturation vapour Temperature
        const dimensionedScalar& pSat() const //add
        {
        	return pSat_;
        }
        //- Return const-access to the latent heat
        const dimensionedScalar& Hfg() const //add
        {
        	return Hfg_;
        }

        //- Return the mass condensation and vaporisation rates as a
        //  coefficient to multiply (1 - alphal) for the condensation rate
        //  and a coefficient to multiply  alphal for the vaporisation rate
        virtual Pair<tmp<volScalarField> > mDotAlphal() const = 0;


        //- Return the mass condensation and vaporisation rates as an
        //  explicit term for the condensation rate and a coefficient to
        //  multiply (p - pSat) for the vaporisation rate
        virtual Pair<tmp<volScalarField> > mDotP() const = 0;


        virtual Pair<tmp<volScalarField> > mDotT() const = 0;

        //- Return the volumetric condensation and vaporisation rates as a
        //  coefficient to multiply (1 - alphal) for the condensation rate
        //  and a coefficient to multiply  alphal for the vaporisation rate
        Pair<tmp<volScalarField> > vDotAlphal() const;


        //- Return the volumetric condensation and vaporisation rates as an
        //  explicit term for the condensation rate and a coefficient to
        //  multiply (p - pSat) for the vaporisation rate
        Pair<tmp<volScalarField> > vDotP() const;


        Pair<tmp<volScalarField> > vDotT() const;


        //- Return const-access to the saturation vapour Temperature
        volScalarField TSatLocal() const;


        //- Correct the phaseChange model
        virtual void correct() = 0;

        //- Read the transportProperties dictionary and update
        virtual bool read() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Foam::twoPhaseMixture

Description
    A two-phase incompressible transportModel
    + adding ( c (specific heat), k (heat conductivity) , D (diffusity)

SourceFiles
    twoPhaseMixtureI.C


\*---------------------------------------------------------------------------*/

#ifndef twoPhaseMixtureI_H
#define twoPhaseMixtureI_H

#include "incompressible/transportModel/transportModel.H"
#include "incompressible/viscosityModels/viscosityModel/viscosityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class twoPhaseMixtureI Declaration
\*---------------------------------------------------------------------------*/

class twoPhaseMixtureI
:
    public transportModel
{
protected:

    // Protected data


        word phase1Name_;
        word phase2Name_;

        autoPtr<viscosityModel> nuModel1_;
        autoPtr<viscosityModel> nuModel2_;

        dimensionedScalar rho1_;
        dimensionedScalar rho2_;

        dimensionedScalar K1_;
        dimensionedScalar K2_;

        dimensionedScalar C1_;
        dimensionedScalar C2_;

        const volVectorField& U_;
        const surfaceScalarField& phi_;

        const volScalarField& alpha1_;

        volScalarField nu_;


    // Private Member Functions

        //- Calculate and return the laminar viscosity
        void calcNu();


public:

    // Constructors

        //- Construct from components
        twoPhaseMixtureI
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            const word& alpha1Name = "alpha1"
        );


    // Destructor

        ~twoPhaseMixtureI()
        {}


    // Member Functions

        const word phase1Name() const
        {
            return phase1Name_;
        }

        const word phase2Name() const
        {
            return phase2Name_;
        }

        //- Return const-access to phase1 viscosityModel
        const viscosityModel& nuModel1() const
        {
            return nuModel1_();
        }

        //- Return const-access to phase2 viscosityModel
        const viscosityModel& nuModel2() const
        {
            return nuModel2_();
        }

        //- Return const-access to phase1 density
        const dimensionedScalar& rho1() const
        {
            return rho1_;
        }

        //- Return const-access to phase2 density
        const dimensionedScalar& rho2() const
        {
            return rho2_;
        };

        //- Return const-access to phase1 specific heat
        const dimensionedScalar& C1() const //add ci
        {
            return C1_;
        }

        //- Return const-access to phase2 specific heat
        const dimensionedScalar& C2() const
        {
            return C2_;
        };
        //- Return const-access to phase1 heat conductivity
        const dimensionedScalar& K1() const //add ki
        {
        	return K1_;
        }

        //- Return const-access to phase2 heat conductivity
        const dimensionedScalar& K2() const
        {
        	return K2_;
        };

        //- Return the dynamic laminar viscosity
        tmp<volScalarField> mu() const;

        //- Return the face-interpolated dynamic laminar viscosity
        tmp<surfaceScalarField> muf() const;

        //- Return the kinematic laminar viscosity
        virtual tmp<volScalarField> nu() const
        {
            return nu_;
        }

        //- Return the face-interpolated dynamic laminar viscosity
        tmp<surfaceScalarField> nuf() const;

        //- Correct the laminar viscosity
        virtual void correct()
        {
            calcNu();
        }

        //- Read base transportProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

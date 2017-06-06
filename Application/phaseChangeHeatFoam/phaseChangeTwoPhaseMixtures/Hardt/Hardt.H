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
    under the terms of the GNU Generac License as published by the
    Free Software Foundation; either  2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the ho it will be useful, but WITHOUT
    ANY WARRANTY; without even the imarranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.he GNU General Public License
    for more details.

    You should have received a copy oNU General Public License
    along with OpenFOAM; if not, write Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floorn, MA 02110-1301 USA

Class
    Foam::phaseChangeTwoPhaseMixtures::Hardt

Description
    Hardt evaporation model.

    Reference:

SourceFiles
    Hardt.C

\*--------------------------------------------------------------------*/

#ifndef Hardt_H
#define Hardt_H

#include "phaseChangeTwoPhaseMixture.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseChangeTwoPhaseMixtures
{

/*--------------------------------------------------------------------*\
                              Class Hardt
\*--------------------------------------------------------------------*/

class Hardt
:
    public phaseChangeTwoPhaseMixture
{
    // Private data

    	dimensionedScalar rv_;

    	dimensionedScalar rc_;

    	dimensionedScalar Cv_;

        //temperature criteria
        dimensionedScalar T0_;

        dimensionedScalar Cm1_;

      // interface Area/Cell volume * normalization factor
        volScalarField AbyV() const;

public:

    //- Runtime type information
    TypeName("Hardt");


    // Constructors

        //- construct from components
        Hardt
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            const word& alpha1Name = "alpha1"
        );


    // Destructor

        virtual ~Hardt()
        {}


    // Member Functions
        dimensionedScalar R() {

        	return R_;
        }


        //multiply  alphal for the vaporization rate
        virtual Pair<tmp<volScalarField> > mDotAlphal() const;


        //add in p.Eqn for the vaporization rate
        virtual Pair<tmp<volScalarField> > mDotP() const;

        //add in T.Eqn for the vaporization rate
        virtual Pair<tmp<volScalarField> > mDotT() const;

        //- Correct the Hardt phaseChange model
        virtual void correct();

        //- Read the transportProperties dictionary and update
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace phaseChangeTwoPhaseMixtures
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

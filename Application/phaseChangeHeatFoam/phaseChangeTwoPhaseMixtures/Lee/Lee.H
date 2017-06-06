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
    Foam::phaseChangeTwoPhaseMixtures::Lee

Description
    Lee evaporation model.

    Reference:

SourceFiles
    Lee.C

\*--------------------------------------------------------------------*/

#ifndef Lee_H
#define Lee_H

#include "phaseChangeTwoPhaseMixture.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseChangeTwoPhaseMixtures
{

/*--------------------------------------------------------------------*\
                              Class Lee
\*--------------------------------------------------------------------*/

class Lee
:
    public phaseChangeTwoPhaseMixture
{
    // Private data

    	//mass transfer intensity factor [s^-1]
        dimensionedScalar rv_;
        //mass transfer intensity factor [s^-1]
        dimensionedScalar rc_;

        dimensionedScalar T0_;


public:

    //- Runtime type information
    TypeName("Lee");


    // Constructors

        //- construct from components
        Lee
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            const word& alpha1Name = "alpha1"
        );


    // Destructor

        virtual ~Lee()
        {}


    // Member Functions

        //multiply  alphal for the vaporization rate
        virtual Pair<tmp<volScalarField> > mDotAlphal() const;


        //add in p.Eqn for the vaporization rate
        virtual Pair<tmp<volScalarField> > mDotP() const;

        //add in T.Eqn for the vaporization rate
        virtual Pair<tmp<volScalarField> > mDotT() const;


        //- Correct the Lee phaseChange model
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

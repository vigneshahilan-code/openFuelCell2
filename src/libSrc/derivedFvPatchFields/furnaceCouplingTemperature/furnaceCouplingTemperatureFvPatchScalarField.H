/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | openFuelCell2 (sofc-designer patches)
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

Class
    Foam::furnaceCouplingTemperatureFvPatchScalarField

Description
    Stack-to-furnace coupling for the parent-mesh temperature (sofc-designer
    F6): a convective boundary q = h (T_f - Ta), as a mixed condition

        refValue = Ta, refGrad = 0, valueFraction = h/(h + k_p deltaCoeffs)

    so that h (T_f - Ta) = k_p deltaCoeffs (T_c - T_f). k_p is the patch value
    of the conductivity field (the value the energy equation's Laplacian uses),
    looked up in updateCoeffs(), not at construction (the parent T is read
    before k exists). A missing field is fatal.

        terminal_air
        {
            type    furnaceCouplingTemperature;
            h       uniform 84.7;   // W/(m2 K), per patch
            Ta      1023.15;        // K
            kappa   k;              // optional, default k
            log     true;           // optional
            value   uniform 1023.15;
        }

    With log on, each update prints the h read, sum h (T_f - Ta) A and, for
    face 0, T_f, T_c, k, deltaCoeffs and the convective and conductive fluxes
    of the last evaluated face value (computed with the k it was set with).

SourceFiles
    furnaceCouplingTemperatureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef furnaceCouplingTemperatureFvPatchScalarField_H
#define furnaceCouplingTemperatureFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class furnaceCouplingTemperatureFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private data

        //- Heat transfer coefficient [W/(m2 K)]
        scalarField h_;

        //- Furnace temperature [K]
        scalar Ta_;

        //- Name of the conductivity field
        word kappaName_;

        //- Log switch
        Switch log_;

        //- Patch k used for the current valueFraction (for the log)
        scalarField kUsed_;


public:

    //- Runtime type information
    TypeName("furnaceCouplingTemperature");

    using mixedFvPatchScalarField::operator=;


    // Constructors

        furnaceCouplingTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        furnaceCouplingTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        furnaceCouplingTemperatureFvPatchScalarField
        (
            const furnaceCouplingTemperatureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        furnaceCouplingTemperatureFvPatchScalarField
        (
            const furnaceCouplingTemperatureFvPatchScalarField&
        );

        furnaceCouplingTemperatureFvPatchScalarField
        (
            const furnaceCouplingTemperatureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new furnaceCouplingTemperatureFvPatchScalarField(*this)
            );
        }

        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new furnaceCouplingTemperatureFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        virtual void updateCoeffs();

        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //

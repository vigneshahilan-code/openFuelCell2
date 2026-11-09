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
    Foam::totalFluxInletFvPatchScalarField

Description
    Total-flux (Danckwerts) inflow condition (sofc-designer F9) for a scalar
    transported by  div(phi, s) - laplacian(Gamma, s): a mixed condition
        refValue = inletValue, refGrad = 0,
        valueFraction = -F/(-F + Gamma |Sf| deltaCoeffs)   (inflow, F < 0)
    with F the patch value of the flux field `phi` (per face, incl. |Sf|) and
    Gamma the product of the patch values of the fields in `coeffs`. Then the
    outward convective flux F s_f plus the outward diffusive flux
    Gamma |Sf| deltaCoeffs (s_c - s_f) equals F inletValue: only the inlet
    stream enters. Outflow faces (F >= 0) get valueFraction 0 (zeroGradient)
    and are counted in a log line. Before fuelCellSystem construction is
    complete (marker fuelCellSystemConstructed in the parent mesh registry,
    found through Time for region fields) one update per patch is deferred;
    a second is fatal, as is a missing field afterwards.

        fuel_inlet      // parent T
        {
            type        totalFluxInlet;
            phi         rhoCpPhi;
            coeffs      (k);
            inletValue  1023.15;
            refValue    uniform 1023.15;
            refGradient uniform 0;
            valueFraction uniform 1;
            value       uniform 1023.15;
        }

        fuel_inlet      // region species
        {
            type        totalFluxInlet;
            phi         alphaRhoPhi;
            coeffs      (alpha thermo:rho diffSp.H2);
            inletValue  0.3092039219;
            ...
        }

SourceFiles
    totalFluxInletFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef totalFluxInletFvPatchScalarField_H
#define totalFluxInletFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "Switch.H"
#include "wordList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class totalFluxInletFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private data

        //- Name of the convective flux field (surfaceScalarField)
        word phiName_;

        //- Names of the fields whose patch values multiply to Gamma
        wordList coeffNames_;

        //- Inlet value
        scalar inletValue_;

        //- Log switch
        Switch log_;

        //- Deferred updates (at most one, before fuelCellSystem construction)
        label nDeferred_;


    // Private member functions

        //- True once fuelCellSystemConstructed is in the parent registry
        bool systemConstructed() const;


public:

    //- Runtime type information
    TypeName("totalFluxInlet");

    using mixedFvPatchScalarField::operator=;


    // Constructors

        totalFluxInletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        totalFluxInletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        totalFluxInletFvPatchScalarField
        (
            const totalFluxInletFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        totalFluxInletFvPatchScalarField
        (
            const totalFluxInletFvPatchScalarField&
        );

        totalFluxInletFvPatchScalarField
        (
            const totalFluxInletFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new totalFluxInletFvPatchScalarField(*this)
            );
        }

        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new totalFluxInletFvPatchScalarField(*this, iF)
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

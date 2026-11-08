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

\*---------------------------------------------------------------------------*/

#include "furnaceCouplingTemperatureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

#include <sstream>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::furnaceCouplingTemperatureFvPatchScalarField::
furnaceCouplingTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    h_(p.size(), Zero),
    Ta_(0),
    kappaName_("k"),
    log_(false),
    kUsed_(p.size(), Zero)
{
    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = Zero;
}


Foam::furnaceCouplingTemperatureFvPatchScalarField::
furnaceCouplingTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    h_("h", dict, p.size()),
    Ta_(dict.get<scalar>("Ta")),
    kappaName_(dict.getOrDefault<word>("kappa", "k")),
    log_(dict.getOrDefault<Switch>("log", false)),
    kUsed_(p.size(), Zero)
{
    refValue() = Ta_;
    refGrad() = Zero;
    valueFraction() = Zero;

    if (dict.found("value"))
    {
        fvPatchScalarField::operator=(scalarField("value", dict, p.size()));
    }
    else
    {
        fvPatchScalarField::operator=(refValue());
    }

    std::ostringstream os;
    os.precision(12);
    os  << "furnaceCouplingTemperature " << p.name() << ": h read = ["
        << gMin(h_) << " " << gMax(h_) << "] W/(m2 K), Ta = " << Ta_
        << " K, kappa = " << kappaName_;
    Info<< os.str().c_str() << endl;
}


Foam::furnaceCouplingTemperatureFvPatchScalarField::
furnaceCouplingTemperatureFvPatchScalarField
(
    const furnaceCouplingTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    h_(ptf.h_, mapper),
    Ta_(ptf.Ta_),
    kappaName_(ptf.kappaName_),
    log_(ptf.log_),
    kUsed_(ptf.kUsed_, mapper)
{}


Foam::furnaceCouplingTemperatureFvPatchScalarField::
furnaceCouplingTemperatureFvPatchScalarField
(
    const furnaceCouplingTemperatureFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf, ptf.internalField()),
    h_(ptf.h_),
    Ta_(ptf.Ta_),
    kappaName_(ptf.kappaName_),
    log_(ptf.log_),
    kUsed_(ptf.kUsed_)
{}


Foam::furnaceCouplingTemperatureFvPatchScalarField::
furnaceCouplingTemperatureFvPatchScalarField
(
    const furnaceCouplingTemperatureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    h_(ptf.h_),
    Ta_(ptf.Ta_),
    kappaName_(ptf.kappaName_),
    log_(ptf.log_),
    kUsed_(ptf.kUsed_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::furnaceCouplingTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalarField& delta = patch().deltaCoeffs();
    const scalarField Tf(*this);
    const scalarField Tc(patchInternalField());

    // Log the last evaluated face value with the k it was computed with
    if (log_ && gMax(kUsed_) > 0)
    {
        const scalarField qConv(h_*(Tf - Ta_));
        const scalar Q = gSum(qConv*patch().magSf());

        std::ostringstream os;
        os.precision(12);
        os  << "furnaceCouplingTemperature " << patch().name()
            << ": h = " << gMin(h_) << ", sum(h (Tf - Ta) A) = " << Q << " W";
        if (patch().size())
        {
            os  << ", face 0: Tf = " << Tf[0] << ", Tc = " << Tc[0]
                << ", k = " << kUsed_[0] << ", delta = " << delta[0]
                << ", q_conv = " << qConv[0]
                << ", q_cond = " << kUsed_[0]*delta[0]*(Tc[0] - Tf[0]);
        }
        Info<< os.str().c_str() << endl;
    }

    // Patch value of the conductivity: looked up here (fatal if missing)
    const fvPatchScalarField& kp =
        patch().lookupPatchField<volScalarField, scalar>(kappaName_);

    kUsed_ = kp;
    refValue() = Ta_;
    refGrad() = Zero;
    valueFraction() = h_/(h_ + kp*delta);

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::furnaceCouplingTemperatureFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchField<scalar>::write(os);
    h_.writeEntry("h", os);
    os.writeEntry("Ta", Ta_);
    os.writeEntry("kappa", kappaName_);
    os.writeEntry("log", log_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        furnaceCouplingTemperatureFvPatchScalarField
    );
}

// ************************************************************************* //

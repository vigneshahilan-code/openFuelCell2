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

#include "totalFluxInletFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "IOdictionary.H"
#include "polyMesh.H"

#include <sstream>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::totalFluxInletFvPatchScalarField::totalFluxInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    phiName_("phi"),
    coeffNames_(),
    inletValue_(0),
    log_(false),
    nDeferred_(0)
{
    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = Zero;
}


Foam::totalFluxInletFvPatchScalarField::totalFluxInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict),
    phiName_(dict.get<word>("phi")),
    coeffNames_(dict.get<wordList>("coeffs")),
    inletValue_(dict.get<scalar>("inletValue")),
    log_(dict.getOrDefault<Switch>("log", false)),
    nDeferred_(0)
{
    // refValue, refGradient, valueFraction and value come from the dictionary
    // (mixed dictionary constructor) and are held until the first update
    if (coeffNames_.empty())
    {
        FatalIOErrorInFunction(dict)
            << "totalFluxInlet " << p.name() << ": empty coeffs list"
            << exit(FatalIOError);
    }

    std::ostringstream os;
    os.precision(12);
    os  << "totalFluxInlet " << internalField().name() << " " << p.name()
        << ": phi = " << phiName_ << ", coeffs = (";
    forAll(coeffNames_, i)
    {
        os  << (i ? " " : "") << coeffNames_[i];
    }
    os  << "), inletValue = " << inletValue_;
    Info<< os.str().c_str() << endl;
}


Foam::totalFluxInletFvPatchScalarField::totalFluxInletFvPatchScalarField
(
    const totalFluxInletFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    coeffNames_(ptf.coeffNames_),
    inletValue_(ptf.inletValue_),
    log_(ptf.log_),
    nDeferred_(ptf.nDeferred_)
{}


Foam::totalFluxInletFvPatchScalarField::totalFluxInletFvPatchScalarField
(
    const totalFluxInletFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf, ptf.internalField()),
    phiName_(ptf.phiName_),
    coeffNames_(ptf.coeffNames_),
    inletValue_(ptf.inletValue_),
    log_(ptf.log_),
    nDeferred_(ptf.nDeferred_)
{}


Foam::totalFluxInletFvPatchScalarField::totalFluxInletFvPatchScalarField
(
    const totalFluxInletFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    phiName_(ptf.phiName_),
    coeffNames_(ptf.coeffNames_),
    inletValue_(ptf.inletValue_),
    log_(ptf.log_),
    nDeferred_(ptf.nDeferred_)
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::totalFluxInletFvPatchScalarField::systemConstructed() const
{
    // The marker is registered on the parent mesh (fuelCellSystem.C, F6).
    // A parent-mesh field finds it in db(); a region field's db() is the
    // region registry, so the parent mesh is looked up in Time (region0).
    if (db().foundObject<IOdictionary>("fuelCellSystemConstructed"))
    {
        return true;
    }

    const objectRegistry& time = db().time();
    if (!time.foundObject<objectRegistry>(polyMesh::defaultRegion))
    {
        FatalErrorInFunction
            << "totalFluxInlet " << patch().name() << ": parent mesh "
            << polyMesh::defaultRegion << " not found in Time"
            << exit(FatalError);
    }

    return
        time.lookupObject<objectRegistry>(polyMesh::defaultRegion)
       .foundObject<IOdictionary>("fuelCellSystemConstructed");
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::totalFluxInletFvPatchScalarField::updateCoeffs()
{
    // No early return on updated(): every call recomputes from the current
    // flux and coefficient fields (as furnaceCouplingTemperature, F6b)
    const word fieldName = internalField().name();

    if (!systemConstructed())
    {
        ++nDeferred_;
        if (nDeferred_ > 1)
        {
            FatalErrorInFunction
                << "totalFluxInlet " << fieldName << " " << patch().name()
                << ": more than one update before fuelCellSystem construction"
                << " completed" << exit(FatalError);
        }

        Info<< "totalFluxInlet " << fieldName << " " << patch().name()
            << ": deferred (fuelCellSystem not yet constructed), deferral n = "
            << nDeferred_ << endl;
        return;
    }

    const scalarField& F =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    scalarField Gamma(patch().size(), 1.0);
    forAll(coeffNames_, i)
    {
        Gamma *=
            patch().lookupPatchField<volScalarField, scalar>(coeffNames_[i]);
    }

    const scalarField& magSf = patch().magSf();
    const scalarField& delta = patch().deltaCoeffs();

    scalarField& vf = valueFraction();
    label nReversed = 0;
    forAll(vf, facei)
    {
        if (F[facei] < 0)
        {
            vf[facei] =
                -F[facei]
               /(-F[facei] + Gamma[facei]*magSf[facei]*delta[facei]);
        }
        else
        {
            vf[facei] = 0;
            ++nReversed;
        }
    }
    refValue() = inletValue_;
    refGrad() = Zero;

    const label nRev = returnReduce(nReversed, sumOp<label>());
    if (nRev > 0)
    {
        Info<< "totalFluxInlet " << fieldName << " " << patch().name()
            << ": n reversed faces = " << nRev << " (zeroGradient)" << endl;
    }

    if (log_)
    {
        // Net inflow at the last evaluated face value s_f with the current
        // coefficients: outward F s_f + Gamma |Sf| delta (s_c - s_f)
        const scalarField sf(*this);
        const scalarField sc(patchInternalField());
        const scalar net =
            gSum(F*sf + Gamma*magSf*delta*(sc - sf));

        std::ostringstream os;
        os.precision(12);
        os  << "totalFluxInlet " << fieldName << " " << patch().name()
            << ": sum F inletValue = " << gSum(F)*inletValue_
            << ", sum (F s_f + diffusive out) = " << net;
        if (patch().size())
        {
            os  << ", face 0: F = " << F[0] << ", Gamma = " << Gamma[0]
                << ", delta = " << delta[0] << ", valueFraction = " << vf[0];
        }
        Info<< os.str().c_str() << endl;
    }

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::totalFluxInletFvPatchScalarField::write(Ostream& os) const
{
    // All mixed entries first (the dictionary constructor requires them, F6d)
    mixedFvPatchScalarField::write(os);
    os.writeEntry("phi", phiName_);
    os.writeEntry("coeffs", coeffNames_);
    os.writeEntry("inletValue", inletValue_);
    os.writeEntry("log", log_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        totalFluxInletFvPatchScalarField
    );
}

// ************************************************************************* //

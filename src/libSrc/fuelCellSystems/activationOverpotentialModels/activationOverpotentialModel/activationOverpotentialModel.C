/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by the original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "activationOverpotentialModel.H"
#include "phaseSystem.H"
#include "constants.H"

#include <sstream>

// * * * * * * * * * * * * * * * * Static members  * * * * * * * * * * * * * * //

const Foam::scalar F = Foam::constant::physicoChemical::F.value();
const Foam::word Foam::activationOverpotentialModel::phiName("phi");
const Foam::word Foam::activationOverpotentialModel::jName("J");

namespace Foam
{
    defineTypeNameAndDebug(activationOverpotentialModel, 0);
    defineRunTimeSelectionTable(activationOverpotentialModel, dictionary);
}

// * * * * * * * * * * * * * * * * Protected function  * * * * * * * * * * * * * * //

void Foam::activationOverpotentialModel::createNernst()
{
    //- Create nernst model
    if (!nernst_.valid())
    {
        const phaseSystem& phaseSys = phase_.mesh().template
            lookupObject<phaseSystem>(phaseSystem::propertiesName);

        if (phaseChange_)
        {
            //- Get the name of the other phase
            const word name = Pair<word>
            (
                phaseSys.phases()[0].name(),
                phaseSys.phases()[1].name()
            ).other(phase_.name());

            const phaseModel& phase1 = phaseSys.phases()[name];
            nernst_ = nernstModel::New
            (
                phase_, phase1, dict_.subDict("nernstModel")
            );
        }
        else
        {
            nernst_ = nernstModel::New
            (
                phase_, phase_, dict_.subDict("nernstModel")
            );
        }
    }
}

void Foam::activationOverpotentialModel::correctCouplingSources()
{
    if (!nernst_.valid())
    {
        return;
    }

    const regionType& fluidPhase =
        region(word(regions_.subDict("fluid").lookup("name")));
    const regionType& electronPhase =
        region(word(regions_.subDict("electron").lookup("name")));
    const regionType& ionPhase =
        region(word(regions_.subDict("ion").lookup("name")));

    const scalarField& SE =
        electronPhase.lookupObject<volScalarField>(jName);
    const scalarField& SI =
        ionPhase.lookupObject<volScalarField>(jName);
    const scalarField& jPrimeE =
        electronPhase.lookupObject<volScalarField>("jPrime");
    const scalarField& jPrimeI =
        ionPhase.lookupObject<volScalarField>("jPrime");
    scalarField& jSourceE = const_cast<volScalarField&>
    (
        electronPhase.lookupObject<volScalarField>("jSource")
    );
    scalarField& jSourceI = const_cast<volScalarField&>
    (
        ionPhase.lookupObject<volScalarField>("jSource")
    );
    const scalarField& phiE =
        electronPhase.lookupObject<volScalarField>(phiName);
    const scalarField& phiI =
        ionPhase.lookupObject<volScalarField>(phiName);

    const scalarField& eta = eta_;
    const scalarField& nernst = nernst_()();

    label znId = fluidPhase.cellZones().findZoneID(zoneName_);
    const labelList& cells = fluidPhase.cellZones()[znId];

    forAll(cells, cellI)
    {
        const label fluidId = cells[cellI];
        const label electronId =
            electronPhase.cellMap()[fluidPhase.cellMapIO()[fluidId]];
        const label ionId =
            ionPhase.cellMap()[fluidPhase.cellMapIO()[fluidId]];

        jSourceE[electronId] =
            SE[electronId]
          + jPrimeE[electronId]*(phiI[ionId] + nernst[fluidId] + eta[fluidId]);
        jSourceI[ionId] =
            SI[ionId]
          + jPrimeI[ionId]*(phiE[electronId] - nernst[fluidId] - eta[fluidId]);
    }
}


void Foam::activationOverpotentialModel::checkCouplingSources()
{
    if (!checkCouplingSources_ || !nernst_.valid())
    {
        return;
    }

    const regionType& electronPhase =
        region(word(regions_.subDict("electron").lookup("name")));
    const regionType& ionPhase =
        region(word(regions_.subDict("ion").lookup("name")));

    scalarField& jSourceE = const_cast<volScalarField&>
    (
        electronPhase.lookupObject<volScalarField>("jSource")
    );
    scalarField& jSourceI = const_cast<volScalarField&>
    (
        ionPhase.lookupObject<volScalarField>("jSource")
    );

    //- jSource as written by correct(), at unchanged phi
    const scalarField jE0(jSourceE);
    const scalarField jI0(jSourceI);

    correctCouplingSources();

    const scalar dE = gMax(mag(jSourceE - jE0)());
    const scalar dI = gMax(mag(jSourceI - jI0)());
    const label nE = returnReduce(label(Foam::sum(pos(mag(jSourceE - jE0) - VSMALL))), sumOp<label>());
    const label nI = returnReduce(label(Foam::sum(pos(mag(jSourceI - jI0) - VSMALL))), sumOp<label>());

    //- leave the run exactly as correct() left it (the check is log only)
    jSourceE = jE0;
    jSourceI = jI0;

    Info<< "P5 source-form identity " << zoneName_
        << ": max |d jSourceE| = " << dE << ", max |d jSourceI| = " << dI
        << ", differing cells E/I = " << nE << "/" << nI << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::activationOverpotentialModel::activationOverpotentialModel
(
    const phaseModel& phase,
    const dictionary& dict
)
:
    phase_(phase),
    dict_(dict),
    eta_
    (
        IOobject
        (
            IOobject::groupName("eta", phase.mesh().name()),
            phase.mesh().time().timeName(),
            phase.mesh(),
            //- F13 (FP-1): readEta true starts from a written relaxed eta
            //  (restarts); default false = upstream (eta starts at zero)
            dict.getOrDefault<Switch>("readEta", false)
          ? IOobject::READ_IF_PRESENT
          : IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phase.mesh(),
        dimensionedScalar
        (
            "eta",
            dimensionSet(1, 2, -3, 0, 0, -1, 0),
            0.0
        )
    ),
    phaseChange_(dict_.lookupOrDefault<Switch>("phaseChange", true)),
    zoneName_(dict_.lookup("zoneName")),
    regions_(dict_.subDict("regions")),
    species_(dict_.subDict("species")),
    phiNames_(dict_.lookup("phiNames")),
    relax_(dict_.lookupOrDefault<scalar>("relax", 1.0)),
    alpha_(dict_.lookupOrDefault<scalar>("alpha", 0.5)),
    gamma_(dict_.lookupOrDefault<scalar>("gamma", 1.0)),
    j0_("j0", dimCurrent/dimVolume, dict_),

    j_
    (
        IOobject
        (
            "j",
            phase_.mesh().time().timeName(),
            phase_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phase_.mesh(),
        dimensionedScalar
        (
            "j",
            dimCurrent/dimVolume,
            0.0
        )
    ),
    includeDj0dT_(dict_.lookupOrDefault<Switch>("includeDj0dT", false)),
    checkCell_(dict_.getOrDefault<label>("checkDQdotdT", -1))
{
    if (includeDj0dT_)
    {
        Info<< "Activation model at " << zoneName_
            << ": includeDj0dT on, dQdotdT += (d ln j0/dT) j"
            << " (sign eta - T dS/(|n| F)); zero for a constant j0" << endl;
    }
    if (dict_.getOrDefault<Switch>("readEta", false))
    {
        Info<< "Activation model at " << zoneName_ << ": readEta on, "
            << eta_.name() << " read from time " << eta_.instance() << ": "
            << (eta_.template typeHeaderOk<volScalarField>(true) ? "yes" : "no (absent, zero start)")
            << ", min/max " << gMin(eta_.primitiveField()) << " / "
            << gMax(eta_.primitiveField()) << endl;
    }
    checkCouplingSources_ = dict_.getOrDefault<Switch>("checkCouplingSources", false);
    if (checkCouplingSources_)
    {
        Info<< "Activation model at " << zoneName_
            << ": checkCouplingSources on (P5 source-form identity; log only)"
            << endl;
    }
    if (checkCell_ >= 0)
    {
        Info<< "Activation model at " << zoneName_
            << ": checkDQdotdT at cell " << checkCell_
            << ", central difference of j0(T) with h = " << checkStep << " K"
            << " (eta, dS and composition frozen); log only" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::activationOverpotentialModel::~activationOverpotentialModel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::tmp<Foam::volScalarField>
Foam::activationOverpotentialModel::Qdot() const
{
    tmp<volScalarField> q
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("q", phase_.name()),
                eta_.mesh().time().timeName(),
                eta_.mesh()
            ),
            eta_.mesh(),
            dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0)
        )
    );

    if (!nernst_.valid())
    {
        return q;
    }

    //- get the sub models
    const regionType& fluidPhase = region
    (
        word(regions_.subDict("fluid").lookup("name"))
    );

    const scalarField& j = j_;
    const scalarField& T = phase_.thermo().T();

    const scalarField& eta = eta_;

    //- only consider catalyst zone
    label znId = fluidPhase.cellZones().findZoneID(zoneName_);
    const labelList& cells = fluidPhase.cellZones()[znId];

    scalar n = nernst_->rxnList()["e"];
    scalar sign = n/mag(n);

    forAll(cells, cellI)
    {
        label fluidId = cells[cellI];

        q.ref()[fluidId] =
            j[fluidId]*(sign*eta[fluidId] - T[fluidId]*nernst_->deltaS()[fluidId]
            /mag(nernst_->rxnList()["e"])/F);
    }

    if (nernst_->mixingEntropy())
    {
        //- F10 log (only with mixingEntropy on): per zone sums [W]
        const scalarField& V = fluidPhase.V();
        scalar qRev = 0, qMix = 0, qAct = 0;
        forAll(cells, cellI)
        {
            const label c = cells[cellI];
            qRev -= j[c]*T[c]*nernst_->deltaS()[c]/mag(n)/F*V[c];
            qMix -= j[c]*T[c]*nernst_->deltaSmix()[c]/mag(n)/F*V[c];
            qAct += j[c]*sign*eta[c]*V[c];
        }
        reduce(qRev, sumOp<scalar>());
        reduce(qMix, sumOp<scalar>());
        reduce(qAct, sumOp<scalar>());

        std::ostringstream os;
        os.precision(12);
        os  << "Qdot " << zoneName_ << ": sum reversible heat = " << qRev
            << " W, composition part = " << qMix
            << " W, sum activation heat = " << qAct << " W";
        Info<< os.str().c_str() << endl;
    }

    return q;
}


Foam::tmp<Foam::volScalarField>
Foam::activationOverpotentialModel::dQdotdT() const
{
    tmp<volScalarField> q
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("dQdotdT", phase_.name()),
                eta_.mesh().time().timeName(),
                eta_.mesh()
            ),
            eta_.mesh(),
            dimensionedScalar
            (
                "zero",
                dimEnergy/dimVolume/dimTime/dimTemperature,
                0
            )
        )
    );

    if (!nernst_.valid())
    {
        return q;
    }

    const regionType& fluidPhase = region
    (
        word(regions_.subDict("fluid").lookup("name"))
    );

    const scalarField& j = j_;

    label znId = fluidPhase.cellZones().findZoneID(zoneName_);
    const labelList& cells = fluidPhase.cellZones()[znId];

    const scalar n = nernst_->rxnList()["e"];

    forAll(cells, cellI)
    {
        label fluidId = cells[cellI];

        q.ref()[fluidId] =
          - j[fluidId]*nernst_->deltaS()[fluidId]/mag(n)/F;
    }

    if (includeDj0dT_)
    {
        //- d(Qdot)/dT through j0 at fixed eta:
        //  Qdot = j (sign eta - T dS/(|n| F)), dj/dT = (d ln j0/dT) j
        const scalarField& T = phase_.thermo().T();
        const scalarField& eta = eta_;
        const scalar sign = n/mag(n);

        forAll(cells, cellI)
        {
            label fluidId = cells[cellI];

            q.ref()[fluidId] +=
                dlnj0dT(T[fluidId])*j[fluidId]
               *(
                    sign*eta[fluidId]
                  - T[fluidId]*nernst_->deltaS()[fluidId]/mag(n)/F
                );
        }
    }

    if (checkCell_ >= 0)
    {
        //- Finite-difference check of the includeDj0dT term (log only):
        //  analytic (d ln j0/dT) j B against
        //  [j0(T + h) - j0(T - h)]/(2 h j0(T)) j B, with B, eta, dS and the
        //  composition frozen at T (j is proportional to j0 at fixed eta).
        if (cells.find(checkCell_) < 0)
        {
            FatalErrorInFunction
                << "checkDQdotdT cell " << checkCell_
                << " is not in cellZone " << zoneName_
                << exit(FatalError);
        }

        const label c = checkCell_;
        const scalar Tc = phase_.thermo().T()[c];
        const scalar sign = n/mag(n);
        const scalar B =
            sign*eta_[c] - Tc*nernst_->deltaS()[c]/mag(n)/F;
        const scalar j0 = j0At(Tc);

        std::ostringstream os;
        os.precision(12);
        os  << "checkDQdotdT " << zoneName_ << " cell " << c
            << ": T = " << Tc << ", j = " << j[c] << ", B = " << B
            << ", includeDj0dT " << (includeDj0dT_ ? "on" : "off");

        if (j0 > 0)
        {
            const scalar h = checkStep;
            const scalar analytic = dlnj0dT(Tc)*j[c]*B;
            const scalar fd =
                (j0At(Tc + h) - j0At(Tc - h))/(2*h*j0)*j[c]*B;

            os  << ", analytic = " << analytic
                << ", fd(h = " << h << " K) = " << fd
                << ", rel = "
                << (mag(fd) > VSMALL ? (analytic - fd)/fd : scalar(0));
        }
        else
        {
            os  << ", no Arrhenius j0: nothing to compare";
        }

        Info<< os.str().c_str() << endl;
    }

    return q;
}

// ************************************************************************* //

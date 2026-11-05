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

#include "ButlerVolmer.H"
#include "phaseModel.H"


// Default bound on the magnitude of the Boltzmann exponent argument.
// exp(709.78) is the double-precision overflow limit; 50 leaves a wide
// margin while still being far outside the physically meaningful range
// (|n alpha F eta / R T| = 50 corresponds to |eta| ~ 1.3 V at 330 K).
static const Foam::scalar butlerVolmerDefaultMaxExpArg = 50;

// Default bound on the magnitude of the activation overpotential [V]. A real
// cell operates at activation overpotentials well under 1 V; 5 V is an order
// of magnitude above any converged solution, so the bound is inactive for
// physical states and only arrests an unconverged potential field.
static const Foam::scalar butlerVolmerDefaultMaxEta = 5;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::activationOverpotentialModels::ButlerVolmer<Thermo>::ButlerVolmer
(
    const phaseModel& phase,
    const dictionary& dict
)
:
    ActivationOverpotentialModel<Thermo>(phase, dict),
    maxExpArg_
    (
        dict.lookupOrDefault<scalar>("maxExpArg", butlerVolmerDefaultMaxExpArg)
    ),
    maxEta_
    (
        dict.lookupOrDefault<scalar>("maxEta", butlerVolmerDefaultMaxEta)
    ),
    moleFractionBasis_(false),
    arrhenius_(dict.found("A") || dict.found("Ea")),
    A_(arrhenius_ ? dict.get<scalar>("A") : 0),
    Ea_(arrhenius_ ? dict.get<scalar>("Ea") : 0)
{
    if (maxExpArg_ <= 0)
    {
        FatalErrorInFunction
            << "maxExpArg must be positive, got " << maxExpArg_
            << exit(FatalError);
    }

    if (maxEta_ <= 0)
    {
        FatalErrorInFunction
            << "maxEta must be positive, got " << maxEta_
            << exit(FatalError);
    }

    const word basis
    (
        dict.lookupOrDefault<word>("activityBasis", "concentration")
    );
    if (basis == "moleFraction")
    {
        moleFractionBasis_ = true;
        Info<< "ButlerVolmer at " << this->zoneName_
            << ": activityBasis moleFraction (a = X^ksi, cRef not read)"
            << endl;
    }
    else if (basis != "concentration")
    {
        FatalIOErrorInFunction(dict)
            << "activityBasis must be concentration or moleFraction, got "
            << basis << exit(FatalIOError);
    }

    if (arrhenius_)
    {
        if (A_ <= 0)
        {
            FatalIOErrorInFunction(dict)
                << "A must be positive, got " << A_
                << exit(FatalIOError);
        }

        Info<< "ButlerVolmer at " << this->zoneName_
            << ": Arrhenius j0 = A*T*exp(-Ea/(R*T)), A = " << A_
            << " A/(m3 K), Ea = " << Ea_ << " J/mol, R = " << Rgas
            << " J/(mol K); j0 entry ignored" << endl;
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Thermo>
Foam::activationOverpotentialModels::ButlerVolmer<Thermo>::~ButlerVolmer()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class Thermo>
void Foam::activationOverpotentialModels::ButlerVolmer<Thermo>::correct()
{
    //- Create Nernst model
    this->createNernst();

    //- Correct
    this->nernst_->correct();

    //- Get sub regions
    //- Refer to regionType
    //- Including: fluid, electron (BPP + GDL + CL), ion (CLs + membrane)
    const regionType& fluidPhase = this->region
    (
        word(this->regions_.subDict("fluid").lookup("name"))
    );
    const regionType& electronPhase = this->region
    (
        word(this->regions_.subDict("electron").lookup("name"))
    );
    const regionType& ionPhase = this->region
    (
        word(this->regions_.subDict("ion").lookup("name"))
    );

    //- Source/sink terms for electron and ion fields
    //- See names in regions/electronIon/electronIon.C
    scalarField& SE = const_cast<volScalarField&>
    (
        electronPhase.template lookupObject<volScalarField>("J")
    );
    scalarField& SI = const_cast<volScalarField&>
    (
        ionPhase.template lookupObject<volScalarField>("J")
    );
    scalarField& jPrimeE = const_cast<volScalarField&>
    (
        electronPhase.template lookupObject<volScalarField>("jPrime")
    );
    scalarField& jPrimeI = const_cast<volScalarField&>
    (
        ionPhase.template lookupObject<volScalarField>("jPrime")
    );
    scalarField& jSourceE = const_cast<volScalarField&>
    (
        electronPhase.template lookupObject<volScalarField>("jSource")
    );
    scalarField& jSourceI = const_cast<volScalarField&>
    (
        ionPhase.template lookupObject<volScalarField>("jSource")
    );

    //- Potential fields
    const scalarField& phiE = electronPhase.template
        lookupObject<volScalarField>
        (
            word(this->phiNames_["electron"])
        );
    const scalarField& phiI = ionPhase.template
        lookupObject<volScalarField>
        (
            word(this->phiNames_["ion"])
        );

    //- Reference
    scalarField& eta = this->eta_;
    //- Nernst field
    scalarField& nernst = this->nernst_()();
    //- Current density
    scalarField& j = this->j_;

    //- Access temperature and reactant
    const scalarField& T = this->thermo_.T();
    const scalarField& s = this->phase_;

    //- Mixture mole fraction
    const scalarField W(this->thermo_.W()/1000.0);
    //- Mixture density
    const scalarField& rho= this->phase_.thermo().rho();
    //- Store values for all species

    scalarField coeff(this->phase_.mesh().nCells(), 1.0);
    //- Loop for all relative species
    forAllConstIter(dictionary, this->species_, iter)
    {
        if (iter().isDict())
        {
            const word& name = iter().keyword();
            const dictionary& dict0 = iter().dict();

            scalar ksi = readScalar(dict0.lookup("ksi"));

            const scalarField& X = this->phase_.X(name);

            if (moleFractionBasis_)
            {
                coeff *= Foam::pow(X, ksi);
            }
            else
            {
                scalar cRef = readScalar(dict0.lookup("cRef"));
                coeff *= Foam::pow(X*rho/W/cRef, ksi);
            }
        }
    }

    //- Only consider catalyst zone
    label znId = fluidPhase.cellZones().findZoneID(this->zoneName_);
    const labelList& cells = fluidPhase.cellZones()[znId];

    //- electron transfer
    //- cathode side, < 0
    //- anode side, >0
    scalar n = this->nernst_->rxnList()["e"];
    scalar sign = n/mag(n);

    //- The total current: volume averaged
    scalar Rj(0.0);

    //- Number of catalyst cells whose Boltzmann exponent had to be clipped in
    //  this call, and the largest unclipped magnitude seen. A non-zero count
    //  means the electrochemical coupling is being driven outside its
    //  converged range; it is reported so the clamp is never silent.
    //- Arrhenius j0: range over the catalyst cells, logged per call
    const scalar j0Const = this->j0_.value();
    scalar j0Min(GREAT);
    scalar j0Max(-GREAT);
    scalar TMin(GREAT);
    scalar TMax(-GREAT);

    label nClippedLocal = 0;
    label nClipped = 0;
    scalar maxExpArgSeen = 0;

    forAll(cells, cellI)
    {
        //- get cell IDs
        label fluidId = cells[cellI];
        label electronId = electronPhase.cellMap()[fluidPhase.cellMapIO()[fluidId]];
        label ionId = ionPhase.cellMap()[fluidPhase.cellMapIO()[fluidId]];

        //- activation overpotential
        eta[fluidId] =
        (
            eta[fluidId]*(scalar(1) - this->relax_)
          + (phiE[electronId] - phiI[ionId] - nernst[fluidId])
          * this->relax_
        );

        // Bound the overpotential itself, not only the Boltzmann exponent.
        // The activation heat source posted to the energy equation is j*eta,
        // so a finite clamped j combined with an unconverged eta of order 1e42
        // V still overflows the source and destroys the temperature field.
        // Clamping eta here keeps the heat source bounded and is inactive for
        // any converged solution.
        if (mag(eta[fluidId]) > maxEta_)
        {
            eta[fluidId] = Foam::sign(eta[fluidId])*maxEta_;
        }

        //- Exchange current: constant, or Arrhenius in the local T.
        //- j, jPrime and jSource below use this same j0.
        scalar j0 = j0Const;
        if (arrhenius_)
        {
            j0 = A_*T[fluidId]*Foam::exp(-Ea_/(Rgas*T[fluidId]));
            j0Min = Foam::min(j0Min, j0);
            j0Max = Foam::max(j0Max, j0);
            TMin = Foam::min(TMin, T[fluidId]);
            TMax = Foam::max(TMax, T[fluidId]);
        }

        const scalar prefactor =
            j0
           *coeff[fluidId]
           *Foam::pow(s[fluidId], this->gamma_);
        // Boltzmann exponents, bounded by maxExpArg_. The argument is
        // proportional to eta/T, so during an unconverged electrochemical
        // coupling iteration a transiently large overpotential drives it past
        // the double-precision overflow limit (~709): Foam::exp returns inf,
        // the inf propagates into the exchange-current sources J, jPrime and
        // jSource of the electron and ion regions, and the next potential
        // solve traps a hardware floating point exception. Clamping the
        // argument keeps j and its Jacobian finite, which lets the coupling
        // iteration recover instead of aborting. The Jacobian below is built
        // from the clamped exponentials, so the linearisation stays
        // consistent with the current that is actually reported.
        const scalar RT = Rgas*T[fluidId];

        scalar argForward = n*this->alpha_*F*eta[fluidId]/RT;
        scalar argReverse = -n*(scalar(1) - this->alpha_)*F*eta[fluidId]/RT;

        const scalar argMag = max(mag(argForward), mag(argReverse));
        maxExpArgSeen = max(maxExpArgSeen, argMag);
        if (argMag > maxExpArg_)
        {
            ++nClippedLocal;
        }

        argForward = min(max(argForward, -maxExpArg_), maxExpArg_);
        argReverse = min(max(argReverse, -maxExpArg_), maxExpArg_);

        const scalar expForward = Foam::exp(argForward);
        const scalar expReverse = Foam::exp(argReverse);

        // Signed Butler-Volmer current and its positive equation Jacobian.
        j[fluidId] = prefactor*(expForward - expReverse);
        const scalar jPrime =
            mag(n)*F/RT*prefactor
           *(
                this->alpha_*expForward
              + (scalar(1) - this->alpha_)*expReverse
            );

        //- anode side: SE = Rj, SI = -Rj
        //- cathode side: SE = -Rj, SI = Rj
        SE[electronId] = -sign*j[fluidId];
        SI[ionId] = -SE[electronId];
        jPrimeE[electronId] = jPrime;
        jPrimeI[ionId] = jPrime;
        jSourceE[electronId] =
            SE[electronId]
          + jPrime*(phiI[ionId] + nernst[fluidId] + eta[fluidId]);
        jSourceI[ionId] =
            SI[ionId]
          + jPrime*(phiE[electronId] - nernst[fluidId] - eta[fluidId]);

        Rj += fluidPhase.V()[fluidId] * j[fluidId];
    }

    reduce(Rj, sumOp<scalar>());

    // Report clipping rather than letting it pass silently: a saturated
    // Boltzmann factor is no longer the Butler-Volmer law, so the reported
    // current is only a numerical holding action and the coupling has not
    // converged. The warning identifies the zone and the exponent magnitude
    // that had to be bounded.
    // Reduced unconditionally: a conditional reduce would be a collective
    // called by only some ranks and would deadlock under MPI.
    nClipped = nClippedLocal;
    reduce(nClipped, sumOp<label>());
    reduce(maxExpArgSeen, maxOp<scalar>());

    if (nClipped)
    {
        WarningInFunction
            << "zone " << this->zoneName_ << ": bounded the Boltzmann exponent"
            << " in " << nClipped << " catalyst cells"
            << " (max |n alpha F eta / R T| = " << maxExpArgSeen
            << " > maxExpArg = " << maxExpArg_ << ')'
            << ". The current density is saturated and the electrochemical"
            << " coupling has not converged." << endl;
    }

    if (arrhenius_)
    {
        reduce(j0Min, minOp<scalar>());
        reduce(j0Max, maxOp<scalar>());
        reduce(TMin, minOp<scalar>());
        reduce(TMax, maxOp<scalar>());

        const int oldPrecision = Sout.precision(10);
        Info<< "Arrhenius j0 at " << this->zoneName_
            << ": min = " << j0Min << ", max = " << j0Max
            << " (T min = " << TMin << ", max = " << TMax << ")" << endl;
        Sout.precision(oldPrecision);
    }

    Info << "Total current (A) at " << this->zoneName_
         << ": " << Rj << endl;
}

// ************************************************************************* //

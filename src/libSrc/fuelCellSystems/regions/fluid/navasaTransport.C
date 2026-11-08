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

#include "navasaTransport.H"
#include "volFields.H"

#include <sstream>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::navasaTransport::navasaTransport
(
    const word& region,
    const dictionary& dict,
    const basicSpecieMixture& composition
)
:
    region_(region),
    muCoeffs_(composition.species().size()),
    kappaCoeffs_(composition.species().size()),
    logCell_(dict.getOrDefault<label>("logCell", -1)),
    logPatch_(dict.getOrDefault<word>("logPatch", word::null))
{
    const dictionary& spDict = dict.subDict("species");

    forAll(composition.species(), i)
    {
        const word& sp = composition.species()[i];

        if (!spDict.found(sp))
        {
            FatalIOErrorInFunction(dict)
                << "navasaTransport in region " << region_
                << ": species " << sp << " has no coefficients"
                << exit(FatalIOError);
        }

        const dictionary& c = spDict.subDict(sp);
        muCoeffs_[i] = c.get<FixedList<scalar, 7>>("mu");
        kappaCoeffs_[i] = c.get<FixedList<scalar, 7>>("kappa");
    }

    Info<< "Fluid region " << region_ << ": navasaTransport on"
        << " (Navasa A1/A3 per species, Wilke mu, Wassiljewa-Mason-Saxena"
        << " z = 1 kappa; Cp janaf)" << endl;
}


// * * * * * * * * * * * * * * * Private functions  * * * * * * * * * * * * //

void Foam::navasaTransport::moleFractions
(
    const scalarList& Y,
    const scalarList& W,
    scalarList& X
)
{
    scalar sum = 0;
    forAll(Y, i)
    {
        X[i] = Y[i]/W[i];
        sum += X[i];
    }
    forAll(X, i)
    {
        X[i] /= sum;
    }
}


void Foam::navasaTransport::mix
(
    const scalarList& X,
    const scalarList& W,
    const scalar T,
    scalar& mu,
    scalar& kappa
) const
{
    if (T < 273.15 || T > 1473.15)
    {
        FatalErrorInFunction
            << "navasaTransport in region " << region_ << ": T = " << T
            << " K outside the Navasa range 0-1200 C"
            << exit(FatalError);
    }

    const label n = X.size();
    scalarList mui(n), ki(n);
    const scalar t = T/1000.0;

    for (label i = 0; i < n; ++i)
    {
        scalar m = 0, k = 0, tp = 1;
        for (label p = 0; p < 7; ++p)
        {
            m += muCoeffs_[i][p]*tp;
            k += kappaCoeffs_[i][p]*tp;
            tp *= t;
        }
        mui[i] = 1e-7*m;
        ki[i] = 0.01*k;
    }

    mu = 0;
    kappa = 0;
    for (label i = 0; i < n; ++i)
    {
        scalar den = 0;
        for (label j = 0; j < n; ++j)
        {
            const scalar phi =
                sqr(1 + Foam::sqrt(mui[i]/mui[j])*Foam::pow(W[j]/W[i], 0.25))
               /Foam::sqrt(8*(1 + W[i]/W[j]));
            den += X[j]*phi;
        }
        mu += X[i]*mui[i]/den;
        kappa += X[i]*ki[i]/den;
    }
}


void Foam::navasaTransport::logLine
(
    const word& where,
    const basicSpecieMixture& comp,
    const scalarList& X,
    const scalar T,
    const scalar mu,
    const scalar kappa,
    const scalar Cp
) const
{
    std::ostringstream os;
    os.precision(12);
    os  << region_ << ": navasaTransport " << where << ": T = " << T;
    forAll(X, i)
    {
        os << ", X " << comp.species()[i] << " = " << X[i];
    }
    os  << ", mu = " << mu << ", kappa = " << kappa << ", Cp = " << Cp
        << ", alpha = " << kappa/Cp;
    Info<< os.str().c_str() << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::navasaTransport::apply(rhoReactionThermo& thermo) const
{
    const basicSpecieMixture& comp = thermo.composition();
    const fvMesh& mesh = thermo.T().mesh();
    const label n = comp.species().size();

    volScalarField& mu =
        mesh.lookupObjectRef<volScalarField>(IOobject::scopedName("thermo", "mu"));
    volScalarField& alpha =
        mesh.lookupObjectRef<volScalarField>(IOobject::scopedName("thermo", "alpha"));

    const volScalarField& T = thermo.T();
    const tmp<volScalarField> tCp = thermo.Cp();
    const volScalarField& Cp = tCp();

    scalarList W(n), Y(n), X(n);
    forAll(W, i)
    {
        W[i] = comp.W(i);
    }

    scalar muMin = GREAT, muMax = -GREAT, kMin = GREAT, kMax = -GREAT;

    // Cells
    forAll(T, celli)
    {
        forAll(Y, i)
        {
            Y[i] = comp.Y(i)[celli];
        }
        moleFractions(Y, W, X);

        scalar m, k;
        mix(X, W, T[celli], m, k);
        mu[celli] = m;
        alpha[celli] = k/Cp[celli];

        muMin = Foam::min(muMin, m);
        muMax = Foam::max(muMax, m);
        kMin = Foam::min(kMin, k);
        kMax = Foam::max(kMax, k);

        if (celli == logCell_)
        {
            logLine("cell " + Foam::name(celli), comp, X, T[celli], m, k, Cp[celli]);
        }
    }

    // Boundary faces
    forAll(T.boundaryField(), patchi)
    {
        const fvPatchScalarField& Tp = T.boundaryField()[patchi];
        const fvPatchScalarField& Cpp = Cp.boundaryField()[patchi];
        fvPatchScalarField& mup = mu.boundaryFieldRef()[patchi];
        fvPatchScalarField& alphap = alpha.boundaryFieldRef()[patchi];
        const word& pName = mesh.boundary()[patchi].name();

        forAll(Tp, facei)
        {
            forAll(Y, i)
            {
                Y[i] = comp.Y(i).boundaryField()[patchi][facei];
            }
            moleFractions(Y, W, X);

            scalar m, k;
            mix(X, W, Tp[facei], m, k);
            mup[facei] = m;
            alphap[facei] = k/Cpp[facei];

            if (facei == 0 && pName == logPatch_)
            {
                logLine("patch " + pName + " face 0", comp, X, Tp[facei], m, k, Cpp[facei]);
            }
        }
    }

    reduce(muMin, minOp<scalar>());
    reduce(muMax, maxOp<scalar>());
    reduce(kMin, minOp<scalar>());
    reduce(kMax, maxOp<scalar>());

    std::ostringstream os;
    os.precision(10);
    os  << region_ << ": navasaTransport mu [" << muMin << " " << muMax
        << "], kappa [" << kMin << " " << kMax << "]";
    Info<< os.str().c_str() << endl;
}


// ************************************************************************* //

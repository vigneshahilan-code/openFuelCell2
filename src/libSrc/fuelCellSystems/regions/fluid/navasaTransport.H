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
    Foam::navasaTransport

Description
    Optional per-species gas transport for a fluid region (sofc-designer F7):
    Navasa (2018) Tables A1 / A3 polynomials per species,

        mu_i    = 1e-7 sum_k a_k (T/1000)^k    [Pa s]   (micro-poise x 1e-7)
        kappa_i = 0.01 sum_k c_k (T/1000)^k    [W/(m K)]

    mixed with Wilke (mu) and Wassiljewa-Mason-Saxena with z = 1 (kappa),
    with the same Phi_ij:

        Phi_ij = (1 + (mu_i/mu_j)^0.5 (W_j/W_i)^0.25)^2 / sqrt(8 (1 + W_i/W_j))
        mix    = sum_i X_i f_i / sum_j X_j Phi_ij

    The thermo's mixture fields thermo:mu and thermo:alpha (= kappa/Cp, Cp
    from the janaf mixture, unchanged) are overwritten, cells and boundary
    faces. Selected by the region dictionary:

        navasaTransport
        {
            active      true;
            species { H2 { mu (7 coeffs); kappa (7 coeffs); } ... }
            logCell     <label>;    // optional
            logPatch    <name>;     // optional, face 0
        }

    Absent or inactive: nothing is constructed, behaviour unchanged.

SourceFiles
    navasaTransport.C

\*---------------------------------------------------------------------------*/

#ifndef navasaTransport_H
#define navasaTransport_H

#include "rhoReactionThermo.H"
#include "dictionary.H"
#include "FixedList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class navasaTransport
{
    // Private data

        //- Region name (for the log)
        word region_;

        //- Coefficients per thermo species index
        List<FixedList<scalar, 7>> muCoeffs_;
        List<FixedList<scalar, 7>> kappaCoeffs_;

        //- Optional log cell (-1: none) and log patch ("" : none)
        label logCell_;
        word logPatch_;


    // Private member functions

        //- Mixture mu and kappa at T for mole fractions X
        void mix
        (
            const scalarList& X,
            const scalarList& W,
            const scalar T,
            scalar& mu,
            scalar& kappa
        ) const;

        //- Mole fractions from mass fractions
        static void moleFractions
        (
            const scalarList& Y,
            const scalarList& W,
            scalarList& X
        );

        //- One log line
        void logLine
        (
            const word& where,
            const basicSpecieMixture& comp,
            const scalarList& X,
            const scalar T,
            const scalar mu,
            const scalar kappa,
            const scalar Cp
        ) const;


public:

    // Constructors

        //- Construct from the region name, the navasaTransport dictionary
        //  and the region's species (every species must be listed)
        navasaTransport
        (
            const word& region,
            const dictionary& dict,
            const basicSpecieMixture& composition
        );


    // Member functions

        //- Overwrite thermo:mu and thermo:alpha (cells and boundary faces)
        void apply(rhoReactionThermo& thermo) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //

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
#include "fvCFD.H"
#include "fuelCellSystem.H"
#include "phaseSystem.H"
#include "HashPtrTable.H"
#include "patchToPatchInterpolation.H"
#include "regionTypeList.H"


// * * * * * * * * * * * * * * * Private functions * * * * * * * * * * * * * //

Foam::IOobject Foam::fuelCellSystem::createIOobject
(
    const fvMesh& mesh
) const
{
    IOobject io
    (
        "cellProperties",
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (io.typeHeaderOk<IOdictionary>(true))
    {
        Info<< "Get cell properties from " << io.name() << nl << endl;

        io.readOpt(IOobject::MUST_READ_IF_MODIFIED);
    }
    else
    {
        Info<< "No cell properties are provided... " << nl << endl;

        io.readOpt(IOobject::NO_READ);
    }

    return io;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fuelCellSystem::fuelCellSystem
(
    const fvMesh& mesh
)
:
    IOdictionary(createIOobject(mesh)),

    mesh_(mesh),

    T_
    (
        IOobject
        (
            "T",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    Qdot_
    (
        T_,
        dimEnergy/dimTime
    ),

    dQdotdT_
    (
        IOobject
        (
            "dQdotdT",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimEnergy/dimVolume/dimTime/dimTemperature, 0)
    ),

    QdotJoule_
    (
        IOobject
        (
            "QdotJoule",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0)
    ),

    k_
    (
        IOobject
        (
            "k",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    contErrCp_
    (
        IOobject
        (
            "contErrCp",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            dimDensity/dimTime*dimSpecificHeatCapacity,
            Zero
        )
    ),

    Cp_
    (
        IOobject
        (
            "Cp",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            dimSpecificHeatCapacity,
            Zero
        )
    ),

    rhoCp_
    (
        IOobject
        (
            "rhoCp",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            dimDensity*dimSpecificHeatCapacity,
            Zero
        )
    ),

    rhoCpPhi_
    (
        IOobject
        (
            "rhoCpPhi",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            dimVelocity*dimDensity*dimSpecificHeatCapacity*dimArea,
            Zero
        )
    ),

    phi_
    (
        IOobject
        (
            "phi",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimVelocity*dimArea, Zero)
    ),
    nPotentialCorr_
    (
        mesh.time().controlDict().lookupOrDefault<label>("nPotentialCorr", 1)
    ),
    potentialCorrTol_
    (
        mesh.time().controlDict().lookupOrDefault<scalar>
        (
            "potentialCorrTol",
            0
        )
    )
{
    if (nPotentialCorr_ < 1)
    {
        FatalErrorInFunction
            << "nPotentialCorr must be at least 1, got " << nPotentialCorr_
            << exit(FatalError);
    }

    if (nPotentialCorr_ > 1)
    {
        Info<< "P5: up to " << nPotentialCorr_ << " potential-coupling passes"
            << " per electric solve (kinetics frozen), tolerance "
            << potentialCorrTol_ << " V" << endl;
    }

    Info << nl << "Creating regions: " << endl;
    regions_.reset
    (
        new regionTypeList
        (
            mesh_
        )
    );

    //- F6: construction complete (registry marker, nothing written)
    constructedMarker_.reset
    (
        new IOdictionary
        (
            IOobject
            (
                "fuelCellSystemConstructed",
                mesh_.time().constant(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            )
        )
    );

    //- Logged only when a boundary that uses the marker is present, so the
    //- log of every other case is unchanged
    bool usesMarker = false;
    forAll(T_.boundaryField(), patchi)
    {
        if (T_.boundaryField()[patchi].type() == "furnaceCouplingTemperature")
        {
            usesMarker = true;
        }
    }
    if (usesMarker)
    {
        Info<< "fuelCellSystem constructed (marker fuelCellSystemConstructed set)"
            << endl;
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fuelCellSystem::~fuelCellSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::fuelCellSystem::correct()
{
    //- regions
    regions_->correct();
}


void Foam::fuelCellSystem::setRDeltaT()
{
    regions_->setRDeltaT();
}


void Foam::fuelCellSystem::potentialCorrection(const scalar firstChange)
{
    //- Deferred correction of the potential coupling: each pass rewrites
    //  jSource from the latest potentials with eta, Nernst, J and jPrime
    //  frozen, so at convergence the pair solves the coupled linearised
    //  system. The outer relaxation of eta stays in correct().
    //  Galvanostatic one-update rule (user, FP-1 plan review): the
    //  terminal potential is set once per outer pass by the first solve;
    //  the inner passes hold it.
    label pass = 1;
    scalar change = firstChange;

    regions_->holdControl(true);

    while (pass < nPotentialCorr_)
    {
        const scalar dJ = regions_->correctCouplingMaxChange();
        change = regions_->solveElectricMaxChange();
        ++pass;

        //- Liveness (user D1): jSource must move when phi has moved
        Info<< "P5 pass " << pass << ": max |d jSource| = " << dJ
            << ", max |dphi| = " << change << endl;

        if (change < potentialCorrTol_)
        {
            break;
        }
    }

    regions_->holdControl(false);

    Info<< "potential coupling: " << pass << " passes, max |dphi| = "
        << change << endl;
}


void Foam::fuelCellSystem::solve()
{
    //- solve for regions
    regions_->solve();

    if (nPotentialCorr_ > 1)
    {
        potentialCorrection(GREAT);
    }
}


void Foam::fuelCellSystem::solveElectric()
{
    if (nPotentialCorr_ > 1)
    {
        potentialCorrection(regions_->solveElectricMaxChange());
    }
    else
    {
        regions_->solveElectric();
    }
}


void Foam::fuelCellSystem::mapToCell()
{
    //- T source, = 0

    Qdot_ *= 0.0;
    dQdotdT_ *= 0.0;
    QdotJoule_ *= 0.0;

    Info<< "\nMap to cell \n" << endl;

    regions_->mapToCell(*this);
}


void Foam::fuelCellSystem::mapFromCell()
{
    Info<< "\nMap from cell \n" << endl;

    regions_->mapFromCell(*this);
}


Foam::regionTypeList& Foam::fuelCellSystem::regions()
{
    return regions_();
}

// ************************************************************************* //

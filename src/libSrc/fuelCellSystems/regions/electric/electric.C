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
#include "electric.H"
#include "IOdictionary.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

#include "sigmaModelList.H"
#include "dissolvedModel.H"
#include "activationOverpotentialModel.H"

#include "fuelCellSystem.H"

#include "surfaceInterpolationScheme.H"
#include "snGradScheme.H"
#include "coupledFvPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{
    defineTypeNameAndDebug(electric, 0);

    addToRunTimeSelectionTable
    (
        regionType,
        electric,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionTypes::electric::electric
(
    const fvMesh& mesh,
    const word& regionName
)
:
    regionType(mesh, regionName),

    mesh_(mesh),
    i_
    (
        IOobject
        (
            "i",
            this->time().timeName(),
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        *this,
        dimensionedVector(dimCurrent/dimArea, Zero),
        zeroGradientFvPatchVectorField::typeName
    ),
    phi_
    (
        IOobject
        (
            activationOverpotentialModel::phiName,
            this->time().timeName(),
            *this,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        *this
    ),
    j_
    (
        IOobject
        (
            activationOverpotentialModel::jName,
            this->time().timeName(),
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        *this,
        dimensionedScalar(dimCurrent/dimVol, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    jPrime_
    (
        IOobject
        (
            "jPrime",
            this->time().timeName(),
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        *this,
        dimensionedScalar(j_.dimensions()/phi_.dimensions(), Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    jSource_
    (
        IOobject
        (
            "jSource",
            this->time().timeName(),
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        *this,
        dimensionedScalar(j_.dimensions(), Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    sigmaField_
    (
        IOobject
        (
            "sigma",
            this->time().timeName(),
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        *this,
        dimensionedScalar(j_.dimensions()/phi_.dimensions()*dimArea, 1.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    T_
    (
        IOobject
        (
            "T",
            this->time().timeName(),
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        *this,
        dimensionedScalar(dimTemperature, dict_.getOrDefault<scalar>("T0", 293.15)),
        zeroGradientFvPatchScalarField::typeName
    ),
    control_(dict_.lookupOrDefault<Switch>("control", false)),
    dissolveOnOff_(dict_.lookupOrDefault<Switch>("dissolveOnOff", false)),
    patchName_(word::null),
    zoneName_(word::null),
    cellZoneIDs_(),
    galvanostatic_(true),
    ibar_(nullptr),
    voltage_(nullptr),
    maxPhiStep_(0.2),
    faceBasedJoule_(false),
    jouleHeat_(nullptr)
{
    // F8 (main chat 2026-11-09 ~06:55/~08:05): read once at start-up; absent
    // key = upstream cell-centred Joule source, code path unchanged
    {
        const word mode
        (
            this->time().controlDict().getOrDefault<word>
            (
                "jouleHeat",
                "cellCentred"
            )
        );

        if (mode == "faceBased")
        {
            faceBasedJoule_ = true;
        }
        else if (mode != "cellCentred")
        {
            FatalIOErrorInFunction(this->time().controlDict())
                << "jouleHeat " << mode << " unknown; valid: cellCentred, "
                << "faceBased" << exit(FatalIOError);
        }

        if (faceBasedJoule_)
        {
            Info<< "jouleHeat faceBased (" << regionName << "): Joule source "
                << "from the two-point face current (harmonic sigma_f), "
                << "apportioned by d/sigma" << endl;

            jouleHeat_.reset
            (
                new volScalarField
                (
                    IOobject
                    (
                        "jouleHeat",
                        this->time().timeName(),
                        *this,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    *this,
                    dimensionedScalar(dimEnergy/dimVolume/dimTime, Zero),
                    zeroGradientFvPatchScalarField::typeName
                )
            );
        }
    }

    if (dissolveOnOff_)
    {
        dissolved_ = dissolvedModel::New(*this, dict_.subDict("dissolved"));
    }

    sigma_.reset
    (
        new sigmaModelList
        (
            *this,
            dict_.subDict("sigma")
        )
    );

    if (control_)
    {
        const dictionary& controlDict = dict_.subDict("galvanostatic");
        patchName_ = controlDict.get<word>("patchName");
        galvanostatic_ = controlDict.get<Switch>("active");

        if (galvanostatic_)
        {
            ibar_.reset(Function1<scalar>::New("ibar", controlDict, this));
            maxPhiStep_ =
                controlDict.lookupOrDefault<scalar>("maxPhiStep", 0.2);

            if (maxPhiStep_ <= 0)
            {
                FatalIOErrorInFunction(controlDict)
                    << "maxPhiStep must be greater than zero"
                    << exit(FatalIOError);
            }
        }
        else
        {
            voltage_.reset(Function1<scalar>::New("voltage", controlDict, this));
        }

        if (controlDict.found("zoneName"))
        {
            zoneName_ = controlDict.get<word>("zoneName");
            cellZoneIDs_ = this->cellZones().indices(zoneName_);
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionTypes::electric::~electric()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionTypes::electric::solvePotential()
{
    tmp<fvScalarMatrix> phiEqn
    (
      - fvm::laplacian(sigmaField_, phi_, "laplacian(sigma,phi)")
      + fvm::Sp(jPrime_, phi_)
     == jSource_
    );

    // A single reference is sufficient.  The positive Newton diagonal
    // normally removes the Laplacian null space in reacting regions.
    if (phi_.needReference() && nCells())
    {
        phiEqn->setReference(0, phi_[0], true);
    }

    phiEqn->solve();
}


Foam::scalar Foam::regionTypes::electric::integratedCurrent() const
{
    const scalarField source
    (
        jSource_.primitiveField()
      - jPrime_.primitiveField()*phi_.primitiveField()
    );

    scalar current = 0;

    if (cellZoneIDs_.size())
    {
        forAll(cellZoneIDs_, zoneI)
        {
            const labelList& cells = cellZones()[cellZoneIDs_[zoneI]];

            forAll(cells, i)
            {
                const label cellI = cells[i];
                current += source[cellI]*V()[cellI];
            }
        }
        reduce(current, sumOp<scalar>());
    }
    else
    {
        current = Foam::gSum(source*V());
    }

    return current;
}


void Foam::regionTypes::electric::solve()
{
    Info << "\nSolve for region " << name() << ":\n" << endl;

    if (control_ && galvanostatic_)
    {
        const label patchID = boundaryMesh().findPatchID(patchName_);

        if (patchID == -1)
        {
            FatalErrorInFunction
                << "Cannot find " << patchName_
                << "; please check patchName" << exit(FatalError);
        }

        const scalar area = Foam::gSum(magSf().boundaryField()[patchID]);
        const scalar oldPhi = Foam::gAverage(phi_.boundaryField()[patchID]);

        // The Newton-linearized equation is affine in the controlled
        // boundary value.  Two unit-response solves therefore give the
        // boundary potential required for the prescribed total current.
        phi_.boundaryFieldRef()[patchID] == scalar(0);
        solvePotential();
        const scalar current0 = integratedCurrent();
        const scalarField phi0(phi_.primitiveField());

        phi_.boundaryFieldRef()[patchID] == scalar(1);
        solvePotential();
        const scalar current1 = integratedCurrent();
        const scalarField phi1(phi_.primitiveField());

        const scalar conductance = current1 - current0;
        const scalar targetCurrent = ibar_->value(time().value())*area;
        scalar targetPhi = oldPhi;

        if (mag(conductance) > VSMALL)
        {
            targetPhi = (targetCurrent - current0)/conductance;
        }
        else
        {
            WarningInFunction
                << "The galvanostatic unit response is too small in region "
                << name() << "; retaining boundary potential " << oldPhi
                << endl;
        }

        const scalar deltaPhi = Foam::max
        (
            -maxPhiStep_,
            Foam::min(maxPhiStep_, targetPhi - oldPhi)
        );
        const scalar boundaryPhi = oldPhi + deltaPhi;

        phi_.primitiveFieldRef() =
            phi0 + boundaryPhi*(phi1 - phi0);
        phi_.boundaryFieldRef()[patchID] == boundaryPhi;
        phi_.correctBoundaryConditions();

        const scalar current = current0 + boundaryPhi*conductance;
        Info<< "ibar: " << current/area
            << '\t' << "voltage: " << boundaryPhi << endl;
    }
    else
    {
        solvePotential();

        if (control_)
        {
            const label patchID = boundaryMesh().findPatchID(patchName_);
            const scalar area = Foam::gSum(magSf().boundaryField()[patchID]);
            Info<< "ibar: " << integratedCurrent()/area
                << '\t' << "voltage: "
                << Foam::gAverage(phi_.boundaryField()[patchID]) << endl;
        }
    }

    i_ = -sigmaField_ * fvc::grad(phi_);
    i_.correctBoundaryConditions();

    if (dissolved_.valid())
    {
        dissolved_->solve();
    }
}


void Foam::regionTypes::electric::setRDeltaT()
{
    //- Do nothing, add if necessary.
}


void Foam::regionTypes::electric::correct()
{
    sigma_->correct(sigmaField_);
    sigmaField_.correctBoundaryConditions();

    // Models that do not provide a Jacobian retain the original explicit
    // source treatment.
    if (Foam::gMax(jPrime_.primitiveField()) <= VSMALL)
    {
        jSource_ = j_;
    }

    if (control_ && !galvanostatic_)
    {
        const label patchID = boundaryMesh().findPatchID(patchName_);

        if (patchID == -1)
        {
            FatalErrorInFunction
                << "Cannot find " << patchName_
                << "; please check patchName" << exit(FatalError);
        }

        phi_.boundaryFieldRef()[patchID] == voltage_->value(time().value());
    }

    if (dissolved_.valid())
    {
        dissolved_->correct();
    }
}


void Foam::regionTypes::electric::mapToCell
(
    fuelCellSystem& fuelCell
)
{
    Info << "Map " << name() << " to Cell " << nl << endl;

    //- heat source; F8: faceBased replaces the cell-centred form only when
    //- selected (controlDict jouleHeat), otherwise the upstream expression
    volScalarField heatSource
    (
        faceBasedJoule_ ? faceJouleHeat() : (i_ & i_)/sigmaField_
    );

    if (faceBasedJoule_)
    {
        jouleHeat_->primitiveFieldRef() = heatSource.primitiveField();
        jouleHeat_->correctBoundaryConditions();

        Info<< "Joule " << name() << " (faceBased): sum q V = "
            << gSum(heatSource.primitiveField()*this->V()) << " W" << endl;
    }

    volScalarField heatSource0
    (
        IOobject
        (
            "heatSource",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar(dimEnergy/dimVolume/dimTime, Zero)
    );

    heatSource0.rmap(heatSource, cellMapIO_);

    fuelCell.Qdot() += heatSource0;

    //- F8 / review of ea54274: the Joule part alone (not written; logged by
    //- the heatSourceScale block only when heatSourceScale is present)
    fuelCell.QdotJoule() += heatSource0;
}


Foam::tmp<Foam::volScalarField>
Foam::regionTypes::electric::faceJouleHeat() const
{
    const fvMesh& m = *this;

    // the phi matrix's own scheme: "Gauss <interpolation> <snGrad>", read as
    // laplacianScheme reads it (finiteVolume laplacianScheme.H)
    ITstream is(m.laplacianScheme("laplacian(sigma,phi)"));
    is.rewind();
    const word gauss(is);

    if (gauss != "Gauss" || is.eof())
    {
        FatalErrorInFunction
            << "jouleHeat faceBased (" << name() << "): laplacian(sigma,phi) "
            << "must be 'Gauss <interpolation> <snGrad>', got "
            << m.laplacianScheme("laplacian(sigma,phi)") << exit(FatalError);
    }

    tmp<surfaceInterpolationScheme<scalar>> interp
    (
        surfaceInterpolationScheme<scalar>::New(m, is)
    );
    tmp<fv::snGradScheme<scalar>> sng(fv::snGradScheme<scalar>::New(m, is));

    // two-point face current along Sf [A] (outward on boundary faces)
    const surfaceScalarField sigmaf(interp().interpolate(sigmaField_));
    const surfaceScalarField F(-sigmaf*m.magSf()*sng().snGrad(phi_));

    tmp<volScalarField> tq
    (
        new volScalarField
        (
            IOobject
            (
                "jouleHeatFace",
                this->time().timeName(),
                m,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            m,
            dimensionedScalar(dimEnergy/dimVolume/dimTime, Zero),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    // accumulate W per cell, divide by V at the end
    scalarField& q = tq.ref().primitiveFieldRef();

    const labelUList& own = m.owner();
    const labelUList& nei = m.neighbour();
    const volVectorField& C = m.C();
    const surfaceVectorField& Cf = m.Cf();
    const surfaceVectorField& Sf = m.Sf();
    const surfaceScalarField& magSf = m.magSf();
    const scalarField& sig = sigmaField_.primitiveField();
    const scalarField& ph = phi_.primitiveField();

    forAll(own, facei)
    {
        const label P = own[facei];
        const label N = nei[facei];
        const vector n = Sf[facei]/magSf[facei];
        const scalar rP = mag(n & (Cf[facei] - C[P]))/sig[P];
        const scalar rN = mag(n & (C[N] - Cf[facei]))/sig[N];
        const scalar D = -F[facei]*(ph[N] - ph[P]);

        q[P] += D*rP/(rP + rN);
        q[N] += D*rN/(rP + rN);
    }

    forAll(F.boundaryField(), patchi)
    {
        const fvPatch& p = m.boundary()[patchi];
        const fvsPatchScalarField& Fp = F.boundaryField()[patchi];
        const fvPatchScalarField& php = phi_.boundaryField()[patchi];
        const labelUList& fc = p.faceCells();

        if (p.coupled())
        {
            // processor faces: own half only, neighbour from the other side
            const scalarField phN(php.patchNeighbourField());
            const scalarField sigN
            (
                sigmaField_.boundaryField()[patchi].patchNeighbourField()
            );
            const vectorField nf(p.nf());
            const vectorField Cfp(p.Cf());
            const vectorField d(p.delta());

            forAll(fc, i)
            {
                const label P = fc[i];
                const scalar dP = mag(nf[i] & (Cfp[i] - C[P]));
                const scalar dN = mag(nf[i] & d[i]) - dP;
                const scalar rP = dP/sig[P];
                const scalar rN = dN/sigN[i];
                const scalar D = -Fp[i]*(phN[i] - ph[P]);

                q[P] += D*rP/(rP + rN);
            }
        }
        else
        {
            // physical boundary: F (phi_P - phi_b), zero on zeroGradient
            forAll(fc, i)
            {
                q[fc[i]] += Fp[i]*(ph[fc[i]] - php[i]);
            }
        }
    }

    q /= m.V();
    tq.ref().correctBoundaryConditions();

    return tq;
}


void Foam::regionTypes::electric::mapFromCell
(
    fuelCellSystem& fuelCell
)
{
    Info << "Map " << name() << " from Cell " << nl << endl;

    scalarField& T = T_;

    forAll(T, cellI)
    {
        T[cellI] = fuelCell.T()[cellMapIO_[cellI]];
    }

    T_.correctBoundaryConditions();
}

// ************************************************************************* //

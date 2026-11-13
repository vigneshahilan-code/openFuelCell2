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

#include "regionTypeList.H"
#include "volFields.H"
#include "IOdictionary.H"
#include "fuelCellSystem.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionTypeList::regionTypeList
(
    const fvMesh& mesh
)
:
    PtrList<regionType>(),
    mesh_(mesh),
    region_(mesh.time())
{
    reset(region_);

    active(true);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionTypeList::~regionTypeList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::regionTypeList::active(const bool warn) const
{
    bool a = false;
    forAll(*this, i)
    {
        a = a || this->operator[](i).active();
    }

    if (warn && this->size() && !a)
    {
        Info<< "Cannot find active regions..." << endl;
    }

    return a;
}


void Foam::regionTypeList::reset(const regionProperties& rp)
{
    wordList regionNames;

    forAllConstIter(HashTable<wordList>, rp, iter)
    {
        const wordList& regions = iter();

        forAll(regions, regionI)
        {
            if (!regionNames.found(regions[regionI]))
            {
                regionNames.append(regions[regionI]);
            }
        }
    }

    this->setSize(regionNames.size());

    //- Creating the regions in order:
    //- fluid-->solid-->electric
    //- It is found solving in this order
    //- makes the solution more stable
    //
    wordList modelTypes = {"fluid", "solid", "electric"};

    //- Add other types in needed..
    for (auto& typeI : regionNames)
    {
        modelTypes.push_uniq(typeI);
    }

    label i = 0;

    for (auto& modelType : modelTypes)
    {
        if (rp.found(modelType))
        {
            for (auto& region : rp[modelType])
            {
                Info << "Creating " << region << endl;

                this->set
                (
                    i++,
                    regionType::New
                    (
                        mesh_,
                        region,
                        modelType
                    )
                );
            }
        }
    }
}


void Foam::regionTypeList::correct()
{
    forAll(*this, i)
    {
        this->operator[](i).correct();
    }
}


void Foam::regionTypeList::setRDeltaT()
{
    forAll(*this, i)
    {
        this->operator[](i).setRDeltaT();
    }
}


void Foam::regionTypeList::solve()
{
    forAll(*this, i)
    {
        this->operator[](i).solve();
    }
}


void Foam::regionTypeList::solveElectric()
{
    forAll(*this, i)
    {
        if (this->operator[](i).type() == "electric")
        {
            this->operator[](i).solve();
        }
    }
}


void Foam::regionTypeList::correctCoupling()
{
    forAll(*this, i)
    {
        this->operator[](i).correctCoupling();
    }
}


Foam::scalar Foam::regionTypeList::correctCouplingMaxChange()
{
    PtrList<scalarField> before(this->size());

    forAll(*this, i)
    {
        regionType& region = this->operator[](i);

        if (region.type() == "electric")
        {
            before.set
            (
                i,
                new scalarField
                (
                    region.lookupObject<volScalarField>("jSource").primitiveField()
                )
            );
        }
    }

    correctCoupling();

    scalar change = 0;

    forAll(*this, i)
    {
        if (before.set(i))
        {
            const scalarField& now =
                this->operator[](i).lookupObject<volScalarField>("jSource").primitiveField();

            const scalar dJ = gMax(mag(now - before[i])());

            //- Per electric region (user D1, packet 08:04): the anode
            //  (electron), cathode (electron) and ion regions each logged,
            //  so one dead electrode cannot hide behind the other
            Info<< "P5 liveness " << this->operator[](i).name()
                << ": max |d jSource| = " << dJ << endl;

            change = max(change, dJ);
        }
    }

    return change;
}


void Foam::regionTypeList::holdControl(const bool hold)
{
    forAll(*this, i)
    {
        this->operator[](i).holdControl(hold);
    }
}


Foam::scalar Foam::regionTypeList::solveElectricMaxChange()
{
    scalar change = 0;

    forAll(*this, i)
    {
        regionType& region = this->operator[](i);

        if (region.type() == "electric")
        {
            const volScalarField& phi =
                region.lookupObject<volScalarField>("phi");
            const scalarField phiOld(phi.primitiveField());

            region.solve();

            change = max(change, gMax(mag(phi.primitiveField() - phiOld)()));
        }
    }

    return change;
}


void Foam::regionTypeList::mapToCell
(
    fuelCellSystem& fuelCell
)
{
    forAll(*this, i)
    {
        this->operator[](i).mapToCell(fuelCell);
    }
}


void Foam::regionTypeList::mapFromCell
(
    fuelCellSystem& fuelCell
)
{
    forAll(*this, i)
    {
        this->operator[](i).mapFromCell(fuelCell);
    }
}

// ************************************************************************* //

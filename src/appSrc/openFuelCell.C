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

Application
    openFuelCell

Description
    Steady (transient) solver for multiple region heat and mass transfer problems:
    e.g. fuel cells and electrolyzers...
    The solver works for several types of fuel cell and electrolyzer applications.
    The solver also provides a framework for multiple region and multiphysical modeling and simulations.

Contributors
    Shidong Zhang (s.zhang@fz-juelich.de)
    Steven B. Beale (s.beale@fz-juelich.de)
    Steffen Hess (s.hess@fz-juelich.de)

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "regionTypeList.H"
#include "fuelCellSystem.H"
#include "phaseSystem.H"
#include "regionCourantNo.H"
#include "phaseCompressibleTurbulenceModel.H"
#include "pimpleControl.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Steady (transient) solver for multiple region heat and mass transfer problems: "
        "e.g. fuel cells and electrolyzers..."
        " The solver works for several types of fuel cell and electrolyzer applications."
        " The solver also provides a framework for multiple region and multiphysical modeling and simulations."
    );

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "createTimeControls.H"

    if (!LTS)
    {
        #include "multiRegionCourantNo.H"
        #include "setInitialMultiRegionDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    const label nCouplingIter =
        runTime.controlDict().lookupOrDefault<label>("nCouplingIter", 3);

    if (nCouplingIter < 1)
    {
        FatalIOErrorInFunction(runTime.controlDict())
            << "nCouplingIter must be at least 1"
            << exit(FatalIOError);
    }

    // Read once at start-up (not run-time modifiable). Default true keeps
    // the upstream behaviour; false skips the energy equation so T stays at
    // its initial field (isothermal runs).
    const bool solveEnergy =
        runTime.controlDict().lookupOrDefault<bool>("solveEnergy", true);

    if (!solveEnergy)
    {
        Info<< "solveEnergy false: energy equation not solved,"
            << " T held at its initial field\n" << endl;
    }

    fuelCell.mapToCell();

    while (runTime.run())
    {
        runTime++;
        Info<< "\nTime = " << runTime.timeName() << nl << endl;

        if (!LTS)
        {
            #include "createTimeControls.H"
            #include "multiRegionCourantNo.H"
            #include "setMultiRegionDeltaT.H"
        }
        else
        {
            fuelCell.setRDeltaT();
        }

        fuelCell.mapFromCell();

        for (label couplingIter = 0; couplingIter < nCouplingIter; ++couplingIter)
        {
            Info<< "\nElectrochemical coupling iteration "
                << couplingIter + 1 << '/' << nCouplingIter << nl << endl;

            fuelCell.correct();

            if (couplingIter == 0)
            {
                fuelCell.solve();
            }
            else
            {
                fuelCell.solveElectric();
            }
        }

        fuelCell.mapToCell();

        if (solveEnergy)
        {
            #include "EEqns.H"
        }
        else
        {
            Info<< "T (not solved)"
                << ": min = " << Foam::min(TCell.primitiveField())
                << ", mean = " << Foam::average(TCell.primitiveField())
                << ", max = " << Foam::max(TCell.primitiveField())
                << endl;
        }

        runTime.write();

        Info<< "ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n\n" << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

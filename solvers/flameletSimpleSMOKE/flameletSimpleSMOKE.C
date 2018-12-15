/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    rhoSimpleFoam

Description
    Steady-state SIMPLE solver for laminar or turbulent RANS flow of
    compressible fluids.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "flameletSMOKEThermo.H"

#if OPENFOAM_VERSION < 40
#include "RASModel.H"
#else
#include "turbulentFluidThermoModel.H"
#endif
#include "bound.H"
#include "simpleControl.H"
#if OPENFOAM_VERSION < 40
#include "fvIOoptionList.H"
#else
#if OPENFOAM_VERSION >= 50
#include "pressureControl.H"
#endif
#include "fvOptions.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#endif
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readMassFlowProperties.H"
    #include "readGravitationalAcceleration.H"

    #if OPENFOAM_VERSION < 40
    simpleControl simple(mesh);
    #else
    #include "createControl.H"
    #endif

    #include "createFields.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"

    #if OPENFOAM_VERSION >= 40
    turbulence->validate();
    #endif

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
    #if OPENFOAM_VERSION >= 60
    while (simple.loop(runTime))
    #else
    while (simple.loop())
    #endif  
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Pressure-velocity SIMPLE corrector
        {
		#include "UEqn.H"

    		#if OPENFOAM_VERSION >= 40
		if (simple.consistent())
		{
			#include "pcEqn.H"
		}
		else
		#endif
		{
			#include "pEqn.H"
		}

		#include "ZEqn.H"
		#include "HEqn.H"
        }

        turbulence->correct();

        runTime.write();

	#include "writeMassFlow.H"
	#include "outputVariables.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

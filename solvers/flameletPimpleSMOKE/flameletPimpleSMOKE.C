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
    rhoPimpleFoam

Description
    Transient solver for laminar or turbulent flow of compressible fluids
    for HVAC and similar applications.

    Uses the flexible PIMPLE (PISO-SIMPLE) solution for time-resolved and
    pseudo-transient simulations.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "flameletSMOKEThermo.H"

#if OPENFOAM_VERSION < 40
#include "RASModel.H"
#else
#include "turbulentFluidThermoModel.H"
#endif
#include "bound.H"
#include "pimpleControl.H"
#if OPENFOAM_VERSION < 40
#include "fvIOoptionList.H"
#else
#include "pressureControl.H"
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
    pimpleControl pimple(mesh);
    #else
    #include "createControl.H"
    #include "createTimeControls.H"
    #endif

    #include "createFields.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"
    #include "readTimeControls.H"

    #if OPENFOAM_VERSION >= 40
    turbulence->validate();
    if (!LTS)
    #endif
    {
        #include "compressibleCourantNo.H"
        #include "setInitialDeltaT.H"
    }
 
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
	#include "readTimeControls.H"

	#if OPENFOAM_VERSION >= 40
        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
	#endif
        {
            #include "compressibleCourantNo.H"
            #include "setDeltaT.H"
        }

        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

	#if OPENFOAM_VERSION >= 40
	if (pimple.nCorrPIMPLE() <= 1)
	#endif
        {
            #include "rhoEqn.H"
        }

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "UEqn.H"
            #include "HEqn.H"
	    #include "ZEqn.H"

	    // --- Pressure corrector loop
            while (pimple.correct())
            {
		#if OPENFOAM_VERSION >= 40
                if (pimple.consistent())
                {
                    #include "pcEqn.H"
                }
                else
		#endif
                {
                    #include "pEqn.H"
                }
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }

	    rho = thermo.rho();
        }

        runTime.write();
	
	#include "writeMassFlow.H"
	#include "outputVariables.H"

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //


{
    volScalarField rAU(1.0/UEqn.A());
    surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
    tUEqn.clear();

    bool closedVolume = false;

    surfaceScalarField phiHbyA("phiHbyA", fvc::flux(rho*HbyA));
    MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p, rho, U, phiHbyA, rhorAUf, MRF);

    if (simple.transonic())
    {
        surfaceScalarField phid
        (
            "phid",
            (fvc::interpolate(psi)/fvc::interpolate(rho))*phiHbyA
        );
        phiHbyA -= fvc::interpolate(p)*phid;

        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix pEqn
            (
                fvc::div(phiHbyA)
              + fvm::div(phid, p)
              - fvm::laplacian(rhorAUf, p)
              ==
                fvOptions(psi, p, rho.name())
            );

            // Relax the pressure equation to ensure diagonal-dominance
            pEqn.relax();
	
	    #if OPENFOAM_VERSION >= 50
            pEqn.setReference
            (
                pressureControl.refCell(),
                pressureControl.refValue()
            );
	    #endif

            pEqn.solve();

            if (simple.finalNonOrthogonalIter())
            {
                phi = phiHbyA + pEqn.flux();
            }
        }
    }
    else
    {
        closedVolume = adjustPhi(phiHbyA, U, p);

        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix pEqn
            (
                fvc::div(phiHbyA)
              - fvm::laplacian(rhorAUf, p)
              ==
                fvOptions(psi, p, rho.name())
            );

            #if OPENFOAM_VERSION >=50
	    pEqn.setReference
            (
                pressureControl.refCell(),
                pressureControl.refValue()
            );
	    #endif

            pEqn.solve();

            if (simple.finalNonOrthogonalIter())
            {
                phi = phiHbyA + pEqn.flux();
            }
        }
    }


    #include "incompressible/continuityErrs.H"

    // Explicitly relax pressure for momentum corrector
    p.relax();

    U = HbyA - rAU*fvc::grad(p);
    U.correctBoundaryConditions();
    fvOptions.correct(U);

    #if OPENFOAM_VERSION >= 50
    pressureControl.limit(p);
    #endif

    // For closed-volume cases adjust the pressure and density levels
    // to obey overall mass continuity
    if (closedVolume)
    {
        p += (initialMass - fvc::domainIntegrate(psi*p))
            /fvc::domainIntegrate(psi);
    }

    p.correctBoundaryConditions();

    rho = thermo.rho();

    if (!simple.transonic())
    {
        rho.relax();
    }


    if(monitoring != true)
    {
        Info<< "rho max/min : " << max(rho).value() << " " << min(rho).value() << endl;
    }
}

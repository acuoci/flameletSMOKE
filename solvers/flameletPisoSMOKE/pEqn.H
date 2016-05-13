rho = thermo.rho();

volScalarField rUA = 1.0/UEqn.A();
U = rUA*UEqn.H();

if (transonic)
{

    surfaceScalarField phid
    (
        "phid",
        fvc::interpolate(psi)
       *(
            	(fvc::interpolate(U) & mesh.Sf())
		#if OPENFOAM_VERSION == 22
			+ fvc::ddtPhiCorr(rUA, rho, U, phi)
		#elif OPENFOAM_VERSION == 23
			+ fvc::interpolate(rUA)*fvc::ddtCorr(rho, U, phi)
		#elif OPENFOAM_VERSION == 24
			+ fvc::interpolate(rUA)*fvc::ddtCorr(rho, U, phi)
		#else
			Compilation_Error_You_Did_Not_Define_The_OPENFOAM_VERSION_Macro
		#endif
        )
    );

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi, p)
          + fvm::div(phid, p)
          - fvm::laplacian(rho*rUA, p)
        );

        pEqn.solve();

        if (nonOrth == nNonOrthCorr)
        {
            phi == pEqn.flux();
        }
    }
}
else
{
    phi =
        fvc::interpolate(rho)*
        (
            	(fvc::interpolate(U) & mesh.Sf())
		#if OPENFOAM_VERSION == 22
			+ fvc::ddtPhiCorr(rUA, rho, U, phi)
		#elif OPENFOAM_VERSION == 23
			+ fvc::interpolate(rUA)*fvc::ddtCorr(rho, U, phi)
		#elif OPENFOAM_VERSION == 24
			+ fvc::interpolate(rUA)*fvc::ddtCorr(rho, U, phi)
		#else
			Compilation_Error_You_Did_Not_Define_The_OPENFOAM_VERSION_Macro
		#endif   
        );

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi, p)
          + fvc::div(phi)
          - fvm::laplacian(rho*rUA, p)
        );

        pEqn.solve();

        if (nonOrth == nNonOrthCorr)
        {
            phi += pEqn.flux();
        }
    }
}

//#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"

U -= rUA*fvc::grad(p);
U.correctBoundaryConditions();

dpdt = fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);

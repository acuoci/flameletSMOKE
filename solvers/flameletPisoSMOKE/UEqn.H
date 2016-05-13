    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(phi, U)
      + turbulence->divDevRhoReff(U)
     ==
        rho*g
      + fvOptions(rho, U)
    );

    if (momentumPredictor)
    {
        solve(UEqn == -fvc::grad(p));
    }

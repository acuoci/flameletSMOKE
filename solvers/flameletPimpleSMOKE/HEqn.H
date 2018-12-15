{
	if (adiabaticMode == true)
	{
		fvScalarMatrix hEqn
		(
			fvm::ddt(rho, H)
	              + fvm::div(phi, H)
	              - fvm::laplacian(turbulence->muEff()/sigmat, H)
		);

		
		hEqn.relax();
    		fvOptions.constrain(hEqn);	
		#if OPENFOAM_VERSION >=60
		hEqn.solve("H");
		#else
		hEqn.solve(mesh.solver("H"));
		#endif
	}

	else
	{
		if (radiationMode == true)
		{
			Qrad = 4. * sigmaSB * as * (pow(thermo.T(),4) - pow(Tenv,4));
		}

		else if(radiationMode == false && switchQoff == false)
		{
			Qrad = 0 * 4. * sigmaSB * as * (pow(thermo.T(),4) - pow(Tenv,4));
			switchQoff = true;
		}

		fvScalarMatrix hEqn
		(
                        fvm::ddt(rho,H)
	              + fvm::div(phi, H)
	              - fvm::laplacian(turbulence->muEff()/sigmat, H) + Qrad
		);

		hEqn.relax();
    		fvOptions.constrain(hEqn);
		#if OPENFOAM_VERSION >=60
		hEqn.solve("H");
		#else
		hEqn.solve(mesh.solver("H"));
		#endif
	}

	thermo.correct();

        Info<< " * Tgas min/max = " << min(thermo.T()).value() << ", " << max(thermo.T()).value() << endl;
}

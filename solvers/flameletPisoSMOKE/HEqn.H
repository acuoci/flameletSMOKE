{
	if (adiabaticMode == true)
	{
		fvScalarMatrix hEqn
		(
			fvm::ddt(rho, H)
	              + fvm::div(phi, H)
	              - fvm::laplacian(turbulence->muEff()/sigmat, H)
		);

		hEqn.solve(mesh.solver("H"));
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

		hEqn.solve(mesh.solver("H"));
	}

	thermo.correct();
}

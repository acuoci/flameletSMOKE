	// Solve the Momentum equation

	#if OPENFOAM_VERSION < 40

		tmp<fvVectorMatrix> UEqn
		(
		fvm::ddt(rho, U)
		+ fvm::div(phi, U)
		+ turbulence->divDevRhoReff(U)
		==
		rho*g
		+ fvOptions(rho, U)
		);

		UEqn().relax();

		fvOptions.constrain(UEqn());

		if (pimple.momentumPredictor())
		{
			solve(UEqn() == -fvc::grad(p));

			fvOptions.correct(U);
			K = 0.5*magSqr(U);
		}

	#else

		MRF.correctBoundaryVelocity(U);

		tmp<fvVectorMatrix> tUEqn
		(
			fvm::ddt(rho, U) + fvm::div(phi, U)
			+ MRF.DDt(rho, U)
			+ turbulence->divDevRhoReff(U)
			==
                        rho*g + 
			fvOptions(rho, U)
		);
		fvVectorMatrix& UEqn = tUEqn.ref();

		UEqn.relax();

		fvOptions.constrain(UEqn);

		if (pimple.momentumPredictor())
		{
			solve(UEqn == -fvc::grad(p));

			fvOptions.correct(U);
			K = 0.5*magSqr(U);
		}

	#endif    


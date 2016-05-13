//-  Solving equation of variance of mixture fraction
fvScalarMatrix ZEqn
(
	(
		fvm::div(phi, Z)
//	      - fvm::Sp(fvc::div(phi), Z)
	      - fvm::laplacian(turbulence->muEff()/sigmat, Z)
	)
);

ZEqn.relax();
ZEqn.solve(mesh.solver("Z"));

if (varianceEquation == false)
{
	//- Solving equation of variance of mixture fraction
	fvScalarMatrix ZvarEqn
	(
		(
		    fvm::div(phi, Zvar)
//		  - fvm::Sp(fvc::div(phi), Zvar)
		  - fvm::laplacian(turbulence->mut()/sigmat, Zvar)
		  - Cg*turbulence->mut()*magSqr(fvc::grad(Z))
		  + Cd*rho*turbulence->epsilon()/(turbulence->k()+k_small)*Zvar
		)
	);

	ZvarEqn.relax();
	ZvarEqn.solve(mesh.solver("Zvar"));
}
else
{
	Info << "Algebraic equation for Zvar" << endl;
	scalarField& ZvarCells = Zvar.internalField();
	ZvarCells = Cg*turbulence->mut()*magSqr(fvc::grad(Z))*turbulence->k()/(Cd*rho*(turbulence->epsilon()+epsilon_small));
}

//- Updating stoichiometric scalar dissipation rate
chi_st=Cx*turbulence->epsilon()/(turbulence->k()+k_small)*Zvar;



	Info<< "Reading thermophysical properties\n" << endl;

	autoPtr<flameletSMOKEThermo> pThermo
	(
		flameletSMOKEThermo::New(mesh)
	);
	flameletSMOKEThermo& thermo = pThermo();

	volScalarField rho
	(
		IOobject
		(
			"rho",
			runTime.timeName(),
			mesh,
			IOobject::READ_IF_PRESENT,
			IOobject::AUTO_WRITE
		),
		thermo.rho()
	);

	volScalarField& p = thermo.p();
	const volScalarField& T = thermo.T();
	const volScalarField& psi = thermo.psi();

	Info << "Preparing additional scalar fields" << endl;
	Info << " * thermo.Z()" 	<< endl;
	Info << " * thermo.Zvar()" 	<< endl;
	Info << " * thermo.chi_st()" 	<< endl;
	Info << " * thermo.H()" 	<< endl;
	Info << " * thermo.as()" 	<< endl;
	Info << endl;

	volScalarField& Z     	= thermo.Z();		// mixture fraction [-]
	volScalarField& Zvar  	= thermo.Zvar();	// variance of mixture fraction [-]
	volScalarField& chi_st 	= thermo.chi_st();	// scalar dissipation rate [1/s]
	volScalarField& H  	= thermo.H();    	// enthalpy [J/kg]
	volScalarField& as  	= thermo.as();   	// mean absorption coefficient [1/m]

	Info<< "Reading field U\n" << endl;
	volVectorField U
	(
		IOobject
		(
			"U",
			runTime.timeName(),
			mesh,
			IOobject::MUST_READ,
			IOobject::AUTO_WRITE
		),
		mesh
	);

    	#include "compressibleCreatePhi.H"

    	label pRefCell = 0;
    	scalar pRefValue = 0.0;
    	setRefCell(p, simple.dict(), pRefCell, pRefValue);

	#if OPENFOAM_VERSION == 30

	mesh.setFluxRequired(p.name());

	dimensionedScalar rhoMax
	(
		dimensionedScalar::lookupOrDefault
		(
		"rhoMax",
		simple.dict(),
		dimDensity,
		GREAT
		)
	);

	dimensionedScalar rhoMin
	(
		dimensionedScalar::lookupOrDefault
		(
		"rhoMin",
		simple.dict(),
		dimDensity,
		0
		)
	);

	#else

    	dimensionedScalar rhoMax(simple.dict().lookup("rhoMax"));
    	dimensionedScalar rhoMin(simple.dict().lookup("rhoMin"));

	dimensionedScalar pMin
	(
		mesh.solutionDict().subDict("SIMPLE").lookup("pMin")
	);

	#endif

	Info<< "Creating turbulence model\n" << endl;
	autoPtr<compressible::RASModel> turbulence
	(
		compressible::RASModel::New
		(
			rho,
			U,
			phi,
			thermo
		)
	);

	Info<< "Reading flameletProperties dictionary\n" << endl;
	IOdictionary flameletsProperties_
	(
		IOobject
		(
			"flameletsProperties",
			U.time().constant(),
			U.db(),
			IOobject::MUST_READ,
			IOobject::NO_WRITE
		)
	);

    	dimensionedScalar initialMass = fvc::domainIntegrate(rho);

    	//- Monitoring
    	Switch monitoring(flameletsProperties_.lookup("monitoring"));

    	//- Semi-empirical constants in mixture fraction variance equation
	scalar Cg     = readScalar(flameletsProperties_.lookup("Cg"));
	scalar Cd     = readScalar(flameletsProperties_.lookup("Cd"));

	//- Semi-empirical constant for scalar dissipation rate evaluation
	scalar Cx     = readScalar(flameletsProperties_.lookup("Cx"));

	//- Prandtl number for mixture fraction, variance of mixture fraction and enthalpy diffusion terms
	scalar sigmat_ = readScalar(flameletsProperties_.lookup("sigmat"));
	dimensionedScalar sigmat("sigmat", dimensionSet(0,0,0,0,0,0,0), sigmat_);

	//- Minimum values (to be used for scalar dissipation rate evaluation)
	dimensionedScalar k_small( flameletsProperties_.lookup("k_small") );
	dimensionedScalar epsilon_small( flameletsProperties_.lookup("epsilon_small") );

	//- Environment temperature (non adiabatic consitions)
	dimensionedScalar Tenv ( flameletsProperties_.lookup("Tenv") );

	//- Variance of mixture fraction: on=transport equation | off=algebraic equation
	Switch varianceEquation(flameletsProperties_.lookup("varianceEquation"));

	//- Adiabatic simulation: on=adiabatic | off=non adiabatic
	Switch adiabaticMode(flameletsProperties_.lookup("adiabaticMode"));
	Switch radiationMode(flameletsProperties_.lookup("radiationMode"));

	bool switchQoff = false;

	Info<< "Preparing field Qrad (radiative heat transfer)\n" << endl;
	volScalarField Qrad
	(
		IOobject
		(
			"Qrad",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		mesh,
		dimensionedScalar("qradiation", dimensionSet(1,-1,-3,0,0,0,0) , 0.0)
	);

	dimensionedScalar sigmaSB("stefan_bolzmann", dimensionSet(1,0,-3,-4,0,0,0), 5.6704e-8);

	// Check the gravity
	if ( (mag(g)).value() != 0.)
	{
		Info << endl;
		Info << "------------------------------------------------------------------" << endl;
		Info << " WARNING: gravity is turned on                                    " << endl;
		Info << "          Please consider that this solver is not conceived       " << endl;
 		Info << "          for bouyancy driven flows. If you think gravity         " << endl;
 		Info << "          plays a non-negligible role, consider the possibility   " << endl;
  		Info << "          to use the laminarBouyantSimpleSMOKE solver             " << endl;
		Info << "------------------------------------------------------------------" << endl;
		Info << endl;
		abort();
	}

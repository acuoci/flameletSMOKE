/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "pdfFlameletSMOKEThermo.H"
#include "Time.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicFlameletSMOKEThermo, class MixtureType>
void Foam::pdfFlameletSMOKEThermo<BasicFlameletSMOKEThermo, MixtureType>::calculate()
{
	const scalarField& pCells	= this->p_.internalField();
	const scalarField& HCells 	= this->H_.internalField();


	const scalarField& ZCells 	= this->Z_.internalField();
	const scalarField& RhoReynolds 	= this->density_reynolds_.internalField();
	const scalarField& muFavre 	= this->mu_favre_.internalField();
	const scalarField& alphaFavre 	= this->alpha_favre_.internalField();

	scalarField& psiCells 		= this->psi_.internalField();
	scalarField& muCells 		= this->mu_.internalField();
	scalarField& alphaCells 	= this->alpha_.internalField();
	scalarField& defectCells 	= this->defect_.internalField();

	forAll(ZCells, celli)
	{
		psiCells[celli] 	= RhoReynolds[celli]/pCells[celli];
		muCells[celli] 		= muFavre[celli];
		alphaCells[celli]	= alphaFavre[celli];
		defectCells[celli]	= HCells[celli] - (HOxidizer+ZCells[celli]*(HFuel-HOxidizer));
	}

	// Boundaries
	forAll(this->T_.boundaryField(), patchi)
	{
		const fvPatchScalarField& pp			= this->p_.boundaryField()[patchi];
		const fvPatchScalarField& pH 			= this->H_.boundaryField()[patchi];
		const fvPatchScalarField& pZ 			= this->Z_.boundaryField()[patchi];
		const fvPatchScalarField& pRhoReynolds 		= this->density_reynolds_.boundaryField()[patchi];
		const fvPatchScalarField& pmuFavre 		= this->mu_favre_.boundaryField()[patchi];
		const fvPatchScalarField& palphaFavre 		= this->alpha_favre_.boundaryField()[patchi];

		fvPatchScalarField& pT 		= this->T_.boundaryField()[patchi];
		fvPatchScalarField& ppsi 	= this->psi_.boundaryField()[patchi];
		fvPatchScalarField& pdefect	= this->defect_.boundaryField()[patchi];
		fvPatchScalarField& pmu 	= this->mu_.boundaryField()[patchi];
		fvPatchScalarField& palpha 	= this->alpha_.boundaryField()[patchi];


		if (pT.fixesValue())
		{
			forAll(pT, facei)
			{
				ppsi[facei] 	= pRhoReynolds[facei]/pp[facei];
				pmu[facei] 	= pmuFavre[facei];
				palpha[facei] 	= palphaFavre[facei];
				pdefect[facei]  = pH[facei] - (HOxidizer+pZ[facei]*(HFuel-HOxidizer));
			}
		}
		else
		{
			forAll(pT, facei)
			{

				ppsi[facei] 	= pRhoReynolds[facei]/pp[facei];
				pmu[facei] 		= pmuFavre[facei];
				palpha[facei] 	= palphaFavre[facei];
				pdefect[facei]  = pH[facei] - (HOxidizer+pZ[facei]*(HFuel-HOxidizer));
			}
		}
	}
}

template<class BasicFlameletSMOKEThermo, class MixtureType>
void Foam::pdfFlameletSMOKEThermo<BasicFlameletSMOKEThermo, MixtureType>::update()
{
	std::vector<double> extracted(7);
	double Zvar_normalized = 0.;
	double defect = 0.;


	const scalarField& Z 		= this->Z_.internalField();
	const scalarField& Zvar 	= this->Zvar_.internalField();
	const scalarField& chi_st 	= this->chi_st_.internalField();
	const scalarField& HCells 	= this->H_.internalField();

	scalarField& TCells 		= this->T_.internalField();
	scalarField& RhoCells 		= this->density_reynolds_.internalField();
	scalarField& asCells 		= this->as_.internalField();
	scalarField& muCells 		= this->mu_favre_.internalField();
	scalarField& alphaCells 	= this->alpha_favre_.internalField();

	double small_eps = 1.e-6;
	double small_chi_st = 1.e-8;

	//- Internal cells
	forAll(Z, celli)
	{
		double max_chi = max(small_chi_st,chi_st[celli]);

		if (adiabaticMode == false)
			defect = HCells[celli] - (HOxidizer+Z[celli]*(HFuel-HOxidizer));

		//- Pure oxidizer
		if (Z[celli]<=small_eps)
		{
			flamelets_library.GetMeanValues(0., 0., max_chi, defect, extracted);
		}

		//- Pure fuel
		else if (Z[celli]>=(1.-small_eps))
		{
			flamelets_library.GetMeanValues(1., 0., max_chi, defect, extracted);
		}

		//- Mixture
		else
		{
			Zvar_normalized = Zvar[celli] / (Z[celli]*(1.-Z[celli]));
			if (Zvar_normalized >= 0.98)
				flamelets_library.GetMeanValues(Z[celli], 0.98, max_chi, defect, extracted);
			else if (Zvar_normalized < 0.)
				flamelets_library.GetMeanValues(Z[celli], 0.00, max_chi, defect, extracted);
			else
				flamelets_library.GetMeanValues(Z[celli], Zvar_normalized, max_chi, defect, extracted);
		}

		TCells[celli] 		= extracted[1];
		RhoCells[celli]		= extracted[2];
		asCells[celli]		= extracted[3];
		muCells[celli]		= extracted[4];
		alphaCells[celli]	= extracted[5];
	}

	//- Boundary conditions
	if (adiabaticMode == true)
	{
		forAll(Z_.boundaryField(), patchi)
		{
			const fvPatchScalarField& pcsi 		= this->Z_.boundaryField()[patchi];
			const fvPatchScalarField& pcsiv2 	= this->Zvar_.boundaryField()[patchi];
			const fvPatchScalarField& pchi_st 	= this->chi_st_.boundaryField()[patchi];

			fvPatchScalarField& pt		= this->T_.boundaryField()[patchi];
			fvPatchScalarField& prho 	= this->density_reynolds_.boundaryField()[patchi];
			fvPatchScalarField& pas 	= this->as_.boundaryField()[patchi];
			fvPatchScalarField& pmu 	= this->mu_favre_.boundaryField()[patchi];
			fvPatchScalarField& palpha 	= this->alpha_favre_.boundaryField()[patchi];

			forAll(pcsi, facei)
			{

				double max_chi = max(small_chi_st, pchi_st[facei]);

				//- Pure oxidizer
				if (pcsi[facei]<=0.)
				{
					flamelets_library.GetMeanValues(0., 0., max_chi, defect, extracted);
				}

				//- Pure fuel
				else if (pcsi[facei]>=1.)
				{
					flamelets_library.GetMeanValues(1., 0., max_chi, defect, extracted);
				}

				//- Mixture
				else
				{
					Zvar_normalized = pcsiv2[facei] / (pcsi[facei]*(1.-pcsi[facei]));

					if (Zvar_normalized >= 0.98)
						flamelets_library.GetMeanValues(pcsi[facei], 0.98, max_chi, defect, extracted);
					else if (Zvar_normalized < 0.)
						flamelets_library.GetMeanValues(pcsi[facei], 0.00, max_chi, defect, extracted);
					else
						flamelets_library.GetMeanValues(pcsi[facei], Zvar_normalized, max_chi,defect, extracted);
				}

				pt[facei] 	= extracted[1];
				prho[facei]	= extracted[2];
				pas[facei]	= extracted[3];
				pmu[facei]	= extracted[4];
				palpha[facei]	= extracted[5];

			}
		}
	}
	else
	{
		forAll(Z_.boundaryField(), patchi)
		{

			if (patch_type_T[patchi] == 0)
			{
				const fvPatchScalarField& pcsi 		= this->Z_.boundaryField()[patchi];
				const fvPatchScalarField& pcsiv2 	= this->Zvar_.boundaryField()[patchi];
				const fvPatchScalarField& pchi_st 	= this->chi_st_.boundaryField()[patchi];
				const fvPatchScalarField& ph		= this->H_.boundaryField()[patchi];

				fvPatchScalarField& pt		= this->T_.boundaryField()[patchi];
				fvPatchScalarField& prho 	= this->density_reynolds_.boundaryField()[patchi];
				fvPatchScalarField& pas 	= this->as_.boundaryField()[patchi];
				fvPatchScalarField& pmu 	= this->mu_favre_.boundaryField()[patchi];
				fvPatchScalarField& palpha 	= this->alpha_favre_.boundaryField()[patchi];

				forAll(pcsi, facei)
				{

					double max_chi = max(small_chi_st, pchi_st[facei]);

					defect = ph[facei] - (HOxidizer+pcsi[facei]*(HFuel-HOxidizer));

					//- Pure oxidizer
					if (pcsi[facei]<=0.)
					{
						flamelets_library.GetMeanValues(0., 0., max_chi, defect, extracted);
					}

					//- Pure fuel
					else if (pcsi[facei]>=1.)
					{
						flamelets_library.GetMeanValues(1., 0., max_chi, defect, extracted);
					}

					//- Mixture
					else
					{
						Zvar_normalized = pcsiv2[facei] / (pcsi[facei]*(1.-pcsi[facei]));

						if (Zvar_normalized >= 0.98)
							flamelets_library.GetMeanValues(pcsi[facei], 0.98, max_chi, defect, extracted);
						else if (Zvar_normalized < 0.)
							flamelets_library.GetMeanValues(pcsi[facei], 0.00, max_chi, defect, extracted);
						else
							flamelets_library.GetMeanValues(pcsi[facei], Zvar_normalized, max_chi,defect, extracted);
					}

					pt[facei] 	= extracted[1];
					prho[facei]	= extracted[2];
					pas[facei]	= extracted[3];
					pmu[facei]	= extracted[4];
					palpha[facei]	= extracted[5];
				}
			}

			//- Added for friendly handling
			//- fixed temperature and fixed enhalpy :: INLETS
			else if (patch_type_T[patchi] == 1 && patch_type_H[patchi] == 1 && patch_type_Z[patchi] == 1)
			{
				const fvPatchScalarField& pcsi 		= this->Z_.boundaryField()[patchi];
				const fvPatchScalarField& pcsiv2 	= this->Zvar_.boundaryField()[patchi];
				const fvPatchScalarField& pchi_st 	= this->chi_st_.boundaryField()[patchi];

				fvPatchScalarField& prho 		= this->density_reynolds_.boundaryField()[patchi];
				fvPatchScalarField& pas 		= this->as_.boundaryField()[patchi];
				fvPatchScalarField& pmu 		= this->mu_favre_.boundaryField()[patchi];
				fvPatchScalarField& palpha 		= this->alpha_favre_.boundaryField()[patchi];

				forAll(pcsi, facei)
				{
					defect = 0.;

					double max_chi = max(small_chi_st, pchi_st[facei]);

					//- Pure oxidizer
					if (pcsi[facei]<=0.)
					{
						flamelets_library.GetMeanValues(0., 0., max_chi, defect, extracted);
					}

					//- Pure fuel
					else if (pcsi[facei]>=1.)
					{
						flamelets_library.GetMeanValues(1., 0., max_chi, defect, extracted);
					}

					//- Mixture
					else
					{
						Zvar_normalized = pcsiv2[facei] / (pcsi[facei]*(1.-pcsi[facei]));

						if (Zvar_normalized >= 0.98)
						{
							flamelets_library.GetMeanValues(pcsi[facei], 0.98, max_chi, defect, extracted);
						}
						else if (Zvar_normalized < 0.)
						{
							flamelets_library.GetMeanValues(pcsi[facei], 0.00, max_chi, defect, extracted);
						}
						else
						{
							flamelets_library.GetMeanValues(pcsi[facei], Zvar_normalized, max_chi,defect, extracted);
						}
					}

					prho[facei]	= extracted[2];
					pas[facei]	= extracted[3];
					pmu[facei]	= extracted[4];
					palpha[facei]	= extracted[5];
				}
			}

			//- Added for enthalpydefect due to fixedTemperature on walls
			//- fixed temperature and fixedEnthalpy :: no inlet
			else if (patch_type_T[patchi] == 1 && patch_type_H[patchi] == 1 && patch_type_Z[patchi] == 0)
			{

				const fvPatchScalarField& pcsi 		= this->Z_.boundaryField()[patchi];
				const fvPatchScalarField& pcsiv2 	= this->Zvar_.boundaryField()[patchi];
				const fvPatchScalarField& pchi_st 	= this->chi_st_.boundaryField()[patchi];

				fvPatchScalarField& ph			= this->H_.boundaryField()[patchi];
				fvPatchScalarField& pt			= this->T_.boundaryField()[patchi];
				fvPatchScalarField& prho 		= this->density_reynolds_.boundaryField()[patchi];
				fvPatchScalarField& pas 		= this->as_.boundaryField()[patchi];
				fvPatchScalarField& pmu 		= this->mu_favre_.boundaryField()[patchi];
				fvPatchScalarField& palpha 		= this->alpha_favre_.boundaryField()[patchi];

				forAll(pcsi, facei)
				{

					double max_chi = max(small_chi_st, pchi_st[facei]);

					//- Pure oxidizer
					if (pcsi[facei]<=0.)
					{
						defect = flamelets_library.GetEnthalpyDefectFromTemperature(0., 0., max_chi, pt[facei]);
						ph[facei] = defect + HOxidizer;
						flamelets_library.GetMeanValues(0., 0., max_chi, defect, extracted);
					}

					//- Pure fuel
					else if (pcsi[facei]>=1.)
					{
						defect = flamelets_library.GetEnthalpyDefectFromTemperature(1., 0., max_chi, pt[facei]);
						ph[facei] = defect + HFuel;
						flamelets_library.GetMeanValues(1., 0., max_chi, defect, extracted);
					}

					//- Mixture
					else
					{
						Zvar_normalized = pcsiv2[facei] / (pcsi[facei]*(1.-pcsi[facei]));

						if (Zvar_normalized >= 0.98)
						{
							defect = flamelets_library.GetEnthalpyDefectFromTemperature(pcsi[facei], 0.98, max_chi, pt[facei]);
							ph[facei] = defect + (HOxidizer+pcsi[facei]*(HFuel-HOxidizer));
							flamelets_library.GetMeanValues(pcsi[facei], 0.98, max_chi, defect, extracted);
						}
						else if (Zvar_normalized < 0.)
						{
							defect = flamelets_library.GetEnthalpyDefectFromTemperature(pcsi[facei], 0.00, max_chi, pt[facei]);
							ph[facei] = defect + (HOxidizer+pcsi[facei]*(HFuel-HOxidizer));
							flamelets_library.GetMeanValues(pcsi[facei], 0.00, max_chi, defect, extracted);
						}
						else
						{
							defect = flamelets_library.GetEnthalpyDefectFromTemperature(pcsi[facei], Zvar_normalized, max_chi, pt[facei]);
							ph[facei] = defect + (HOxidizer+pcsi[facei]*(HFuel-HOxidizer));
							flamelets_library.GetMeanValues(pcsi[facei], Zvar_normalized, max_chi,defect, extracted);
						}
					}

					prho[facei]	= extracted[2];
					pas[facei]	= extracted[3];
					pmu[facei]	= extracted[4];
					palpha[facei]	= extracted[5];
				}
			}

			//- Not implemented boundary condition :: give an error
			else
			{
				FatalErrorIn
				(
					"pdfFlameletSMOKEThermo<BasicFlameletSMOKEThermo, MixtureType>::update()"
				)
				<< "Boundary conditions are wrong: "
				<< "fixed temperature BC must be fixed enthaplie BC with value 0;"
				<< abort(FatalError);
			}
		}
	}
}

template<class BasicFlameletSMOKEThermo, class MixtureType>
void Foam::pdfFlameletSMOKEThermo<BasicFlameletSMOKEThermo, MixtureType>::updateMassFractions()
{
	std::vector<double> extracted(flamelets_library.number_of_species()+1);

	double Zvar_normalized = 0.;
	double defect = 0.;

	const scalarField& Z 		= this->Z_.internalField();
	const scalarField& Zvar		= this->Zvar_.internalField();
	const scalarField& chi_st 	= this->chi_st_.internalField();
	const scalarField& HCells 	= this->H_.internalField();

	double small_eps = 1.e-6;
	double small_chi_st = 1.e-8;

	forAll(Z, celli)
	{
		double max_chi = max(small_chi_st,chi_st[celli]);

		if (adiabaticMode == false)
		defect = HCells[celli] - (HOxidizer+Z[celli]*(HFuel-HOxidizer));

		//- Pure oxidizer
		if (Z[celli]<=small_eps)
		{
			flamelets_library.ExtractMeanValues(0., 0., max_chi, defect, extracted);
		}

		//- Pure fuel
		else if (Z[celli]>=(1.-small_eps))
		{
			flamelets_library.ExtractMeanValues(1., 0., max_chi, defect, extracted);
		}

		//- Mixture
		else
		{
			Zvar_normalized = Zvar[celli] / (Z[celli]*(1.-Z[celli]));

			if (Zvar_normalized >= 0.98)
				flamelets_library.ExtractMeanValues(Z[celli], 0.98, max_chi, defect, extracted);
			else if (Zvar_normalized < 0.)
				flamelets_library.ExtractMeanValues(Z[celli], 0.00, max_chi, defect, extracted);
			else
				flamelets_library.ExtractMeanValues(Z[celli], Zvar_normalized, max_chi, defect, extracted);
		}

		for(int j=0;j<flamelets_library.number_of_species()+1;j++)
		{
			if(j<flamelets_library.number_of_species())
			{
				omega_[j].internalField()[celli] = extracted[j+1];
			}
		}
	}

	forAll(Z_.boundaryField(), patchi)
	{
		const fvPatchScalarField& pZ 		= this->Z_.boundaryField()[patchi];
		const fvPatchScalarField& pZvar 	= this->Zvar_.boundaryField()[patchi];
		const fvPatchScalarField& pchi_st 	= this->chi_st_.boundaryField()[patchi];
		const fvPatchScalarField& ph		= this->H_.boundaryField()[patchi];

		forAll(pZ, facei)
		{

			double max_chi = max(small_chi_st, pchi_st[facei]);

			if (adiabaticMode == false)
				defect = ph[facei] - (HOxidizer+pZ[facei]*(HFuel-HOxidizer));

			//- Pure oxidizer
			if (pZ[facei]<=small_eps)
			{
				flamelets_library.ExtractMeanValues(0., 0., max_chi, defect, extracted);
			}

			//- Pure fuel
			else if (pZ[facei]>=(1.-small_eps))
			{
				flamelets_library.ExtractMeanValues(1., 0., max_chi, defect, extracted);
			}

			//- Mixture
			else
			{
				Zvar_normalized = pZvar[facei] / (pZ[facei]*(1.-pZ[facei]));

				if (Zvar_normalized >= 0.98)
					flamelets_library.ExtractMeanValues(pZ[facei], 0.98, max_chi, defect, extracted);
				else if (Zvar_normalized < 0.)
					flamelets_library.ExtractMeanValues(pZ[facei], 0.00, max_chi, defect, extracted);
				else
					flamelets_library.ExtractMeanValues(pZ[facei], Zvar_normalized, max_chi,defect, extracted);
			}

			for(int j=0;j<flamelets_library.number_of_species()+1;j++)
			{
				if(j<flamelets_library.number_of_species())
				{
					omega_[j].boundaryField()[patchi][facei] = extracted[j+1];
				}
			}
		}
	}
}

template<class BasicFlameletSMOKEThermo, class MixtureType>
void Foam::pdfFlameletSMOKEThermo<BasicFlameletSMOKEThermo, MixtureType>::errorMessage(const string message)
{
	Info << "Class: pdfFlameletSMOKEThermo" << endl;
	Info << "Error: " << message << endl;
	getchar();
}

template<class BasicFlameletSMOKEThermo, class MixtureType>
void Foam::pdfFlameletSMOKEThermo<BasicFlameletSMOKEThermo, MixtureType>::infoMessage() const
{
	Info << "/*-------------------------------------*\\" << endl <<
			"|    Flamelet Thermo initialization     |" << endl <<
			"|                                       |" << endl <<
			"|   Rebuild by Tobias Holzmann M.Eng.   |" << endl <<
			"|          www.Holzmann-cfd.de          |" << endl <<
			"\\*-------------------------------------*/" << endl << endl;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicFlameletSMOKEThermo, class MixtureType>
Foam::pdfFlameletSMOKEThermo<BasicFlameletSMOKEThermo, MixtureType>::pdfFlameletSMOKEThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    heThermo<BasicFlameletSMOKEThermo, MixtureType>(mesh, phaseName),

    Z_
    (
    	IOobject
    	(
    		"Z",
    		mesh.time().timeName(),
    	    mesh,
    	    IOobject::MUST_READ,
    	    IOobject::AUTO_WRITE
    	),
    	mesh
    ),

    Zvar_
    (
      	IOobject
       	(
       		"Zvar",
       		mesh.time().timeName(),
       	    mesh,
       	    IOobject::MUST_READ,
       	    IOobject::AUTO_WRITE
       	),
       	mesh
    ),

    chi_st_
    (
      	IOobject
       	(
       		"chi_st",
       		mesh.time().timeName(),
       	    mesh,
       	    IOobject::NO_READ,
       	    IOobject::AUTO_WRITE
       	),
       	mesh,
       	dimensionedScalar("chi_st",dimensionSet(0,0,-1,0,0,0,0) , 0.0)
    ),

    H_
    (
      	IOobject
       	(
       		"H",
       		mesh.time().timeName(),
       	    mesh,
       	    IOobject::MUST_READ,
       	    IOobject::AUTO_WRITE
       	),
       	mesh
    ),

    defect_
    (
      	IOobject
       	(
       		"defect",
       		mesh.time().timeName(),
       	    mesh,
       	    IOobject::NO_READ,
       	    IOobject::AUTO_WRITE
       	),
       	mesh,
       	dimensionedScalar("defect",dimensionSet(0,2,-2,0,0,0,0) , 0.0)
    ),

    as_
    (
      	IOobject
       	(
       		"as",
       		mesh.time().timeName(),
       	    mesh,
       	    IOobject::NO_READ,
       	    IOobject::AUTO_WRITE
       	),
       	mesh,
       	dimensionedScalar("as",dimensionSet(0,-1,0,0,0,0,0) , 0.0)
    ),

    density_reynolds_
    (
      	IOobject
       	(
       		"rho_reynolds",
       		mesh.time().timeName(),
       	    mesh,
       	    IOobject::NO_READ,
       	    IOobject::NO_WRITE
       	),
       	mesh,
       	dimensionedScalar("density_reynolds",dimensionSet(1,-3,0,0,0,0,0) , 0.0)
    ),

    mu_favre_
    (
      	IOobject
       	(
       		"mu_lam",
       		mesh.time().timeName(),
       	    mesh,
       	    IOobject::NO_READ,
       	    IOobject::AUTO_WRITE
       	),
       	mesh,
       	dimensionedScalar("mu_lam",dimensionSet(1,-1,-1,0,0,0,0) , 0.0)
    ),

    alpha_favre_
    (
      	IOobject
       	(
       		"alpha_lam",
       		mesh.time().timeName(),
       	    mesh,
       	    IOobject::NO_READ,
       	    IOobject::AUTO_WRITE
       	),
       	mesh,
       	dimensionedScalar("alpha_lam",dimensionSet(0,2,-1,0,0,0,0) , 0.0)
    ),

    adiabaticMode(false),
    showFlamelet(false),
    showFlameletLibrary(false)

{
	//- Info message
	infoMessage();

	//- Get fixedValue enthalpy boundaries
	Info << "Enthalpy:" << endl;
	forAll(this->H_.boundaryField(), patchi)
	{
		if (isA<fixedValueFvPatchScalarField>(this->H_.boundaryField()[patchi]))
		{
			Info << "     + " << mesh.boundary()[patchi].name() << " <fixedValue>" << endl;
			patch_type_H.push_back(1);
		}
		else
		{
			Info << "     + " << mesh.boundary()[patchi].name() << "" << endl;
			patch_type_H.push_back(0);
		}
	}

	//- Get fixedValue temperatur boundaries
	Info << endl << "Temperature: " << endl;
	forAll(this->T_.boundaryField(), patchi)
	{
		if (isA<fixedValueFvPatchScalarField>(this->T_.boundaryField()[patchi]))
		{
			Info << "     + " << mesh.boundary()[patchi].name() << " <fixedValue>" << endl;
			patch_type_T.push_back(1);
		}
		else
		{
			Info << "     + " << mesh.boundary()[patchi].name() << "" << endl;
			patch_type_T.push_back(0);
		}
	}

	//- Get fuel and oxidizer inlets
	Info << endl << "Mixture fraction:" << endl;
	forAll(this->Z_.boundaryField(), patchi)
	{
		if (isA<fixedValueFvPatchScalarField>(this->Z_.boundaryField()[patchi]))
		{
			Info << "     + " << mesh.boundary()[patchi].name() << " <fixedValue>" << endl;
			patch_type_Z.push_back(1);
		}
		else
		{
			Info << "     + " << mesh.boundary()[patchi].name() << "" << endl;
			patch_type_Z.push_back(0);
		}
	}
	Info << endl;

	//- IOFlamelet properties
	IOdictionary flameletsProperties_
	(
		IOobject
		(
		    "flameletsProperties",
		    Z_.time().constant(),
		    Z_.db(),
		    IOobject::MUST_READ,
		    IOobject::NO_WRITE
		)
	);

	//- Get flamelet properties and set fields and variables
	{
		Switch adiabaticMode_(flameletsProperties_.lookup("adiabaticMode"));
		Switch showFlamelet_(flameletsProperties_.lookup("showFlamelet"));
		Switch showFlameletLibrary_(flameletsProperties_.lookup("showFlameletLibrary"));
		propertyUpdate	 	= readLabel(flameletsProperties_.lookup("propertyUpdate"));
		massFractionsUpdate 	= readLabel(flameletsProperties_.lookup("massFractionsUpdate"));

		adiabaticMode 		= adiabaticMode_;
		showFlamelet		= showFlamelet_;
		showFlameletLibrary	= showFlameletLibrary_;
		counter 		= propertyUpdate;
		counter_mass_fractions 	= massFractionsUpdate;

		string libraryPath 	= flameletsProperties_.lookup("libraryPath");
		string chiPDF 		= flameletsProperties_.lookup("pdf");
		string list_of_species 	= flameletsProperties_.lookup("species");

		boost::filesystem::path library_folder = libraryPath;
		flamelets_library.SetLibraryFolder(library_folder);
		flamelets_library.SetSpeciesToExtract(list_of_species);

		//- Set the dissipation mode

			//- Set adiabatic mode
			if (adiabaticMode == true)
			{
				Info << "Flamelet thermo mode is <adiabatic>" << endl;
				flamelets_library.SetAdiabaticMode();
			}
			else
			{
				Info << "Flamelet thermo mode is <non-adiabatic>" << endl;
			}

			//- Scalar dissipation rate distribution
			if (chiPDF == "logNormal")
			{
				Info << "Flamelet thermo use the <log-normal distribution> for the scalar dissipation rate" << endl;

				scalar chi_sigma(readScalar(flameletsProperties_.lookup("sigma")));
				flamelets_library.SetLogNormalChiDistribution(chi_sigma);
			}
			else if (chiPDF == "dirac")
			{
				Info << "Flamelet thermo uses the <delta-dirac distribution> for the scalar dissipation rate" << endl;
			}

			//- Set show flamelet mode
			if (showFlamelet == true)
			{
				Info << "Flamelet thermo shows all single flamelet properties" << endl;
				flamelets_library.SetShowFlamelet();
			}
			else
			{
				Info << "Flamelet thermo does not show the single flamelet properties" << endl;
			}

			//- Set show flamelet library mode
			if (showFlameletLibrary == true)
			{
				Info << "Flamelet thermo shows all flamelet library properties" << endl;
				flamelets_library.SetShowFlameletLibrary();
			}
			else
			{
				Info << "Flamelet thermo does not show the flamelet library properties" << endl;
			}

		//- Initialise the wanted mass fraction fields
		omega_.setSize(flamelets_library.number_of_species());

		for (int j=0;j<flamelets_library.number_of_species();j++)
		{
			if(j < flamelets_library.number_of_species())
			{
				std::string name_of_species = "omega_" + flamelets_library.species()[j];

				omega_.set
				(
					j,
					new volScalarField
					(
						IOobject
						(
							name_of_species,
							mesh.time().timeName(),
							mesh,
							IOobject::NO_READ,
							IOobject::AUTO_WRITE),
							mesh,
							dimensionedScalar("zero", dimensionSet(0,0,0,0,0,0,0),
							0.0
						)
					)
				);
			}
		}
	}

	//- Read all flamelet
	flamelets_library.Read();
	flamelets_library.Summary();

	//- Set the adiabat enthalpy
	Info << "Thermo flamelet set the adiabatic enthalpy for enthalpy defect calculation:" << endl << endl;
	{
		HOxidizer = flamelets_library.enthalpy_f_oxidizer();
		HFuel = flamelets_library.enthalpy_f_fuel();

		Info << "     + Adiabatic enthalpy fuel:     " << HFuel << endl;
		Info << "     + Adiabatic enthalpy oxidizer: " << HOxidizer << endl << endl;
	}

	//- Extract all variables
	update();

	//- Calculate first time
    	calculate();

    // Switch on saving old time
    this->psi_.oldTime();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicFlameletSMOKEThermo, class MixtureType>
Foam::pdfFlameletSMOKEThermo<BasicFlameletSMOKEThermo, MixtureType>::~pdfFlameletSMOKEThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class BasicFlameletSMOKEThermo, class MixtureType>
void Foam::pdfFlameletSMOKEThermo<BasicFlameletSMOKEThermo, MixtureType>::correct()
{
    if (debug)
    {
        Info<< "entering pdfFlameletSMOKEThermo<BasicFlameletSMOKEThermo, MixtureType>::correct()"
            << endl;
    }

    // force the saving of the old-time values
    this->psi_.oldTime();

    if(counter == propertyUpdate || counter_mass_fractions == massFractionsUpdate)
    {
    	Info << endl << "Flamelet thermo: " << endl;
    }

	if (counter == propertyUpdate)
	{
		Info << "     + Flamelet thermo updates all thermo variables from the Look-Up-Table" << endl;

		update();
		counter = 0;
	}

	if (counter_mass_fractions == massFractionsUpdate)
	{
		Info << "     + Flamelet thermo updates all mass fraction from Look-Up-Table" << endl << endl;

		updateMassFractions();
		counter_mass_fractions = 0;
	}
	else if(counter_mass_fractions != massFractionsUpdate && counter == 0)
	{
		Info << endl;
	}
	else
	{
	}

    calculate();

    counter++;
    counter_mass_fractions++;

    if (debug)
    {
        Info<< "exiting pdfFlameletSMOKEThermo<BasicFlameletSMOKEThermo, MixtureType>::correct()"
            << endl;
    }
}


// ************************************************************************* //

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

Class
    Foam::pdfFlameletSMOKEThermo

Description
    Energy for a mixture based on compressibility

SourceFiles
    pdfFlameletSMOKEThermo.C

\*---------------------------------------------------------------------------*/

#ifndef pdfFlameletSMOKEThermo_H
#define pdfFlameletSMOKEThermo_H

#include "flameletSMOKEThermo.H"
#include "heThermo.H"
#include "flames/flamelets/lookuptablereader/NonAdiabaticFlameletLibraryReader.h"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class pdfFlameletSMOKEThermo Declaration
\*---------------------------------------------------------------------------*/

template<class BasicFlameletSMOKEThermo, class MixtureType>
class pdfFlameletSMOKEThermo
:
    public heThermo<BasicFlameletSMOKEThermo, MixtureType>
{

    // Private Member Functions

        //- Calculate the thermo variables (get it from Look-Up-Table)
        void calculate();

        //- Update the thermo variables by extracting from Look-Up-Table
        void update();

        //- Update the species mass fraction by extracting from Look-Up-Table
    	void updateMassFractions();

    	//- Errormessage
    	void errorMessage(const string message);

    	//- Infomessage
    	void infoMessage() const;

        //- Construct as copy (not implemented)
        pdfFlameletSMOKEThermo(const pdfFlameletSMOKEThermo<BasicFlameletSMOKEThermo, MixtureType>&);

public:

    //- Runtime type information
    TypeName("pdfFlameletSMOKEThermo");


    // Constructors

        //- Construct from mesh and phase name
        pdfFlameletSMOKEThermo
        (
            const fvMesh&,
            const word& phaseName
        );


    //- Destructor
    virtual ~pdfFlameletSMOKEThermo();


    // Member functions

        //- Update properties
        virtual void correct();

private:

    //- Variables

        //- Mixture fraction
        volScalarField Z_;

        //- Variance of mixture fraction
        volScalarField Zvar_;

        //- Scalar dissipation rate
        volScalarField chi_st_;

        //- Enthalpy of the mixture
        volScalarField H_;

        //- Enthalpydefect due to radiation and BC (fixed wall temperature)
        volScalarField defect_;

        //- Plank absorbtions coefficient (for radiation model in free stream flames)
        volScalarField as_;

        // - Reynolds density [kg/m³]
        volScalarField density_reynolds_;

        //- Dynamic viscosity (favre)
        volScalarField mu_favre_;

        //- Thermal diffusivity
        volScalarField alpha_favre_;

        //- Species mass fraction
        PtrList<volScalarField> omega_;

        //- adiabat enthalpy Fuel
        double HFuel;

        //- adiabat enthalpy Oxidizer
        double HOxidizer;

        //- variable for set thermo variable extraction
        label counter;

        //- variable for set mass fraction extraction
        label counter_mass_fractions;

        //- User defined numbers for update
        label propertyUpdate;
        label massFractionsUpdate;


        label nDirac;
        Switch adiabaticMode;
        Switch showFlamelet;
        Switch showFlameletLibrary;

        //- stuff (modifyed by Holzmann)
    	std::vector<int> patch_type_T;
    	std::vector<int> patch_type_H;
    	std::vector<int> patch_type_Z;

private:

    //- Non adiabatic flamelet library
    OpenSMOKE::NonAdiabaticFlameletLibraryReader flamelets_library;


public:

	//- Mixture fraction

		//  Non-const access allowed for transport equations
		virtual volScalarField& Z()
		{
			return Z_;
		}

		//- Mixture fraction
		virtual const volScalarField& Z() const
		{
			return Z_;
		}

		//- Variance of Mixture fraction
		//  Non-const access allowed for transport equations
		virtual volScalarField& Zvar()
		{
			return Zvar_;
		}

		//- Variance of Mixture fraction
		virtual const volScalarField& Zvar() const
		{
			return Zvar_;
		}

		//- Scalar dissipation rate [1/s]
		//  Non-const access allowed for transport equations
		virtual volScalarField& chi_st()
		{
			return chi_st_;
		}

		//- Scalar dissipation rate [1/s]
		virtual const volScalarField& chi_st() const
		{
			return chi_st_;
		}

		//- Enthalpy [J/kg]
		//  Non-const access allowed for transport equations
		virtual volScalarField& H()
		{
			return H_;
		}

		//- Enthalpy [J/kg]
		virtual const volScalarField& H() const
		{
			return H_;
		}

		//- Absorption coefficient [1/m]
		//  Non-const access allowed for transport equations
		virtual volScalarField& as()
		{
			return as_;
		}

		//- Absorption coefficient [1/m]
		virtual const volScalarField& as() const
		{
			return as_;
		}

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
#   include "pdfFlameletSMOKEThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

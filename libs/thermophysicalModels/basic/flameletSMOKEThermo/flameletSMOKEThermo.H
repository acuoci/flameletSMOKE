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
    Foam::flameletSMOKEThermo

Description
    Basic thermodynamic properties based on compressibility

SourceFiles
    flameletSMOKEThermo.C

\*---------------------------------------------------------------------------*/

#ifndef flameletSMOKEThermo_H
#define flameletSMOKEThermo_H

#include "fluidThermo.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class flameletSMOKEThermo Declaration
\*---------------------------------------------------------------------------*/

class flameletSMOKEThermo
:
    public fluidThermo
{

protected:

    // Protected data

        //- Compressibility [s^2/m^2]
        volScalarField psi_;

        //- Dynamic viscosity [kg/m/s]
        volScalarField mu_;


    // Protected Member Functions

        //- Construct as copy (not implemented)
        flameletSMOKEThermo(const flameletSMOKEThermo&);


public:

    //- Runtime type information
    TypeName("flameletSMOKEThermo");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        flameletSMOKEThermo,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );


    // Constructors

        //- Construct from mesh and phase name
        flameletSMOKEThermo
        (
            const fvMesh&,
            const word& phaseName
        );


    //- Selector
    static autoPtr<flameletSMOKEThermo> New
    (
        const fvMesh& mesh,
        const word& phaseName=word::null
    );


    //- Destructor
    virtual ~flameletSMOKEThermo();

    //- Member functions

        //- Virtual stuff
        virtual volScalarField& Z();

        virtual const volScalarField& Z() const;

        virtual volScalarField& Zvar();

        virtual const volScalarField& Zvar() const;

        virtual volScalarField& chi_st();

        virtual const volScalarField& chi_st() const;

        virtual volScalarField& H();

        virtual const volScalarField& H() const;

        virtual volScalarField& as();

        virtual const volScalarField& as() const;

        // Fields derived from thermodynamic state variables

        //- Add the given density correction to the density field.
        //  Used to update the density field following pressure solution.
        //  For psiThermo does nothing.
        virtual void correctRho(const volScalarField& deltaRho);


        //- Density [kg/m^3] - uses current value of pressure
        virtual tmp<volScalarField> rho() const;
	virtual tmp<scalarField> rho(const label patchi) const;


        //- Compressibility [s^2/m^2]
        virtual const volScalarField& psi() const;


        // Access to transport state variables

	#if OPENFOAM_VERSION < 40
        //- Dynamic viscosity of mixture [kg/m/s]
        virtual const volScalarField& mu() const;
	virtual const scalarField& mu(const label patchi) const;
	#else
	virtual tmp<volScalarField> mu() const;
	virtual tmp<scalarField> mu(const label patchi) const;
	#endif
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

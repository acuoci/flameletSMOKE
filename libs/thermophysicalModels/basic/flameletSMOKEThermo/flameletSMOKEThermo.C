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
    (at your option) any later version.flam

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "flameletSMOKEThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(flameletSMOKEThermo, 0);
    defineRunTimeSelectionTable(flameletSMOKEThermo, fvMesh);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::flameletSMOKEThermo::flameletSMOKEThermo(const fvMesh& mesh, const word& phaseName)
:
    fluidThermo(mesh, phaseName),

    psi_
    (
        IOobject
        (
            phasePropertyName("thermo:psi"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(0, -2, 2, 0, 0)
    ),

    mu_
    (
        IOobject
        (
            phasePropertyName("thermo:mu"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(1, -1, -1, 0, 0)
    )
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::flameletSMOKEThermo> Foam::flameletSMOKEThermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    return basicThermo::New<flameletSMOKEThermo>(mesh, phaseName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::flameletSMOKEThermo::~flameletSMOKEThermo()
{}


// * * * * * * * * * * * * * * * Virtual Stuff  * * * * * * * * * * * * * * * //

Foam::volScalarField& Foam::flameletSMOKEThermo::Z()
{
    notImplemented("flameletSMOKEThermo::Z()");
    return const_cast<volScalarField&>(volScalarField::null());
}

const Foam::volScalarField& Foam::flameletSMOKEThermo::Z() const
{
    notImplemented("flameletSMOKEThermo::Z() const");
    return volScalarField::null();
}

Foam::volScalarField& Foam::flameletSMOKEThermo::Zvar()
{
    notImplemented("flameletSMOKEThermo::Zvar()");
    return const_cast<volScalarField&>(volScalarField::null());
}

const Foam::volScalarField& Foam::flameletSMOKEThermo::Zvar() const
{
    notImplemented("flameletSMOKEThermo::Zvar() const");
    return volScalarField::null();
}

Foam::volScalarField& Foam::flameletSMOKEThermo::chi_st()
{
    notImplemented("flameletSMOKEThermo::chi_st()");
    return const_cast<volScalarField&>(volScalarField::null());
}

const Foam::volScalarField& Foam::flameletSMOKEThermo::chi_st() const
{
    notImplemented("flameletSMOKEThermo::chi_st() const");
    return volScalarField::null();
}

Foam::volScalarField& Foam::flameletSMOKEThermo::H()
{
    notImplemented("flameletSMOKEThermo::H()");
    return const_cast<volScalarField&>(volScalarField::null());
}

const Foam::volScalarField& Foam::flameletSMOKEThermo::H() const
{
    notImplemented("flameletSMOKEThermo::H() const");
    return volScalarField::null();
}

Foam::volScalarField& Foam::flameletSMOKEThermo::as()
{
    notImplemented("flameletSMOKEThermo::as()");
    return const_cast<volScalarField&>(volScalarField::null());
}

const Foam::volScalarField& Foam::flameletSMOKEThermo::as() const
{
    notImplemented("flameletSMOKEThermo::as() const");
    return volScalarField::null();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::flameletSMOKEThermo::rho() const
{
    return p_*psi_;
}

Foam::tmp<Foam::scalarField> Foam::flameletSMOKEThermo::rho(const label patchi) const
{
    return p_.boundaryField()[patchi]*psi_.boundaryField()[patchi];
}

const Foam::volScalarField& Foam::flameletSMOKEThermo::psi() const
{
    return psi_;
}

const Foam::volScalarField& Foam::flameletSMOKEThermo::mu() const
{
    return mu_;
}

const Foam::scalarField& Foam::flameletSMOKEThermo::mu(const label patchi) const
{
    return mu_.boundaryField()[patchi];
}	

// ************************************************************************* //

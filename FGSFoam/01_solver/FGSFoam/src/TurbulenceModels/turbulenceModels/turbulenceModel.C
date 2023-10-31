/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "turbulenceModel.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(turbulenceModel, 0);
}

const Foam::word Foam::turbulenceModel::propertiesName("turbulenceProperties");


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turbulenceModel::turbulenceModel
(
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const word& propertiesName
)
:
    IOdictionary
    (
        IOobject
        (
            IOobject::groupName(propertiesName, alphaRhoPhi.group()),
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    runTime_(U.time()),
    mesh_(U.mesh()),

    U_(U),
    alphaRhoPhi_(alphaRhoPhi),
    phi_(phi),
    y_(mesh_),

        Z_
    (
        IOobject
        (
            "Z",
            U.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    /*
    varZ_
    (
        IOobject
        (
            "varZ",
            U.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    */
    PV_                                 //--- Added L.Ma, 07-10-2014
    (
        IOobject
        (
            "PV",
            U.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    /*
    varPV_
    (                                   //--- Added L.Ma, 07-10-2014
        IOobject
        (
            "varPV",
            U.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    */
    sourcePV_                          //--- Added L.Ma, 07-10-2014
    (
        IOobject
        (
            "sourcePV",
            U.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    scaledPV_                          //--- Added L.Ma, 07-10-2014
    (
        IOobject
        (
            "scaledPV",
            U.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    )

    /*
    YWI_                          //--- Added L.Ma, 06-03-2015
    (
        IOobject
        (
            "YWI",
            U.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("YWI", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0.0)
    ),
    */
/*
    YuWI_                          //--- Added L.Ma, 06-03-2015
    (
        IOobject
        (
            "YuWI",
            U.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("YuWI", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0.0)
    ),
*/
/*
    YbWI_                          //--- Added L.Ma, 06-03-2015
    (
        IOobject
        (
            "YbWI",
            U.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("YbWI", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0.0)
    ),
*/
/*
    scaledVarPV_                          //--- Added L.Ma, 07-10-2014
    (
        IOobject
        (
            "scaledVarPV",
            U.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE //mb191213
        ),
        mesh_,
        dimensionedScalar("scaledVarPV", dimless, 0.0)
    ),
    */
/*
    Zeta_                          //--- Added L.Ma, 06-03-2015
    (
        IOobject
        (
            "Zeta",
            U.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE //mb191213
        ),
        mesh_,
        dimensionedScalar("Zeta", dimless, 0.0)
    ),
*/
/*
    PVeta_                          //--- Added L.Ma, 06-03-2015
    (
        IOobject
        (
            "PVeta",
            U.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE //mb191213
        ),
        mesh_,
        dimensionedScalar("PVeta", dimless, 0.0)
    ),
*/
/* //200130 mb
    ChiZ_                          //--- Added L.Ma, 06-03-2015
    (
        IOobject
        (
            "chiZ",
            U.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("chiZ", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0.0)
    ),

    ChiPV_                          //--- Added L.Ma, 06-03-2015
    (
        IOobject
        (
            "chiPV",
            U.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("chiZ", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0.0)
    )
    */
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField> Foam::turbulenceModel::phi() const
{
    return phi_;
}


bool Foam::turbulenceModel::read()
{
    return regIOobject::read();
}


void Foam::turbulenceModel::validate()
{}


void Foam::turbulenceModel::correct()
{
    if (mesh_.changing())
    {
        y_.correct();
    }
}

//void Foam::turbulenceModel::correctVarZ()
//{}

//void Foam::turbulenceModel::correctChiZ()
//{}

//void Foam::turbulenceModel::correctChiPV()
//{}

//void Foam::turbulenceModel::correctVarPV() 
//{}
// ************************************************************************* //

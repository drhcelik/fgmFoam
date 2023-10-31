/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/


#include "PVTable.H"
#include <stdio.h>

namespace Foam
{

defineTypeNameAndDebug(PVTable, 0);
defineRunTimeSelectionTable(PVTable, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
//--- PV tables are read here
PVTable::PVTable(const fvMesh& mesh, const word& tablePath, const word& tableName)
:
   IOdictionary
   (
      IOobject
      (
         tableName,
         tablePath,     // --- Changed L.Ma, 22-10-2014, the absolute path of tables
         mesh,
         IOobject::READ_IF_PRESENT,
         IOobject::NO_WRITE
      )
   ),
   tableName_(tableName)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

PVTable::~PVTable()
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::PVTable>
Foam::PVTable::New(const fvMesh& mesh, const word& tablePath, const word& tableName)
{
    const word modelName
    (
        IOdictionary
        (
            IOobject
            (
                "PVtableProperties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        ).lookup("interpolationType")
    );

    dictionaryConstructorTable::iterator cstrIter =
    dictionaryConstructorTablePtr_->find(modelName);

    return cstrIter()(mesh, tablePath, tableName);
}

// * * * * * * * * * * * * * *  Member Functions * * * * * * * * * * * * * * //

} // End Foam namespace

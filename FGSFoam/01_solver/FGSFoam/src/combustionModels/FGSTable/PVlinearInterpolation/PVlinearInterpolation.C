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


#include "PVlinearInterpolation.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const scalarList PVlinearInterpolation::defaultPVList(0.0); //200127 mb

defineTypeNameAndDebug(PVlinearInterpolation, 0);
addToRunTimeSelectionTable(PVTable, PVlinearInterpolation, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

PVlinearInterpolation::PVlinearInterpolation(const fvMesh& mesh, const word& tablePath, const word& tableName)
:
    PVTable(mesh, tablePath, tableName),
    PVtableValues_(this->lookupOrDefault<scalarList>(tableName, defaultPVList)) //201027 mb
{
  if (PVtableValues_ == defaultPVList)
  {
    Info << "*********** table '" << tableName << "' not found! ***********" << endl;
    Info << "'tablePath': " << tablePath << endl;
    Info << "Please specify correct 'tablePath' or check the FGS tables!" << endl;
    Info << "**********************************************************" << endl;
  }

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

PVlinearInterpolation::~PVlinearInterpolation()
{}

// * * * * * * * * * * * * * * * Member Function * * * * * * * * * * * * * * //

//- Linear interpolation for ScalarPV table 
inline scalar PVlinearInterpolation::PVinterpolate(const List<int>& ub, const scalarList& pos) const
{
      scalar c0 = PVtableValues_[ub[0]-1]*(1-pos[0]) + PVtableValues_[ub[0]]*(pos[0]); //200129 mb
      return c0;



}

} // End Foam namespace

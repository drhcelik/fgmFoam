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


#include "linearInterpolation.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const List<scalarList> linearInterpolation::defaultList(0.0); //--- Changed L.Ma, 01-10-2014 to store the 4D integrated FGS talbe. //200127 mb: 2D

defineTypeNameAndDebug(linearInterpolation, 0);
addToRunTimeSelectionTable(FGSTable, linearInterpolation, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

linearInterpolation::linearInterpolation(const fvMesh& mesh, const word& tablePath, const word& tableName)
:
    FGSTable(mesh, tablePath, tableName),
    tableValues_(this->lookupOrDefault<List<scalarList> >(tableName, defaultList)) //--- Changed L.Ma, 01-10-2014 to store the 4D integrated FGS talbe. //200127 mb: 2D
{
  if (tableValues_ == defaultList)
  {
    Info << "*********** table '" << tableName << "' not found!********" << endl;
    Info << "'tablePath': " << tablePath << endl;
    Info << "Please specify correct 'tablePath' or check the FGS tables!" << endl;
    Info << "**********************************************************" << endl;    
  }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

linearInterpolation::~linearInterpolation()
{}

// * * * * * * * * * * * * * * * Member Function * * * * * * * * * * * * * * //

//- 2D interpolation (bilinear integration in Z,C space)  for integrated FGS table
inline scalar linearInterpolation::interpolate(const List<int>& ub, const scalarList& pos) const
{
      scalar c0 = tableValues_[ub[0]-1][ub[1]-1]*(1-pos[0]) + tableValues_[ub[0]][ub[1]-1]*pos[0]; //200128 mb
      scalar c1 = tableValues_[ub[0]-1][ub[1]]*(1-pos[0]) + tableValues_[ub[0]][ub[1]]*pos[0];
      return c0*(1-pos[1]) + c1*pos[1];

}

List<scalarList>  linearInterpolation::tableValues() const //200127 mb
{
  return tableValues_;
}

} // End Foam namespace

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

const List<List<List<scalarList> > > linearInterpolation::defaultList(0.0); //--- Changed L.Ma, 01-10-2014 to store the 4D integrated FGM talbe.

defineTypeNameAndDebug(linearInterpolation, 0);
addToRunTimeSelectionTable(FGMTable, linearInterpolation, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

linearInterpolation::linearInterpolation(const fvMesh& mesh, const word& tablePath, const word& tableName)
:
    FGMTable(mesh, tablePath, tableName),
    tableValues_(this->lookupOrDefault<List<List<List<scalarList> > > >(tableName, defaultList)) //--- Changed L.Ma, 01-10-2014 to store the 4D integrated FGM talbe.
{
  if (tableValues_ == defaultList)
  {
    Info << "*********** table '" << tableName << "' not found!********" << endl;
    Info << "'tablePath': " << tablePath << endl;
    Info << "Please specify correct 'tablePath' or check the FGM tables!" << endl;
    Info << "**********************************************************" << endl;    
  }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

linearInterpolation::~linearInterpolation()
{}

// * * * * * * * * * * * * * * * Member Function * * * * * * * * * * * * * * //

//- 4D interpolation for integrated FGM table
inline scalar linearInterpolation::interpolate(const List<int>& ub, const scalarList& pos) const
{
   // Perform 4d-linear interpolation
   // for varPV
   scalar c000, c100, c010, c110, c001, c101, c011, c111;
   
    if (ub[0] == 0) //- No PV variance
    { 
      c000 = tableValues_[ub[0]][ub[1] -1][ub[2] -1][ub[3] -1];
      c100 = tableValues_[ub[0]][ub[1]][ub[2] -1][ub[3] -1];
      c010 = tableValues_[ub[0]][ub[1] -1][ub[2]][ub[3] -1];
      c110 = tableValues_[ub[0]][ub[1]][ub[2]][ub[3] -1]; 
      c001 = tableValues_[ub[0]][ub[1] -1][ub[2] -1][ub[3]];
      c101 = tableValues_[ub[0]][ub[1]][ub[2] -1][ub[3]];
      c011 = tableValues_[ub[0]][ub[1] -1][ub[2]][ub[3]];
      c111 = tableValues_[ub[0]][ub[1]][ub[2]][ub[3]];
         
    }
    else
    {     
      //- c000 = c0000*(1-pos[0]) + c1000*pos[0]
      c000 = tableValues_[ub[0] -1][ub[1] -1][ub[2] -1][ub[3]-1]*(1-pos[0])+ tableValues_[ub[0]][ub[1] -1][ub[2] -1][ub[3] -1]*pos[0];
      //- c100 = c0100*(1-pos[0]) + c1100*pos[0]
      c100 = tableValues_[ub[0] -1][ub[1]][ub[2] -1][ub[3] -1]*(1-pos[0])+ tableValues_[ub[0]][ub[1]][ub[2] -1][ub[3] -1]*pos[0];
      //- c010 = c0010*(1-pos[0]) + c1010*pos[0]
      c010 = tableValues_[ub[0] -1][ub[1] -1][ub[2]][ub[3] -1]*(1-pos[0])+ tableValues_[ub[0]][ub[1] -1][ub[2]][ub[3] -1]*pos[0];
      //- c110 = c0110*(1-pos[0]) + c1110*pos[0]
      c110 = tableValues_[ub[0] -1][ub[1]][ub[2]][ub[3]-1]*(1-pos[0])+ tableValues_[ub[0]][ub[1]][ub[2]][ub[3] -1]*pos[0]; 
      //- c001 = c0001*(1-pos[0]) + c1001*pos[0]
      c001 = tableValues_[ub[0] -1][ub[1] -1][ub[2] -1][ub[3]]*(1-pos[0])+ tableValues_[ub[0]][ub[1] -1][ub[2] -1][ub[3]]*pos[0];
      //- c101 = c0101*(1-pos[0]) + c1101*pos[0]
      c101 = tableValues_[ub[0] -1][ub[1]][ub[2] -1][ub[3]]*(1-pos[0])+ tableValues_[ub[0]][ub[1]][ub[2] -1][ub[3]]*pos[0];
      //- c011 = c0011*(1-pos[0]) + c1011*pos[0]
      c011 = tableValues_[ub[0] -1][ub[1] -1][ub[2]][ub[3]]*(1-pos[0])+ tableValues_[ub[0]][ub[1] -1][ub[2]][ub[3]]*pos[0];
      //- c111 = c0111*(1-pos[0]) + c1111*pos[0]
      c111 = tableValues_[ub[0] -1][ub[1]][ub[2]][ub[3]]*(1-pos[0])+ tableValues_[ub[0]][ub[1]][ub[2]][ub[3]]*pos[0];
    }
   // for PV
     scalar c00 = c000*(1-pos[1]) + c100*pos[1];
     scalar c10 = c010*(1-pos[1]) + c110*pos[1];
     scalar c01 = c001*(1-pos[1]) + c101*pos[1];
     scalar c11 = c011*(1-pos[1]) + c111*pos[1];
   
   // for varZ 
     scalar c0 = c00*(1-pos[2]) + c10*pos[2];
     scalar c1 = c01*(1-pos[2]) + c11*pos[2];
   
   // for Z
   return c0*(1-pos[3]) + c1*pos[3];
}

List<List<List<scalarList> > >  linearInterpolation::tableValues() const
{
  return tableValues_;
}

} // End Foam namespace

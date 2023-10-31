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


#include "tableSolver.H"

namespace Foam
{
namespace combustionModels
{

defineTypeNameAndDebug(tableSolver, 0);
defineRunTimeSelectionTable(tableSolver, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

tableSolver::tableSolver(const fvMesh& mesh, const wordList& tableNames)
:
   IOdictionary
   (
      IOobject
      (
         "tableProperties",
         mesh.time().constant(),
         mesh,
         IOobject::MUST_READ_IF_MODIFIED,
         IOobject::NO_WRITE
      )
   ),
   tablePath_(this->lookup("tablePath")),
   tableNames_(tableNames),
   paramNames_(2),
   tables_(tableNames_.size() + 1)

{
    forAll(tableNames_, i)
    {
        Info << "Reading table: " << tableNames_[i] << endl;

        tableNames_[i] = tableNames_[i] + "_table";

        tables_.set(i, FGSTable::New(mesh, tablePath_, tableNames_[i]));

    }

    paramNames_[0] = "PV_param";
    paramNames_[1] = "Z_param";

    forAll(paramNames_, i)
    {
    	params_.append(scalarList(this->lookup(paramNames_[i])));
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

tableSolver::~tableSolver()
{}

// * * * * * * * * * * * * * *  Member Functions * * * * * * * * * * * * * * //

List<int> tableSolver::upperBounds(const scalarList& x) const
{
    // Upper, lower bound for table Interpolation and temp Value for bisection method
	List<int> ub(paramNames_.size(), 0);
    int tlb, tub;
	int newVal = 0;

    // Determine upper bounds and interpolation weights for table interpolation
    // Bisection Method
    for (register label j=0; j<paramNames_.size(); j++)
    {
	   tub = params_[j].size() - 1;
   	   tlb = 0;
	   
	   if (tub == tlb)   // In case the last parameter has only one value
	   {
	     ub[j] = tlb;
	     continue;
	   }
	   
   	   while (tub-1 != tlb)
       {
          newVal = (tub + tlb)/2;
          if (x[j] < params_[j][newVal])
     	      tub = newVal;
          else
      	      tlb = newVal;
       }
   	   ub[j] = tub;
    }
   	// Return upper bounds
    return ub;
}

scalarList tableSolver::position(const List<int>& ub, const scalarList& x) const
{
    scalarList pos(paramNames_.size(), 0.0);

    for (register label j=0; j<paramNames_.size(); j++)
    {
        if (ub[j] == 0)
	{
	  pos[j] = 0;
	}
	else
	{
          pos[j] = (x[j] - params_[j][ub[j]-1]) / (params_[j][ub[j]] - params_[j][ub[j]-1]);
	}
    }

    return pos;
}

scalar tableSolver::interpolate(const List<int>& ub , const scalarList& pos, const label& i) const
{
    return tables_[i].interpolate(ub, pos);
}

List<scalarList>  tableSolver::tableValues( const label& i) const //200127 mb
{
  return tables_[i].tableValues();
}

int tableSolver::sizeTableNames() const
{
        return tableNames_.size();
}

} // End combustionModels namespace
} // End Foam namespace

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License

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
    Foam::combustionModels::FGMSolver

Description
    Class that holds the table entries and performs the interpolation.
    based on the flameletFoam: 
    (http://openfoamwiki.net/index.php/Extend-bazaar/solvers/combustion/flameletFoam).

@author Likun Ma, Delft University of Technology
@email  malikun-2005@hotmail.com
@version 14.11.2014

\*---------------------------------------------------------------------------*/

#ifndef FGMTABLE_H_
#define FGMTABLE_H_

#include "fvMesh.H"
#include "Time.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "labelList.H"
#include "scalarList.H"

namespace Foam
{

class fvMesh;

class FGMTable
:
	public IOdictionary
{
private:

	word tableName_;

public:

    //- Runtime type information
    TypeName("FGMTable");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        FGMTable,
        dictionary,
        (
            const fvMesh &mesh,
            const word& tablePath,
            const word& tableName
        ),
        (mesh, tablePath, tableName)
    );

    //- Constructor
    FGMTable(const fvMesh& mesh, const word& tablePath, const word& tableName);

    //- Construct and return a clone
    autoPtr<FGMTable> clone() const;

    //- Destructor
    virtual ~FGMTable();


    // Selectors
    //- Return a reference to the selected tableSolver
    static autoPtr<FGMTable> New(const fvMesh& mesh, const word& tablePath, const word& tableName);

    virtual scalar interpolate(const List<int>& ub, const scalarList& pos) const = 0;
    
    virtual List<List<List<scalarList> > > tableValues() const = 0;
};

} // End Foam namespace

#endif /* FGMTABLE_H_ */

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
    Foam::combustionModels::PVSolver

Description
    Class that holds the table entries and performs the interpolation.
    based on the flameletFoam: 
    (http://openfoamwiki.net/index.php/Extend-bazaar/solvers/combustion/flameletFoam).

@author Likun Ma, Delft University of Technology
@email  malikun-2005@hotmail.com
@version 14.11.2014


\*---------------------------------------------------------------------------*/

#ifndef PVTABLE_H_
#define PVTABLE_H_

#include "fvMesh.H"
#include "Time.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "labelList.H"
#include "scalarList.H"

namespace Foam
{

class fvMesh;

class PVTable
:
	public IOdictionary
{
private:

	word tableName_;

public:

    //- Runtime type information
    TypeName("PVTable");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        PVTable,
        dictionary,
        (
            const fvMesh &mesh,
	    const word& tablePath,
            const word& tableName
        ),
        (mesh, tablePath, tableName)
    );

    //- Constructor
    PVTable(const fvMesh& mesh, const word& tablePath, const word& tableName); 

    //- Construct and return a clone
    autoPtr<PVTable> clone() const;

    //- Destructor
    virtual ~PVTable();


    // Selectors
    //- Return a reference to the selected tableSolver
    static autoPtr<PVTable> New(const fvMesh& mesh, const word& tablePath, const word& tableName);

    virtual scalar PVinterpolate(const List<int>& ub, const scalarList& pos) const = 0;

};

} // End Foam namespace

#endif /* PVTABLE_H_ */

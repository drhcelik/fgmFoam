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
    Foam::combustionModels::linearInter

Description
    Linear Table Interpolation
    based on the flameletFoam: 
    (http://openfoamwiki.net/index.php/Extend-bazaar/solvers/combustion/flameletFoam).

@author Likun Ma, Delft University of Technology
@email  malikun-2005@hotmail.com
@version 14.11.2014

\*---------------------------------------------------------------------------*/

#ifndef LINEARINTERPOLATION_H_
#define LINEARINTERPOLATION_H_

#include "typeInfo.H"
#include "FGSTable.H"

namespace Foam
{

class fvMesh;

class linearInterpolation
:
    public FGSTable
{

private:

    //- Disallow default bitwise copy construct
    linearInterpolation(const linearInterpolation&);

    //- Disallow default bitwise assignment
    void operator=(const linearInterpolation&);

    //- Table Values  
    List<scalarList> tableValues_ ; //--- Changed L.Ma, 01-10-2014 to store the 4D integrated FGS table //200127 mb: 2D

public:

        static const List<scalarList> defaultList; //--- Changed L.Ma, 01-10-2014 to store the 4D integrated FGS table //200127 mb: 2D

        //- Runtime type information
        TypeName("linearInterpolation");

        // Constructors

        //- Construct from components
        linearInterpolation(const fvMesh& mesh, const word& tablePath, const word& tableName);

        //- Destructor
        virtual ~linearInterpolation();

        //- Member functions

        //- Interpolate tableValues_ in a 2D-Table
        virtual inline scalar interpolate(const List<int>& ub, const scalarList& pos) const;
	
	//- Return table values
	virtual List<scalarList> tableValues() const; //201027 mb
};

} // End Foam namespace

#endif /* LINEARINTERPOLATION_H_ */

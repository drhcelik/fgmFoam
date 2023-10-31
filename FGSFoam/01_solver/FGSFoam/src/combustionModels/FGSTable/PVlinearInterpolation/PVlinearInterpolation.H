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

@author Likun Ma, Delft University of Technology
@email  malikun-2005@hotmail.com
@version 14.11.2014

\*---------------------------------------------------------------------------*/

#ifndef LINEARINTERPOLATION_H_
#define LINEARINTERPOLATION_H_

#include "typeInfo.H"
#include "PVTable.H"

namespace Foam
{

class fvMesh;

class PVlinearInterpolation
:
    public PVTable
{

private:

    //- Disallow default bitwise copy construct
    PVlinearInterpolation(const PVlinearInterpolation&);

    //- Disallow default bitwise assignment
    void operator=(const PVlinearInterpolation&);

    //- Tables values for PVmin and PVMax
    scalarList PVtableValues_ ; //200127 mb
    
public:

         static const scalarList defaultPVList; //200127 mb
 
        //- Runtime type information
        TypeName("PVlinearInterpolation");

        // Constructors

        //- Construct from components
        PVlinearInterpolation(const fvMesh& mesh, const word& tablePath, const word& tableName);

       //- Destructor
        virtual ~PVlinearInterpolation();

        //- Member functions

        //- Interpolate tableValues_ in a 1D-Table
        virtual inline scalar PVinterpolate(const List<int>& ub, const scalarList& pos) const; //200127 mb


};

} // End Foam namespace

#endif /* LINEARINTERPOLATION_H_ */

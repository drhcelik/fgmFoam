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
    Foam::combustionModels::FGSModel

Description
    Base class for FGSModel.

@author Likun Ma, Delft University of Technology
@email  malikun-2005@hotmail.com
@version 14.11.2014

SourceFiles
    FGSModel.C

\*---------------------------------------------------------------------------*/

#ifndef FGSModel_H
#define FGSModel_H

#include "ThermoCombustion.H"
#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "zeroGradientFvPatchFields.H"
#include "surfaceFieldsFwd.H"
#include "fvMatricesFwd.H"
#include "tableSolver.H"
#include "PVtableSolver.H"      //--- Added L.Ma, 07-10-2014
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "autoPtr.H"
#include "turbulenceModel.H"   //senbin

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                  Class FGSModel Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionThermo>
class FGSModel
:
    public ThermoCombustion<ReactionThermo>
{

private:

    // Protected data

        //- Tables
        tableSolver solver_;

        //- PV tables
        PVtableSolver PVsolver_;
        
        //- Temperature
        volScalarField& T_;     // Added L.Ma, 02-10-2014   

        //- Pressure
        volScalarField& psi_;     // Added L.Ma, 23-10-2014

        //- Thermal conductivity
        volScalarField& alpha_;     // Added L.Ma, 12-10-2014
        
        //- Viscosity
        volScalarField& mu_;     // Added L.Ma, 12-10-2014
        
        //- Progress variable source term
        volScalarField& sourcePV_;

        //- Scaled Progress variable
        volScalarField& scaledPV_;	
	
        //- Mixture fraction
        const volScalarField& Z_;
	
        //- Progress variable
        volScalarField& PV_;

        //- Disallow copy construct
        FGSModel(const FGSModel&);

        //- Disallow default bitwise assignment
        void operator=(const FGSModel&);


public:

    // Constructors

        //- Runtime type information
        TypeName("FGSModel");
    
        //- Variable names of the integrated FGS talbe
        hashedWordList tables();
        
        //- PVmin and PVmax
        hashedWordList PVtables();

        //- Construct from components
        FGSModel(const word& modelType, 
		 ReactionThermo& thermo,
                 const compressibleTurbulenceModel& turb,
                 const word& combustionProperties);

        //- Destructor
        virtual ~FGSModel();

        //- Correct combustion rate
        virtual void correct();

        //- Correct density during PISO
        Switch correctDensity();

        //- Find unscaled PV from boundary and initial conditions
        virtual void findUscaledPV();

        //- Fuel consumption rate matrix 
        virtual tmp<fvScalarMatrix> R(volScalarField& Y) const;   

	// virtual Qdot, heat release rate [kg/m3/s]
        virtual tmp<volScalarField> Qdot() const;

        // Progress variable source   
        virtual tmp<volScalarField> SourcePV() const;         

        //- Enthalpy source term
        inline tmp<volScalarField> Sh() const;	

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#ifdef NoRepository
#   include "FGSModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

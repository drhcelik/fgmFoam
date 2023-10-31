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

\*---------------------------------------------------------------------------*/

#include "FGSModel.H"
#include "reactingMixture.H"
#include "volFields.H"
#include "hashedWordList.H"

namespace Foam
{
namespace combustionModels
{
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo>
FGSModel<ReactionThermo>::FGSModel
(
    const word& modelType, ReactionThermo& thermo,
    const compressibleTurbulenceModel& turb,
    const word& combustionProperties
)
:
    ThermoCombustion<ReactionThermo>(modelType, thermo, turb),
    solver_(tableSolver(this->mesh(), tables())),
    PVsolver_(PVtableSolver(this->mesh(), PVtables())),        
     
    T_(this->thermo().T()),                 //--- senbin 
    psi_(const_cast<volScalarField&>(this->thermo().psi())),
    alpha_(const_cast<volScalarField&>(this->thermo().alpha())),           //-

    mu_(const_cast<volScalarField&>(this->thermo().mu()())),                  
    sourcePV_(const_cast<volScalarField&>(this->turbulence().sourcePV())),                
    scaledPV_(const_cast<volScalarField&>(this->turbulence().scaledPV())),                  
    Z_(const_cast<volScalarField&>(this->turbulence().Z())),
    PV_(const_cast<volScalarField&>(this->turbulence().PV()))                       
    
{
findUscaledPV();
}

// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

template<class ReactionThermo>
FGSModel<ReactionThermo>::~FGSModel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ReactionThermo>
hashedWordList FGSModel<ReactionThermo>::tables()
{
	// - FGS tables to be read
	hashedWordList tableNames;
	tableNames.clear();
	tableNames.append("T");   //--- Read 'T' table
	tableNames.append("psi");
	tableNames.append("mu");
	tableNames.append("alpha");
	tableNames.append("SourcePV"); //--- Read source term of PV

        return tableNames;
}

template<class ReactionThermo>
hashedWordList FGSModel<ReactionThermo>::PVtables()
{
        // - PV tables to be read
	hashedWordList PVtableNames;
	PVtableNames.clear();
	PVtableNames.append("PVmin");
	PVtableNames.append("PVmax");

        return PVtableNames;
}


template<class ReactionThermo>
void FGSModel<ReactionThermo>::correct()
{
    Info << "Entering FGSModel correct()" << endl;    //--- Added L.Ma, 14-10-2014
    const scalarField& ZCells = Z_.internalField();  
    const scalarField& PVCells = PV_.internalField();
    
    scalarField& TCells = T_.ref();
    scalarField& sourcePVCells = sourcePV_.ref();       
    scalarField& scaledPVCells = scaledPV_.ref();
  
    scalarField& muCells = mu_.ref();                   
    scalarField& alphaCells = alpha_.ref();             
    scalarField& psiCells = psi_.ref();             

    //- Update the species and temperature field
    if(this->active())
    {
       //- 1D table for minimum and maximum progress variable, (varPV, PV)
       scalarList y(1, 0.0); //200128 mb
           
       //- 2D independent variables (PV, Z)
       scalarList x(2, 0.0); //200128 mb
       
       //- Upper bounds and position for table interpolation
       List<int>  ubIF_, ubPVIF_, ubP_, ubPVP_;     //senbin 06-11-18         //--- Added L.Ma, 06-10-2014
       scalarList posIF_, posPVIF_, posP_, posPVP_;  //senbin 06-11-18          //--- Added L.Ma, 06-10-2014

       //- Minimum and Maximum unscaled progress variable
       scalar PVMinCells, PVMaxCells; //200128 mb, works
       scalar pPVMin, pPVMax;
       scalar sPV;

    
       //- A small number to prevent divide by zero
       scalar smallVar(1e-5);

       forAll(ZCells, cellI)    // 1st: for internal field
       {
            //- Calculate scaled progress variable
            y[0] = max(min(ZCells[cellI], 1.0), 0.0);
      
            ubPVIF_ = PVsolver_.upperBounds(y);
            posPVIF_ = PVsolver_.position(ubPVIF_, y);
                 
            //- Update minimum and maximum unscaled progress variable
            PVMinCells = PVsolver_.interpolate(ubPVIF_, posPVIF_, 0);
            PVMaxCells = PVsolver_.interpolate(ubPVIF_, posPVIF_, 1);
                
            //- Update scaled progress variable
            sPV = (PVCells[cellI] - PVMinCells)/max((PVMaxCells - PVMinCells),smallVar);
            scaledPVCells[cellI] =  max(min(sPV, 1.0), 0.0);    // = C

                     
             //- Prepare independent variables for table look up          
	
             x[0] = scaledPVCells[cellI]; // 200128 mb
             x[1] = max(min(ZCells[cellI], 1.0), 0.0);

	     //- find up-bound and pos for FGS table interpolation
             ubIF_ = solver_.upperBounds(x);
             posIF_ = solver_.position(ubIF_, x);
 
             //- Update temperature, sourcePV, flame sensor internal field
             TCells[cellI] = solver_.interpolate(ubIF_, posIF_, 0);               
	     psiCells[cellI] = solver_.interpolate(ubIF_, posIF_, 1);     
             muCells[cellI] = solver_.interpolate(ubIF_, posIF_, 2);
	     alphaCells[cellI] = solver_.interpolate(ubIF_, posIF_, 3);
	     sourcePVCells[cellI] = solver_.interpolate(ubIF_, posIF_, 4);

       }

       // Interpolate for patches
       forAll(T_.boundaryField(), patchi)    // Changed L.Ma, 24-06-2014, 2nd, for patches.
       {
	  const fvPatchScalarField& pPV = PV_.boundaryField()[patchi];
          const fvPatchScalarField& pZ = Z_.boundaryField()[patchi];
	      
          fvPatchScalarField& pT = T_.boundaryFieldRef()[patchi];        // Added L.Ma, 24-06-2014 
          fvPatchScalarField& psourcePV = sourcePV_.boundaryFieldRef()[patchi];        // Added L.Ma, 24-06-2014
          fvPatchScalarField& pscaledPV = scaledPV_.boundaryFieldRef()[patchi];        // Added L.Ma, 24-06-2014	  
	  
          fvPatchScalarField& pMu = mu_.boundaryFieldRef()[patchi];        // Added L.Ma, 24-06-2014
          fvPatchScalarField& pAlpha = alpha_.boundaryFieldRef()[patchi];        // Added L.Ma, 24-06-2014
          fvPatchScalarField& ppsi = psi_.boundaryFieldRef()[patchi];
          
          forAll(pZ, facei)
          {
      
               //- Calculate scalaed progress variable 
               y[0] = max(min(pZ[facei], 1.0), 0.0); //200128 mb
      
               ubPVP_ = PVsolver_.upperBounds(y);
               posPVP_ = PVsolver_.position(ubPVP_, y);
      
               //- Update minimum and maximum unscaled progress variable
               pPVMin = PVsolver_.interpolate(ubPVP_, posPVP_, 0);
               pPVMax = PVsolver_.interpolate(ubPVP_, posPVP_, 1);
      
               //- Update scaled progress variable
	       sPV =  (pPV[facei] - pPVMin)/max((pPVMax - pPVMin),smallVar);
               pscaledPV[facei] =  max(min(sPV, 1.0), 0.0);
	    
		     
	       //- Prepare independent variables for table look up
	       x[0] = pscaledPV[facei]; //200128 mb
	       x[1] = max(min(pZ[facei], 1.0), 0.0);

	       ubP_ = solver_.upperBounds(x);
	       posP_ = solver_.position(ubP_, x);

	       //- Updated temperature patch field
	       pT[facei] = solver_.interpolate(ubP_, posP_, 0);
	       ppsi[facei] = solver_.interpolate(ubP_, posP_, 1);  
	       pMu[facei] = solver_.interpolate(ubP_, posP_, 2);
	       pAlpha[facei] = solver_.interpolate(ubP_, posP_, 3);
	       psourcePV[facei] = solver_.interpolate(ubP_, posP_, 4);
             }
       }
  
    }
}

template<class ReactionThermo>
Switch FGSModel<ReactionThermo>::correctDensity()
{
	return true;
}
//---------------------------------------------------------------
template<class ReactionThermo>
void FGSModel<ReactionThermo>::findUscaledPV()
{
    Info << endl;
    Info << "Initializing unscaled progress variable" << endl;
    Info << endl;

    List<int> uPVIF_, uPVP_;
    scalarList poPVIF_, poPVP_;
    scalar PVMinIF, PVMaxIF, PVMinp, PVMaxp;
    scalarList Ind_PV(1, 0.0);
    
    const scalarField& ZIF = Z_.internalField(); 
    const scalarField& scaledPVIF = scaledPV_.internalField(); 
    scalarField& PVIF = PV_.ref();
    

    forAll(ZIF, cellI)
    {
	Ind_PV[0] = max(min(ZIF[cellI], 1.0), 0.0);
	
        uPVIF_ = PVsolver_.upperBounds(Ind_PV);
        poPVIF_ = PVsolver_.position(uPVIF_, Ind_PV);

        //- Find minimum and maximum unscaled progress variable
        PVMinIF = PVsolver_.interpolate(uPVIF_, poPVIF_, 0);
        PVMaxIF = PVsolver_.interpolate(uPVIF_, poPVIF_, 1);
	
        PVIF[cellI]  =  PVMinIF + scaledPVIF[cellI]*(PVMaxIF-PVMinIF);
    } 

    forAll(Z_.boundaryField(), patchi)    
    {
      const fvPatchScalarField& Zp = Z_.boundaryField()[patchi];      
      const fvPatchScalarField& scaledPVp = scaledPV_.boundaryField()[patchi];  
      fvPatchScalarField& PVp = PV_.boundaryFieldRef()[patchi];        
	
      forAll(Zp , facei)
      {
	  Ind_PV[0] = max(min(Zp[facei], 1.0), 0.0); //200129 mb

	  uPVP_ = PVsolver_.upperBounds(Ind_PV);
	  poPVP_ = PVsolver_.position(uPVP_, Ind_PV);

	  //- Find minimum and maximum unscaled progress variable
	  PVMinp = PVsolver_.interpolate(uPVP_, poPVP_, 0);
	  PVMaxp = PVsolver_.interpolate(uPVP_, poPVP_, 1);

	  PVp[facei] = PVMinp + scaledPVp[facei]*(PVMaxp - PVMinp);
      }
    }
}
//\-------------------------------------------------------------------

template<class ReactionThermo>
Foam::tmp<Foam::fvScalarMatrix>
FGSModel<ReactionThermo>::R
(
    volScalarField& Y                 //---Changed L.Ma, 03-07-2014
) const
{
    tmp<fvScalarMatrix> tSu(new fvScalarMatrix(Y, dimMass/dimTime));
    return tSu;
}

// Progress variable source 
template<class ReactionThermo>
Foam::tmp<Foam::volScalarField>
FGSModel<ReactionThermo>::SourcePV() const          //Added L.Ma, 08-10-2014
{
  return sourcePV_;
}

template<class ReactionThermo>
Foam::tmp<Foam::volScalarField>
FGSModel<ReactionThermo>::Sh() const
{
    tmp<volScalarField> tSh
    (
        new volScalarField
        (
            IOobject
            (
                "Sh",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    return tSh;
}

template<class ReactionThermo>
Foam::tmp<Foam::volScalarField>
FGSModel<ReactionThermo>::Qdot() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "dQ",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("dQ", dimEnergy/dimVolume/dimTime, 0.0), //here the unit is right for HRR
            zeroGradientFvPatchScalarField::typeName
        )
    );

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

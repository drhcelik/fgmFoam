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

#include "FGMModel.H"
#include "reactingMixture.H"
#include "volFields.H"
#include "hashedWordList.H"

namespace Foam
{
namespace combustionModels
{
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo>
FGMModel<ReactionThermo>::FGMModel
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
    YWI_(const_cast<volScalarField&>(this->turbulence().YWI())),
    YuWI_(const_cast<volScalarField&>(this->turbulence().YuWI())),
    YbWI_(const_cast<volScalarField&>(this->turbulence().YbWI())),
    scaledVarPV_(const_cast<volScalarField&>(this->turbulence().scaledVarPV())),                
    Zeta_(const_cast<volScalarField&>(this->turbulence().Zeta())),
    PVeta_(const_cast<volScalarField&>(this->turbulence().PVeta())),   
    Z_(const_cast<volScalarField&>(this->turbulence().Z())),
    varZ_(const_cast<volScalarField&>(this->turbulence().varZ())),
    PV_(const_cast<volScalarField&>(this->turbulence().PV())),                
    varPV_(const_cast<volScalarField&>(this->turbulence().varPV())),          
    
    useProgressVariableVariance_(this->coeffs().lookupOrDefault("useProgressVariableVariance",true))      
{
findUscaledPV();
}

// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

template<class ReactionThermo>
FGMModel<ReactionThermo>::~FGMModel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ReactionThermo>
hashedWordList FGMModel<ReactionThermo>::tables()
{
        // - Read useProgressVariableVariance_  
        read();
	// - FGM tables to be read
	hashedWordList tableNames;
	tableNames.clear();
	tableNames.append("T");   //--- Read 'T' table
	tableNames.append("psi");
	tableNames.append("mu");
	tableNames.append("alpha");
	tableNames.append("SourcePV"); //--- Read source term of PV
        
        if (useProgressVariableVariance_)
	{
           tableNames.append("YWI"); //--- Read source term of PV
           tableNames.append("YuWI"); //--- To calculate scalar dissipation rate for Yc
           tableNames.append("YbWI"); //--- To calculate scalar dissipation rate for Yc   
	}
        return tableNames;
}

template<class ReactionThermo>
hashedWordList FGMModel<ReactionThermo>::PVtables()
{
        // - Read useProgressVariableVariance_  
        read();
        // - PV tables to be read
	hashedWordList PVtableNames;
	PVtableNames.clear();
	PVtableNames.append("PVmin");
	PVtableNames.append("PVmax");
    
	if (useProgressVariableVariance_)
	{
	   PVtableNames.append("Yu2I");
	   PVtableNames.append("YuYbI");
	   PVtableNames.append("Yb2I");
	}
       return PVtableNames;
}


template<class ReactionThermo>
void FGMModel<ReactionThermo>::correct()
{
    Info << "Entering FGMModel correct()" << endl;    //--- Added L.Ma, 14-10-2014
    const scalarField& ZCells = Z_.internalField();  
    const scalarField& varZCells = varZ_.internalField(); 
    const scalarField& PVCells = PV_.internalField();
    const scalarField& varPVCells = varPV_.internalField();
    
    scalarField& TCells = T_.ref();
    scalarField& sourcePVCells = sourcePV_.ref();       
    scalarField& scaledPVCells = scaledPV_.ref();       
    scalarField& YWICells = YWI_.ref();               
    scalarField& YuWICells = YuWI_.ref();              
    scalarField& YbWICells = YbWI_.ref();             
    scalarField& ZetaCells = Zeta_.ref();           
    scalarField& PVetaCells = PVeta_.ref();     
    scalarField& scaledVarPVCells = scaledVarPV_.ref();  
    scalarField& muCells = mu_.ref();                   
    scalarField& alphaCells = alpha_.ref();             
    scalarField& psiCells = psi_.ref();             

    //- Update the species and temperature field
    if(this->active())
    {
       //- 2D table for minimum and maximum progress variable, (varPV, PV)
       scalarList y(2, 0.0); 
           
       //- 4D independent variables (varPV, PV, varZ, Z)
       scalarList x(4, 0.0);
       
       //- Upper bounds and position for table interpolation
       List<int>  ubIF_, ubPVIF_, ubP_, ubPVP_;     //senbin 06-11-18         //--- Added L.Ma, 06-10-2014
       scalarList posIF_, posPVIF_, posP_, posPVP_;  //senbin 06-11-18          //--- Added L.Ma, 06-10-2014

       //- Minimum and Maximum unscaled progress variable
       scalar PVMinCells, PVMaxCells, pPVMin, pPVMax;
       scalar sPV, sVarPV;
       //- Scale progress variable and its variance
       scalar fc, gc, hc;
       scalar Yu2I, YuYbI, Yb2I;
       
       //- A small number to prevent divide by zero
       scalar smallVar(1e-5);

       forAll(ZCells, cellI)    // 1st: for internal field
       {

            //- Scaled varZ (Zeta) is stored in the pre-integrated FGM table as a independent variable
            ZetaCells[cellI] = max(min(varZCells[cellI]/max(ZCells[cellI]*(1.0 - ZCells[cellI]), smallVar),0.99), 0.0);
                   
            //- Calculate scaled progress variable   
	    y[0] = ZetaCells[cellI];
            y[1] = max(min(ZCells[cellI], 1.0), 0.0);
      
            ubPVIF_ = PVsolver_.upperBounds(y);
            posPVIF_ = PVsolver_.position(ubPVIF_, y);
                 
            //- Update minimum and maximum unscaled progress variable
            PVMinCells = PVsolver_.interpolate(ubPVIF_, posPVIF_, 0);
            PVMaxCells = PVsolver_.interpolate(ubPVIF_, posPVIF_, 1);
                
            //- Update scaled progress variable
            sPV = (PVCells[cellI] - PVMinCells)/max((PVMaxCells - PVMinCells),smallVar);
            scaledPVCells[cellI] =  max(min(sPV, 1.0), 0.0);    // = C''2

             //- Update scaled progress variable variance?
	     if (useProgressVariableVariance_)
	     {
	        //- calculate scaled progress variable variance from unscaled propress variance
	        //- Reference 1: A progress variable approach based on premixed flamelets for turbulent combustion modeling, B.A. Albrecht, W.J.S.Ramaekers et al
	        //- Eq.(10) and Eq.(17)  (in numerator should be -2gc)
	        //- Reference 2: A premixed Flamelet-PDF Model for Biomass Combustion in a Grate Furnace, Energy & Fuels, 2008, Albrecht, Oijen et al
	        //- Eq.(10) and Eq.(11) (in denomenator should be (Yb-Yu)^2)
	        //- Reference 3: Private communication, Scaling of the reaction progress variable
		     
	        Yu2I = PVsolver_.interpolate(ubPVIF_, posPVIF_, 2);
	        YuYbI = PVsolver_.interpolate(ubPVIF_, posPVIF_, 3);
	        Yb2I = PVsolver_.interpolate(ubPVIF_, posPVIF_, 4);
		     
	        fc = Yu2I; gc = YuYbI - Yu2I; hc = Yb2I - 2.0*YuYbI + Yu2I;
                sVarPV =  (varPVCells[cellI] + sqr(PVCells[cellI]) - fc - 2.0*gc*scaledPVCells[cellI])/max(hc, smallVar) - sqr(scaledPVCells[cellI]);
	        scaledVarPVCells[cellI] = max(sVarPV, 0.0);
                if (hc < smallVar)        
	            {
	              scaledVarPVCells[cellI] = 0;
		    }     
             }
	     else
	     {
	         scaledVarPVCells[cellI] = 0.0;
	     }
                   
             PVetaCells[cellI] = max(min(scaledVarPVCells[cellI]/max(scaledPVCells[cellI]*(1.0 -scaledPVCells[cellI]), smallVar), 0.99), 0.0);
                                    
             //- Prepare independent variables for table look up          
             x[0] = PVetaCells[cellI];
             x[1] = scaledPVCells[cellI];
	     x[2] = max(min(ZetaCells[cellI], 0.99), 0.0);
             x[3] = max(min(ZCells[cellI], 1.0), 0.0);
         
	     //- find up-bound and pos for FGM table interpolation
             ubIF_ = solver_.upperBounds(x);
             posIF_ = solver_.position(ubIF_, x);
 
             //- Update temperature, sourcePV, flame sensor internal field
             TCells[cellI] = solver_.interpolate(ubIF_, posIF_, 0);               
	     psiCells[cellI] = solver_.interpolate(ubIF_, posIF_, 1);     
             muCells[cellI] = solver_.interpolate(ubIF_, posIF_, 2);
	     alphaCells[cellI] = solver_.interpolate(ubIF_, posIF_, 3);
	     sourcePVCells[cellI] = solver_.interpolate(ubIF_, posIF_, 4);
	     
	     if (useProgressVariableVariance_)
             {
		     YWICells[cellI] = solver_.interpolate(ubIF_, posIF_, 5);
		     YuWICells[cellI] = solver_.interpolate(ubIF_, posIF_, 6);
		     YbWICells[cellI] = solver_.interpolate(ubIF_, posIF_, 7);
	     }
	     else
             {
		     YWICells[cellI] = 0.0;
		     YuWICells[cellI] = 0.0;
		     YbWICells[cellI] = 0.0;
	     }  
       }

       // Interpolate for patches
       forAll(T_.boundaryField(), patchi)    // Changed L.Ma, 24-06-2014, 2nd, for patches.
       {
          const fvPatchScalarField& pvarPV = varPV_.boundaryField()[patchi];
	  const fvPatchScalarField& pPV = PV_.boundaryField()[patchi];
          const fvPatchScalarField& pvarZ = varZ_.boundaryField()[patchi];
          const fvPatchScalarField& pZ = Z_.boundaryField()[patchi];
	      
          fvPatchScalarField& pT = T_.boundaryFieldRef()[patchi];        // Added L.Ma, 24-06-2014 
          fvPatchScalarField& psourcePV = sourcePV_.boundaryFieldRef()[patchi];        // Added L.Ma, 24-06-2014
          fvPatchScalarField& pscaledPV = scaledPV_.boundaryFieldRef()[patchi];        // Added L.Ma, 24-06-2014
          fvPatchScalarField& pYWI = YWI_.boundaryFieldRef()[patchi];
	  fvPatchScalarField& pYuWI = YuWI_.boundaryFieldRef()[patchi];
    	  fvPatchScalarField& pYbWI = YbWI_.boundaryFieldRef()[patchi];
	  fvPatchScalarField& pZeta = Zeta_.boundaryFieldRef()[patchi];	  
	  fvPatchScalarField& pPVeta = PVeta_.boundaryFieldRef()[patchi];	  
	  
          fvPatchScalarField& pscaledVarPV = scaledVarPV_.boundaryFieldRef()[patchi];        // Added L.Ma, 24-06-2014
          fvPatchScalarField& pMu = mu_.boundaryFieldRef()[patchi];        // Added L.Ma, 24-06-2014
          fvPatchScalarField& pAlpha = alpha_.boundaryFieldRef()[patchi];        // Added L.Ma, 24-06-2014
          fvPatchScalarField& ppsi = psi_.boundaryFieldRef()[patchi];
          
          forAll(pZ, facei)
          {
      
               //- Scaled varZ (Zeta) is stored in the pre-integrated FGM table as a independent variable
               pZeta[facei] = max(min(pvarZ[facei]/max(pZ[facei]*(1.0 - pZ[facei]), smallVar), 0.99), 0.0);

               //- Calculate scaled progress variable 
               y[0] = pZeta[facei];
               y[1] = max(min(pZ[facei], 1.0), 0.0);
      
               ubPVP_ = PVsolver_.upperBounds(y);
               posPVP_ = PVsolver_.position(ubPVP_, y);
      
               //- Update minimum and maximum unscaled progress variable
               pPVMin = PVsolver_.interpolate(ubPVP_, posPVP_, 0);
               pPVMax = PVsolver_.interpolate(ubPVP_, posPVP_, 1);
      
               //- Update scaled progress variable
	       sPV =  (pPV[facei] - pPVMin)/max((pPVMax - pPVMin),smallVar);
               pscaledPV[facei] =  max(min(sPV, 1.0), 0.0);

               //- How about the progress variable variance?
	       if (useProgressVariableVariance_)
	       {    
		     Yu2I = PVsolver_.interpolate(ubPVP_, posPVP_, 2);
		     YuYbI = PVsolver_.interpolate(ubPVP_, posPVP_, 3);
		     Yb2I = PVsolver_.interpolate(ubPVP_, posPVP_, 4);
		       
		     fc = Yu2I; gc = YuYbI - Yu2I; hc = Yb2I - 2.0*YuYbI + Yu2I;
                     sVarPV =  (pvarPV[facei] + sqr(pPV[facei]) - fc - 2.0*gc*pscaledPV[facei])/max(hc, smallVar) - sqr(pscaledPV[facei]);
			
		     pscaledVarPV[facei] =  max(sVarPV, 0.0);
		              
		     //- For the situation where Pvmin = PVmax, this normally happens at the fuel or oxidizer inlet
		     //- The scaled progress variable is set to zero in this situation
		     if (hc < smallVar)        
		        {
	                   pscaledVarPV[facei] = 0;
		        }		       
	       }
	       else
	       {
		     pscaledVarPV[facei] = 0.0;
	       }
		     
	       pPVeta[facei] = max(min(pscaledVarPV[facei]/max(pscaledPV[facei]*(1.0 - pscaledPV[facei]), smallVar), 0.99), 0.0);
		     
	       //- Prepare independent variables for table look up
	       x[0] = pPVeta[facei];
	       x[1] = pscaledPV[facei];
	       x[2] = max(min(pZeta[facei], 0.99), 0.0);
	       x[3] = max(min(pZ[facei], 1.0), 0.0);
	       ubP_ = solver_.upperBounds(x);
	       posP_ = solver_.position(ubP_, x);

	       //- Updated temperature patch field
	       pT[facei] = solver_.interpolate(ubP_, posP_, 0);
	       ppsi[facei] = solver_.interpolate(ubP_, posP_, 1);  
	       pMu[facei] = solver_.interpolate(ubP_, posP_, 2);
	       pAlpha[facei] = solver_.interpolate(ubP_, posP_, 3);
	       psourcePV[facei] = solver_.interpolate(ubP_, posP_, 4);
	       
	       if (useProgressVariableVariance_)
	       {     
	             pYWI[facei] = solver_.interpolate(ubP_, posP_, 5);
	             pYuWI[facei] = solver_.interpolate(ubP_, posP_, 6);
	             pYbWI[facei] = solver_.interpolate(ubP_, posP_, 7);
	       }
	       else
	       {
	             pYWI[facei] = 0.0;
	             pYuWI[facei] = 0.0;
	             pYbWI[facei] = 0.0;
	       }
             }
       }
  
    }
}

template<class ReactionThermo>
Switch FGMModel<ReactionThermo>::correctDensity()
{
	return true;
}
//---------------------------------------------------------------
template<class ReactionThermo>
void FGMModel<ReactionThermo>::findUscaledPV()
{
    Info << endl;
    Info << "Initializing unscaled progress variable" << endl;
    Info << endl;

    List<int> uPVIF_, uPVP_;
    scalarList poPVIF_, poPVP_;
    scalar PVMinIF, PVMaxIF, PVMinp, PVMaxp, sZ;
    scalarList Ind_PV(2, 0.0);
    
    const scalarField& ZIF = Z_.internalField();
    const scalarField& varZIF = varZ_.internalField();  
    const scalarField& scaledPVIF = scaledPV_.internalField(); 
    scalarField& PVIF = PV_.ref();
    
    forAll(ZIF, cellI)
    {
        sZ = varZIF[cellI]/max(ZIF[cellI]*(1.0 - ZIF[cellI]), SMALL);
        Ind_PV[0] = max(min(sZ, 0.99), 0.0);
	Ind_PV[1] = max(min(ZIF[cellI], 1.0), 0.0);
	
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
      const fvPatchScalarField& varZp = varZ_.boundaryField()[patchi];      
      const fvPatchScalarField& scaledPVp = scaledPV_.boundaryField()[patchi];  
      fvPatchScalarField& PVp = PV_.boundaryFieldRef()[patchi];        
	
      forAll(Zp , facei)
      {
	  sZ = varZp[facei]/max(Zp[facei]*(1.0 - Zp[facei]), SMALL);
	  //- Mixture fraction
          Ind_PV[0] = max(min(sZ, 0.99), 0.0);
	  Ind_PV[1] = max(min(Zp[facei], 1.0), 0.0);

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
FGMModel<ReactionThermo>::R
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
FGMModel<ReactionThermo>::SourcePV() const          //Added L.Ma, 08-10-2014
{
  return sourcePV_;
}

template<class ReactionThermo>
Foam::tmp<Foam::volScalarField>
FGMModel<ReactionThermo>::Sh() const
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
FGMModel<ReactionThermo>::Qdot() const
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

template<class ReactionThermo>
bool FGMModel<ReactionThermo>::read()
{
    if (ThermoCombustion<ReactionThermo>::read())
    {
        this->coeffs().lookup("useProgressVariableVariance") >> useProgressVariableVariance_; 
        return true;
    }
    else
    {
        return false;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

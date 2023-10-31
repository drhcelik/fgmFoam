/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

Description
    Calculates species mass fractions and thermodynamic properties
    from given Z, varZ , PV and varPV fields

    @author Likun Ma, Delft University of Technology
    @email  malikun-2005@hotmail.com
    @version 14.11.2014
    @version 30.01.2020 mb

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "tableSolver.H"
#include "PVtableSolver.H"

#include "fvCFD.H"
#include "CombustionModel.H"  //changed by senbin
#include "psiReactionThermo.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "addRegionOption.H"
    argList::addOption
    (
        "fields",
        "list",
        "specify a list of fields to be reconstructed. Eg, '(U T p)' - "
        "regular expressions not currently supported"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

    HashSet<word> selectedFields;
    if (args.optionFound("fields"))
    {
        args.optionLookup("fields")() >> selectedFields;
    }

    //reading thermo properties, added by senbin
    autoPtr<psiReactionThermo> pThermo(psiReactionThermo::New(mesh));
    psiReactionThermo& thermo = pThermo();

    const IOdictionary combProps
    (
        IOobject
        (
            "combustionProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const IOdictionary tableProps
    (
        IOobject
        (
            "tableProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const IOdictionary PVtableProps
    (
        IOobject
        (
            "PVtableProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    
    const word combTypeName = combProps.lookup("combustionModel");
    const label tempOpen = combTypeName.find('<');
    const word modelType = combTypeName(0, tempOpen);
    dictionary coeffs_(combProps.subDict(modelType + "Coeffs"));

    Switch useProgressVariableVariance_(coeffs_.lookup("useProgressVariableVariance"));

    //- 2D table for minimum and maximum progress variable
    scalarList y(2, 0.0);
    
    //- 4D independent variables (Z, varZ, PV, varPV)
    scalarList x(4, 0.0);
    
    double Zeta, ZetaC;


    List<List<int> > ubIF(mesh.cells().size()), ubPVIF(mesh.cells().size());
    List<scalarList> posIF(mesh.cells().size()), posPVIF(mesh.cells().size());
    List<List<int> > ubP(mesh.faces().size()), ubPVP(mesh.faces().size());
    List<scalarList> posP(mesh.faces().size()), posPVP(mesh.faces().size());

    //- Minimum and Maximum unscaled progress variable
    scalar PVMinCells, PVMaxCells, pPVMin, pPVMax;

    //- Scaled progress variable   
    scalar sPVCells, psPV;
    scalar sVarPVCells(0.0), psVarPV(0.0);  

    //- Scale progress variable and its variance
    scalar fc, gc, hc;
    scalar Yu2I, YuYbI, Yb2I;

    wordList tableNames(thermo.composition().species());
    tableNames.append("T");   //--- Read 'T' table
    tableNames.append("SourcePV");
    Foam::combustionModels::tableSolver solver(Foam::combustionModels::tableSolver(mesh, tableNames));
    
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

    Foam::combustionModels::PVtableSolver PVsolver(Foam::combustionModels::PVtableSolver(mesh, PVtableNames));
    
    
    PtrList<volScalarField>& Y(thermo.composition().Y());
//    volScalarField& hs(thermo.he());   // Changed l.Ma, 25-06-2014    

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info << nl << "Time = " << runTime.timeName() << nl << endl;

        volScalarField Z
        (
            IOobject
            (
                "Z",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );
        volScalarField varZ
        (
            IOobject
            (
                "varZ",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );
        volScalarField PV
        (
            IOobject
            (
                "PV",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

        volScalarField varPV
        (
            IOobject
            (
                "varPV",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

       // Interpolate for internal Field
       forAll(Y, i)
       {
    	  scalarField& YCells = Y[i].ref();   //internalField(); senbin

          forAll(Z, cellI)
          {
        	 if (i == 0)  // Zeta and scaledPV determined once is enough
        	 {
		   
		   //- Scaled varZ (Zeta) is stored in the pre-integrated FGM table as an independent variable
                   Zeta = varZ[cellI]/max(Z[cellI]*(1 - Z[cellI]), SMALL);

		   //- Calculate scaled progress variable   
		   y[0] = max(min(Zeta, 0.99), 0.0);
                   y[1] = max(min(Z[cellI], 1.0), 0.0);
		   
                   ubPVIF[cellI] = PVsolver.upperBounds(y);
                   posPVIF[cellI] = PVsolver.position(ubPVIF[cellI], y);
      
                   //- Update minimum and maximum unscaled progress variable
                   PVMinCells = PVsolver.interpolate(ubPVIF[cellI], posPVIF[cellI], 0);
                   PVMaxCells = PVsolver.interpolate(ubPVIF[cellI], posPVIF[cellI], 1);
      
                   //- Update scalaed progress variable
                   sPVCells = (PV[cellI] - PVMinCells)/max((PVMaxCells - PVMinCells),SMALL);

                   //- How about the progress variable variance?
		   if (useProgressVariableVariance_)
		   {
		     //- calculate scaled progress variable variance from unscaled propress variance
		     //- Reference 1: A progress variable approach based on premixed flamelets for turbulent combustion modeling, B.A. Albrecht, W.J.S.Ramaekers et al
		     //- Eq.(10) and Eq.(17)
		     //- Reference 2: A premixed Flamelet-PDF Model for Biomass Combustion in a Grate Furnace, Energy & Fuels, 2008, Albrecht, Oijen et al
		     //- Eq.(10) and Eq.(11)
		     
		     Yu2I = PVsolver.interpolate(ubPVIF[cellI], posPVIF[cellI], 2);
		     YuYbI = PVsolver.interpolate(ubPVIF[cellI], posPVIF[cellI], 3);
		     Yb2I = PVsolver.interpolate(ubPVIF[cellI], posPVIF[cellI], 4);
		     
		     fc = Yu2I; gc = YuYbI - Yu2I; hc = Yb2I - 2*YuYbI + Yu2I;
                     sVarPVCells =  (varPV[cellI] + sqr(PV[cellI]) - fc - 2*gc*sPVCells)/max(hc, SMALL) - sqr(sPVCells); 
		   }
		   else
		   {
		     sVarPVCells = 0.0;
		   }

                   ZetaC = sVarPVCells/max(sPVCells*(1 - sPVCells), SMALL);

                   x[0] = max(min(ZetaC, 0.99), 0.0);
                   x[1] = max(min(sPVCells, 1.0), 0.0);
		   x[2] = max(min(Zeta, 0.99), 0.0);
                   x[3] = max(min(Z[cellI], 1.0), 0.0);
		   //- find up-bound and pos for table interpolation
                   ubIF[cellI] = solver.upperBounds(x);
                   posIF[cellI] = solver.position(ubIF[cellI], x);

		   YCells[cellI] = solver.interpolate(ubIF[cellI], posIF[cellI], i);
        	 }
                 //- Update species
        	 YCells[cellI] = solver.interpolate(ubIF[cellI], posIF[cellI], i); 
          }
       }

       // Interpolate for patches
       forAll(Z.boundaryField(), patchi)    // Changed L.Ma, 24-06-2014
       {
          const fvPatchScalarField& pvarPV = varPV.boundaryField()[patchi];
	  const fvPatchScalarField& pPV = PV.boundaryField()[patchi];
          const fvPatchScalarField& pvarZ = varZ.boundaryField()[patchi];
          const fvPatchScalarField& pZ = Z.boundaryField()[patchi];

          forAll(Y, i)
          {
              fvPatchScalarField& pY = Y[i].boundaryFieldRef()[patchi];

              forAll(pY , facei)
              {
             	 if (i == 0)
             	 {
                      Zeta = pvarZ[facei]/max(pZ[facei]*(1 - pZ[facei]), SMALL);
                     //- Calculate scaled progress variable 
                     y[0] = max(min(Zeta, 0.99), 0.0);
                     y[1] = max(min(pZ[facei], 1.0), 0.0);
      
                     ubPVP[facei] = PVsolver.upperBounds(y);
                     posPVP[facei] = PVsolver.position(ubPVP[facei], y);
      
                     //- Update minimum and maximum unscaled progress variable
                     pPVMin = PVsolver.interpolate(ubPVP[facei], posPVP[facei], 0);
                     pPVMax = PVsolver.interpolate(ubPVP[facei], posPVP[facei], 1);
      
                     //- Update scalaed progress variable
                     psPV = (pPV[facei] - pPVMin)/max((pPVMax - pPVMin),SMALL);

                     //- How about the progress variable variance?
		     if (useProgressVariableVariance_)
		     {    
		       Yu2I = PVsolver.interpolate(ubPVP[facei], posPVP[facei], 2);
		       YuYbI = PVsolver.interpolate(ubPVP[facei], posPVP[facei], 3);
		       Yb2I = PVsolver.interpolate(ubPVP[facei], posPVP[facei], 4);
		       
		       fc = Yu2I; gc = YuYbI - Yu2I; hc = Yb2I - 2*YuYbI + Yu2I;
                       psVarPV =  (pvarPV[facei] + sqr(pPV[facei]) - fc - 2*gc*psPV)/max(hc, SMALL) - sqr(psPV); 
		     }
		     else
		     {
		       psVarPV = 0.0;
		     }

                     ZetaC = psVarPV/max(psPV*(1 - psPV), SMALL);

                     x[0] = max(min(ZetaC, 1.0), 0.0);
		     x[1] = max(min(psPV, 1.0), 0.0);
                     x[2] = max(min(Zeta, 0.99), 0.0);
                     x[3] = max(min(pZ[facei], 1.0), 0.0);
		     
                     ubP[facei] = solver.upperBounds(x);
                     posP[facei] = solver.position(ubP[facei], x);

		     pY[facei] = solver.interpolate(ubP[facei], posP[facei], i);
             	 }
                 //- update speces
                 pY[facei] = solver.interpolate(ubP[facei], posP[facei], i);
             }
          }
       }
  

        if (selectedFields.empty())
        {
        	forAll(Y, i)
            {
               Info << "Writing field " << thermo.composition().Y()[i].name() << endl;
        	   thermo.composition().Y()[i].write();
            }
        }
        else
        {
        	forAll(Y, i)
            {
        	   if (selectedFields[thermo.composition().Y()[i].name()])
        	   {
                   Info << "Writing field " << thermo.composition().Y()[i].name() << endl;
        		   thermo.composition().Y()[i].write();
        	   }
            }
        }

    }

	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	Info<< "End\n" << endl;

	return 0;
}

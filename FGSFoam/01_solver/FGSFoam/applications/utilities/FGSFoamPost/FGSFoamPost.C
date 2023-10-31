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
    from given Z and PV fields

    @author Likun Ma, Delft University of Technology
    @email  malikun-2005@hotmail.com
    @version 14.11.2014
    @version 30.1.2020 mb

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

    //- 1D table for minimum and maximum progress variable
    scalarList y(1, 0.0); //200128 mb
    
    //- 2D independent variables (Z, varZ, PV, varPV)
    scalarList x(2, 0.0); //200128 mb
    
    List<List<int> > ubIF(mesh.cells().size()), ubPVIF(mesh.cells().size());
    List<scalarList> posIF(mesh.cells().size()), posPVIF(mesh.cells().size());
    List<List<int> > ubP(mesh.faces().size()), ubPVP(mesh.faces().size());
    List<scalarList> posP(mesh.faces().size()), posPVP(mesh.faces().size());

    //- Minimum and Maximum unscaled progress variable
    scalar PVMinCells, PVMaxCells, pPVMin, pPVMax;

    //- Scaled progress variable   
    scalar sPVCells, psPV;  

    //- Scale progress variable and its variance
    wordList tableNames(thermo.composition().species());
    tableNames.append("T");   //--- Read 'T' table
    tableNames.append("SourcePV");
    Foam::combustionModels::tableSolver solver(Foam::combustionModels::tableSolver(mesh, tableNames));
    
    hashedWordList PVtableNames;
    PVtableNames.clear();
    PVtableNames.append("PVmin");
    PVtableNames.append("PVmax");

    Foam::combustionModels::PVtableSolver PVsolver(Foam::combustionModels::PVtableSolver(mesh, PVtableNames));
    
    
    PtrList<volScalarField>& Y(thermo.composition().Y());   

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

       // Interpolate for internal Field
       forAll(Y, i)
       {
    	  scalarField& YCells = Y[i].ref();   //internalField(); senbin

          forAll(Z, cellI)
          {
        	 if (i == 0)  // Zeta and scaledPV determined once is enough
        	 {
		   
		   //- Calculate scaled progress variable   
                   y[0] = max(min(Z[cellI], 1.0), 0.0); //200128 mb
		   
                   ubPVIF[cellI] = PVsolver.upperBounds(y);
                   posPVIF[cellI] = PVsolver.position(ubPVIF[cellI], y);
      
                   //- Update minimum and maximum unscaled progress variable
                   PVMinCells = PVsolver.interpolate(ubPVIF[cellI], posPVIF[cellI], 0);
                   PVMaxCells = PVsolver.interpolate(ubPVIF[cellI], posPVIF[cellI], 1);
      
                   //- Update scalaed progress variable
                   sPVCells = (PV[cellI] - PVMinCells)/max((PVMaxCells - PVMinCells),SMALL);

                   x[0] = max(min(sPVCells, 1.0), 0.0); //200128 mb
                   x[1] = max(min(Z[cellI], 1.0), 0.0);

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
	  const fvPatchScalarField& pPV = PV.boundaryField()[patchi];
          const fvPatchScalarField& pZ = Z.boundaryField()[patchi];

          forAll(Y, i)
          {
              fvPatchScalarField& pY = Y[i].boundaryFieldRef()[patchi];

              forAll(pY , facei)
              {
             	 if (i == 0)
             	 {
                     //- Calculate scaled progress variable 
                     y[0] = max(min(pZ[facei], 1.0), 0.0); //200128 mb
                     
      
                     ubPVP[facei] = PVsolver.upperBounds(y);
                     posPVP[facei] = PVsolver.position(ubPVP[facei], y);
      
                     //- Update minimum and maximum unscaled progress variable
                     pPVMin = PVsolver.interpolate(ubPVP[facei], posPVP[facei], 0);
                     pPVMax = PVsolver.interpolate(ubPVP[facei], posPVP[facei], 1);
      
                     //- Update scaled progress variable
                     psPV = (pPV[facei] - pPVMin)/max((pPVMax - pPVMin),SMALL);

		     x[0] = max(min(psPV, 1.0), 0.0);
                     x[1] = max(min(pZ[facei], 1.0), 0.0);
		     
                     ubP[facei] = solver.upperBounds(x);
                     posP[facei] = solver.position(ubP[facei], x);

		     pY[facei] = solver.interpolate(ubP[facei], posP[facei], i);
             	 }
                 //- update species
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

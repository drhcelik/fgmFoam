/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

Application
    FGMFoam

Description
    Solver for combustion with the steady laminar FGM model.

    @author Likun Ma, Delft University of Technology
    @email  malikun-2005@hotmail.com
    @version 14.11.2014
    @version 30.01.2020 mb

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "turbulentFluidThermoModel.H"  //senbin
#include "psiReactionThermo.H"
#include "CombustionModel.H"
#include "pressureControl.H"
#include "multivariateScheme.H"
#include "pimpleControl.H"
#include "localEulerDdtScheme.H"
#include "fvOptions.H"
#include "fvcSmooth.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readGravitationalAcceleration.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"

//    pimpleControl pimple(mesh);

    turbulence->validate();

    #include "compressibleCourantNo.H"
    #include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "compressibleCourantNo.H"
        #include "setDeltaT.H"

        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "rhoEqn.H"

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
    	{
            #include "UEqn.H"
            #include "ZEqn.H" 
            #include "PVEqn.H"

	    // --- Pressure corrector loop
            while (pimple.correct())
    	    {
                if (pimple.consistent())
                {
                    #include "pcEqn.H"
                }
                else
                {
                    #include "pEqn.H"
		}
            }
            
	    if (pimple.turbCorr())
            {
                turbulence->correct();
            }
	}

        rho = thermo.rho();  
	runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

    tmp<fvScalarMatrix> tPVEqn
    (
	(
	  fvm::ddt(rho, PV)
	+ fvm::div(phi, PV)
        - fvm::laplacian( (turbulence->muEff()), PV)   //-- changed by senbin, 2019-06-24 muEff=muSgs()/Sct_ + mu()/Sc_, containing rho already
        ==  
         reaction->SourcePV()   //--- added by senbin, 2019-04-26
	 + fvOptions(rho, PV)
        )
    );
    
    fvScalarMatrix& PVEqn = tPVEqn.ref();
    PVEqn.relax();
    fvOptions.constrain(PVEqn);
    
    PVEqn.solve("PV");
    fvOptions.correct(PV);
    PV = max(   min( PV, 1.0 ), 0.0   );
    Info<< "PV min/max = " << min(PV).value() << "/" << max(PV).value() << endl;

    reaction->correct();   //correct SourcePV, T, sensor
    
    Info<< "T gas min/max = " << min(thermo.T()).value() << "/" << max(thermo.T()).value() << endl;

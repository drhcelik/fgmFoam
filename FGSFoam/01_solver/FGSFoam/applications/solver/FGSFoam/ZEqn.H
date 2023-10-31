
Switch IsPremixed(chemistryProperties.lookup("premixed"));
if (!IsPremixed)
{
   tmp<fvScalarMatrix> tZEqn
   (
        (
	  fvm::ddt(rho, Z)
	+ fvm::div(phi, Z)
        - fvm::laplacian( (turbulence->muEff()), Z) //changed by senbin
        )
    );

    fvScalarMatrix& ZEqn = tZEqn.ref();
    ZEqn.relax();

    ZEqn.solve("Z");

    reaction->correct();
}    
    Z = max(   min( Z, 1.0 ), 0.0   );
    Info<< "Z min/max = " << min(Z).value() << "/" << max(Z).value() << endl;
    

    


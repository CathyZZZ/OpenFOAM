{
    Info << nl << "Solving air flow" << endl;

    pAir.storePrevIter();
//    rho1.storePrevIter();
    // Read controls
    const dictionary& piso = mesh.solutionDict().subDict("PISO");
    int nCorr(readInt(piso.lookup("nCorrectors")));
    int nNonOrthCorr(readInt(piso.lookup("nNonOrthogonalCorrectors")));

    fvVectorMatrix UEqn
    (
     
       fvm::div(phiAir, Uair)
      - fvm::laplacian(muAir, Uair,  "laplacian(muAir,Uair)")
    );
//    pZones.addResistance(UEqn);

//    Info << "muAir = " << muAir << endl;


    UEqn.relax();

    solve(UEqn == -fvc::grad(pAir));

    // --- PISO loop

    for (int corr = 0; corr < nCorr; corr++)
    {
//        pAir.boundaryField().updateCoeffs();
        rUAair = 1.0/UEqn.A();
        Uair = UEqn.H()/UEqn.A();
        Uair.correctBoundaryConditions();

        phiAir = fvc::interpolate(rho*Uair) & mesh.Sf();
 //       rho= pAir/fvc:: interpolate(R*T);

       for (int nonOrth = 0; nonOrth <= nNonOrthCorr; nonOrth++)
        {
            fvScalarMatrix pEqn
            (
                fvm::laplacian(rho/UEqn.A(), pAir) == fvc::div(phiAir)// - mO2source
            );

            pEqn.setReference(pAirRefCell, pAirRefValue);
            pEqn.solve();

            if (nonOrth == nNonOrthCorr)
            {
                phiAir -= pEqn.flux();
            }

         }

//        continuityErrs(phiAir);

//        Info << (Foam::min(pAir).value()) << ", " << (Foam::max(pAir).value()) << endl;


        // Explicitly relax pressure for momentum corrector
        pAir.relax();

        Uair -= fvc::grad(pAir)/UEqn.A();
        Uair.correctBoundaryConditions();
    }
}
    Info << "min,max(Uair): "
         << (Foam::min(Uair)).value() <<" , "
         << (Foam::max(Uair)).value() << endl;

    Info << "min,max(pAir): "
         << (Foam::min(pAir)).value() <<" , "
         << (Foam::max(pAir)).value() << endl;

{

    Info << nl << "Solving energy" << endl;

    fvScalarMatrix TEqn
    (
//        fvm::div(linearInterpolate(Cp)*phiAir, T,"div(phi,T)")
//    fvm::SuSp(-fvc::div(linearInterpolate(Cp)*phiAir), T)
//        fvm::div(linearInterpolate(Cp)*phiAir), T)
      - fvm::laplacian(k, T, "laplacian(k,T)")
    );

    TEqn.solve();

/*
fvScalarMatrix rhoEqn
      (
         pAir/(R * T)
    );

    solve(rhoEqn == rho);
    rho.write();
*/

Info << "T min mean max     = "
         << Foam::gMin(T.internalField()) << "    " 
         << Foam::gAverage(T.internalField()) << "    "
         << Foam::gMax(T.internalField())
         << nl <<endl;

    Info<< " sum TEqn.source() =   " << gSum(TEqn.source()) << nl;
    Info<< endl;
}

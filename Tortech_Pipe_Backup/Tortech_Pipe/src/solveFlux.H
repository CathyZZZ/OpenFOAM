{

    Info << nl << "Solving Flux" << endl;

    fvScalarMatrix NEqn
    (
//       fvm::div(linearInterpolate(Cp)*phiAir, T,"div(phi,T)")
//        fvm::div((rho*phiAir), N, "div((rho*phiAir),N)")
//        fvm::ddt(N)
         fvm::div(phiAir, N, "div(phiAir,N)")
//      +  fvm::div(Uair,N)
      - fvm::laplacian(rho*D, N, "laplacian(rho*D,N)")

    );

    NEqn.solve();


    Info << "N min mean max     = "
         << Foam::gMin(N.internalField()) << "    " 
         << Foam::gAverage(N.internalField()) << "    "
         << Foam::gMax(N.internalField())
         << nl <<endl;

    Info<< " sum NEqn.source() =   " << gSum(NEqn.source()) << nl;
    Info<< endl;
}

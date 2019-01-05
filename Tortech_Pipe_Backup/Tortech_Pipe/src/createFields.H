Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

//dimensionedScalar mu (transportProperties.lookup("mu"));
dimensionedScalar R (transportProperties.lookup("R"));
dimensionedScalar Na (transportProperties.lookup("Na"));


Info<< "Creating Fields for Flow Calculation\n" << endl;
//Info<< muAir << endl;
volScalarField pAir
(
 IOobject
 (
  "pAir",
  runTime.timeName(),
  mesh,
  IOobject::MUST_READ,
  IOobject::AUTO_WRITE
  ),
 mesh
 );

volScalarField muAir
(
 IOobject
 (
  "muAir",
  runTime.timeName(),
  mesh,
  IOobject::READ_IF_PRESENT,
  IOobject::AUTO_WRITE
//  IOobject::NO_WRITE
  ),
 mesh
// mu,
// zeroGradientFvPatchScalarField::typeName
 );
//muAir.correctBoundaryConditions();






Info<< "Creating field Velocity" <<endl;
volVectorField Uair
(
 IOobject
 (
  "Uair",
  runTime.timeName(),
  mesh,
  IOobject::MUST_READ,
  IOobject::AUTO_WRITE
  ),
 mesh
 );

volScalarField rUAair
(
 IOobject
 (
  "(1|A(U))Air",
  runTime.timeName(),
  mesh
  ),
 mesh,
 dimensionedScalar("rUAA", dimTime/dimDensity, 1),
 zeroGradientFvPatchScalarField::typeName
 );

label pAirRefCell = 0;
scalar pAirRefValue = 0.0;
setRefCell (pAir,mesh.solutionDict().subDict("PISO"), pAirRefCell, pAirRefValue);
 
Info<< "Reading cell field T\n" << endl;
volScalarField T
(
 IOobject
 (
  "T",
  runTime.timeName(),
  mesh,
  IOobject::MUST_READ,
  IOobject::AUTO_WRITE
  ),
 mesh
 );

volScalarField rho
(
 IOobject
 (
  "rho",
  runTime.timeName(),
  mesh,
  IOobject::READ_IF_PRESENT,
  IOobject::AUTO_WRITE
  ),
mesh
//zeroGradientFvPatchScalarField::typeName
// pAir/linearInterpolate(R *T) 
 );


volScalarField cp
(
 IOobject
 (
  "cp",
  runTime.timeName(),
  mesh,
  IOobject::READ_IF_PRESENT,
  IOobject::AUTO_WRITE
  ),
mesh
 );





//volScalarField rho("rho", pAir / (R * T));

surfaceScalarField phiAir
(
 IOobject
 (
  "phiAir",
  runTime.timeName(),
  mesh,
  IOobject::READ_IF_PRESENT,
  IOobject::AUTO_WRITE
  ),
 linearInterpolate(rho*Uair) & mesh.Sf()
 );

volScalarField k
(
 IOobject
 (
  "k",
  runTime.timeName(),
  mesh,
  IOobject::MUST_READ,
  IOobject::AUTO_WRITE
  ),
 mesh
 );

Info<< "Reading cell field N\n" << endl;
volScalarField N
(
 IOobject
 (
  "N",
  runTime.timeName(),
  mesh,
  IOobject::MUST_READ,
  IOobject::AUTO_WRITE
  ),
 mesh
 );

volScalarField D
(
 IOobject
 (
  "D",
  runTime.timeName(),
  mesh,
  IOobject::MUST_READ,
  IOobject::AUTO_WRITE
  ),
 mesh
 );

/*
Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
Info<< "Reading field T\n" <<endl;
volScalarField T
(
    IOobject
    (
         "T",
         runTime.timeName(),
         mesh,
         IOobject::MUST_READ,
         IOobject::AUTO_WRITE
     ),
     mesh
);

#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());
*/


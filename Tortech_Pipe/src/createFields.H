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
dimensionedScalar dAir (transportProperties.lookup("dAir"));
//dimensionedScalar dp (transportProperties.lookup("dp"));
dimensionedScalar R (transportProperties.lookup("R"));
dimensionedScalar Na (transportProperties.lookup("Na"));
dimensionedScalar BoltzConst  (transportProperties.lookup("BoltzConst"));

dimensionedScalar ACoeff(transportProperties.lookup("ACoeff"));
dimensionedScalar BCoeff(transportProperties.lookup("BCoeff"));
dimensionedScalar CCoeff(transportProperties.lookup("CCoeff"));
dimensionedScalar DCoeff(transportProperties.lookup("DCoeff"));
dimensionedScalar MWDisp(transportProperties.lookup("MWDisp"));
dimensionedScalar rhoDisp(transportProperties.lookup("rhoDisp"));


dimensionedScalar nodes(transportProperties.lookup("nodes"));

//dimensionedScalar d1(transportProperties.lookup("dpMax"));
dimensionedScalar dimCorr(transportProperties.lookup("dimCorr"));

Switch nucleation(transportProperties.lookup("nucleation"));//switch between true and false.
Switch growth(transportProperties.lookup("growth"));
Switch coagulation(transportProperties.lookup("coagulation"));
Switch monomerBal(transportProperties.lookup("monomerBal"));

double kb = BoltzConst.value();
scalar node = nodes.value(); // number of size nodes
scalar pi = constant ::mathematical::pi;
int node1 = nodes.value();
int node2 = nodes.value();
int node3 = nodes.value();
//calcultation of first node mass
dimensionedScalar d1;
dimensionedScalar v1 = MWDisp/(rhoDisp*Na );
dimensionedScalar m1 = MWDisp/Na;
d1.value() = Foam::pow(v1.value()*6/pi,0.333333333333333);
Info<< nl << "v1 " << v1 << nl <<endl;
dimensionedScalar s1 = pi*Foam::pow(d1,2);
Info<< nl << "m1 " << m1 << nl <<endl;

// setup one dimensional arrays for particle source term calculations
PtrList<dimensionedScalar> dp(node);
PtrList<dimensionedScalar> vp(node);

for (label i = 0; i < node; i++) {
	 dp.set (i,new dimensionedScalar ("dp",dimensionSet(0,1,0,0,0,0,0),0));
         vp.set (i,new dimensionedScalar ("vp",dimensionSet(0,3,0,0,0,0,0),0));
}


// read and create an array for particle size nodes where size has dp[i]
std::ifstream read("size.txt");
if (read.is_open()) {
	for ( label i=0; i<node; i++) {
		read >> dp[i].value();
		vp[i] = pi/6* pow(dp[i],3);
		Info << "diameter for node " << i+1 << "" << dp[i].dimensions() << endl;
	}
}
else {
	Info << "cannot find the file" << endl;
}
read.close();
Info << nl;

//calculate coagulation splitting operators(3dimensioned)
PtrList<PtrList<PtrList<dimensionedScalar> > >Xijk(node1);

forAll(Xijk, k) {
	Xijk.set( k, new PtrList<PtrList<dimensionedScalar> >(node2));
	forAll(Xijk[k],j) {
		Xijk[k].set( j, new PtrList<dimensionedScalar>(node3));
		forAll(Xijk[k][j] ,i) {
	  			Xijk[k][j].set(i, new dimensionedScalar ("X", dimensionSet(0,0,0,0,0,0,0),0));
		}
	}
}

Info << "nodes1" << node1 << endl;
Info << "nodes2" << node2 << endl;
Info << "nodes3" << node3 << endl;
for (label k = 1; k < node; k++) {
    for (label j = 0; j < node; j++) {
        for (label i = 0; i < node; i++) {
            if (k != node-1) {
                if ( (vp[k].value() < (vp[i].value() + vp[j].value())) &&
                     (vp[k+1].value() >(vp[i].value() + vp[j].value() ))) {
                    Xijk[k][j][i].value() =
                        (vp[k+1].value()-vp[i].value()-vp[j].value()) / (vp[k+1].value()-vp[k].value());


                   Info << "vpi="<<vp[i].value()<<nl<<"vpj="<<vp[j]<<nl<<"vpk="<<vp[k]<<endl;
                   Info << "i="<<i<<nl<<"j="<<j<<nl<<"k="<<k<<nl<<"Xijk="<<Xijk[k][j][i].value()<<endl;
                }
                else if ( ( vp[k-1].value() < (vp[i].value() + vp[j].value()) ) &&
                          ( vp[k].value() > (vp[i].value() + vp[j].value()) )) {
                   Xijk[k][j][i].value() = (vp[i].value() + vp[j].value() - vp[k-1].value()) / (vp[k].value() - vp[k-1].value());
               
                
                   Info << "vpi="<<vp[i].value()<<nl<<"vpj="<<vp[j]<<nl<<"vpk="<<vp[k]<<endl;
                   Info << "i="<<i<<nl<<"j="<<j<<nl<<"k="<<k<<nl<<"Xijk="<<Xijk[k][j][i].value()<<endl;
                }
                else {
                    Xijk[k][j][i].value() = 0;
                }
                }
            else if (k == node-1) {
                if (vp[k].value() < (vp[i].value() + vp[j].value())) {
                    Xijk[k][j][i].value() = 1;
                }
                else if ((vp[k].value() > (vp[i].value() + vp[j].value())) &&
                        (vp[k-1].value() < (vp[i].value() + vp[j].value()))) {
                    Xijk[k][j][i].value() = (vp[i].value() + vp[j].value() - vp[k-1].value()) / (vp[k].value() - vp[k-1].value());      
                }
                else {
                    Xijk[k][j][i].value() =0;
                }
            }
            
        }
    }
}


Info<< "Reading field nPart\n" << endl;
PtrList<volScalarField> Nair(nodes.value());
for (label i=0; i<nodes.value(); i++) {
    char rnPart[16];
    label indx = i+1;
    sprintf (rnPart, "nPartNode%d", indx);
    Nair.set( i, 
              new volScalarField
              (
                IOobject
                (
                    rnPart,
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
              )
            );
}


PtrList<volScalarField> NSatk(nodes.value());
for (label i = 0; i < node; i++) {
        char rNSatk[16];
        label indx = i+1;
        printf (rNSatk, "NSatk%d", indx);
        NSatk.set ( i, new volScalarField //("NSatk",dimensionSet(0,0,0,0,0,0,0),0));
                (
                    IOobject
                    (
                        rNSatk,
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh,

                    dimensionedScalar(rNSatk,dimensionSet(0,-3,0,0,0,0,0), 0.0),
                    zeroGradientFvPatchScalarField::typeName
               )
                    );
        //gammaNs.add(Ngamma[j][i]);

}



PtrList<volScalarField> Nzeta(nodes.value());
for (label i = 0; i < node; i++) {
        char rNzeta[16];
        label indx = i+1;
        printf (rNzeta, "Nzeta%d", indx);
        Nzeta.set ( i, new volScalarField //("NSatk",dimensionSet(0,0,0,0,0,0,0),0));
                (
                    IOobject
                    (
                        rNzeta,
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh,

                    dimensionedScalar(rNzeta,dimensionSet(0,0,0,0,0,0,0), 0.0),
                    zeroGradientFvPatchScalarField::typeName
                )
                    );
}


PtrList<PtrList<volScalarField> > Beta(nodes.value());

forAll(Beta, j) {
    Beta.set ( j, new PtrList<volScalarField>(nodes.value()));
    forAll(Beta[j],i) {

    char rBeta[16];
    label indxj = j+1;
    label indxi = i+1;
    sprintf (rBeta, "Beta%d_%d", indxj,indxi);
        Beta[j].set ( i, new volScalarField
                (
                    IOobject
                    (
                        rBeta,
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh,
                    dimensionedScalar(rBeta, dimensionSet(0,3,-1,0,0,0,0), 0.0),
                    zeroGradientFvPatchScalarField::typeName
                )
                    );
        //gammaNs.add(Ngamma[j][i]);
    }
}

//*******************************************************************************
volScalarField S//saturation ratio
(
    IOobject
    (
     "S",
     runTime.timeName(),
     mesh,
     IOobject::NO_READ,
     IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("S", dimensionSet(0,0,0,0,0,0,0), 1.001),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField pSat//saturation vapor pressure
(
     IOobject
    (
     "pSat",
     runTime.timeName(),
     mesh,
     IOobject::NO_READ,
     IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("pSat", dimensionSet(1,-1,-2,0,0,0,0), 0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField theta//nondimensional surface tension
(
      IOobject
     (
      "theta",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
     ),
     mesh,
     dimensionedScalar("theta", dimensionSet(0,0,0,0,0,0,0), 0),
     zeroGradientFvPatchScalarField::typeName
);

volScalarField Jk
(
    IOobject
    (
        "Jk",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Jk", dimensionSet(0,-3,-1,0,0,0,0), 0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField Jk2
(
    IOobject
    (
        "Jk2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Jk2", dimensionSet(0,0,-1,0,0,0,0), 0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField sigma//surface tension field
(
    IOobject
    (
        "sigma",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,

    dimensionedScalar("sigma", dimensionSet(1,0,-2,0,0,0,0), 0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField dpstar//critical cluster size
(
    IOobject
    (
        "dpstar",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("dpstar", dimensionSet(0,1,0,0,0,0,0), 0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField vpstar//critical cluster size (volume)
(
    IOobject
    (
        "vpstar",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("vpstar", dimensionSet(0,3,0,0,0,0,0), 0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField kstar//critical cluster size node
(
    IOobject
    (
        "kstar",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("kstar", dimensionSet(0,0,0,0,0,0,0), 0),
    zeroGradientFvPatchScalarField::typeName
);
//*******************************************************************************
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
/*
    volScalarField Nair
    (
      IOobject
      (
       "Nair",
       runTime.timeName(),
       mesh,
       IOobject::NO_READ,
       IOobject::AUTO_WRITE
      ),
    mesh,
    dimensionedScalar("Nair", dimLength, 0.0),
zeroGradientFvPatchScalarField::typeName 
);
*/


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



volScalarField lambda
(
 IOobject
 (
  "lambda",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 //lambda,
 dimensionedScalar("lambda", dimLength, 0.0),
 zeroGradientFvPatchScalarField::typeName
 );
 
 
 volScalarField kn
 (
 IOobject
 (
 "kn",
  runTime.timeName(),
   mesh,
   IOobject::NO_READ,
   IOobject::AUTO_WRITE
   ),
   mesh,
   dimensionedScalar("kn", dimLength, 0),
   zeroGradientFvPatchScalarField::typeName
   );
 

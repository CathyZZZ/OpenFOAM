    Info<< "Reading	 transportProperties\n" << endl;

    volTensorField directions
    (
        IOobject
        (
            "directions",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedTensor("directns", dimensionSet(0, 0, 0, 0, 0, 0, 0), tensor::zero)
    );

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

    dimensionedScalar rhoAir (transportProperties.lookup("rho"));
    dimensionedScalar mu (transportProperties.lookup("mu"));
    dimensionedScalar CpAir (transportProperties.lookup("Cp"));
    dimensionedScalar ACoeff(transportProperties.lookup("ACoeff"));
    dimensionedScalar BCoeff(transportProperties.lookup("BCoeff"));
    dimensionedScalar CCoeff(transportProperties.lookup("CCoeff"));
    dimensionedScalar DCoeff(transportProperties.lookup("DCoeff"));
    dimensionedScalar MWDisp(transportProperties.lookup("MWDisp"));
    dimensionedScalar rhoDisp(transportProperties.lookup("rhoDisp"));
    dimensionedScalar nodes(transportProperties.lookup("nodes"));
    dimensionedScalar d1(transportProperties.lookup("dpMax"));
    dimensionedScalar dimCorr(transportProperties.lookup("dimCorr"));
    Switch nucleation(transportProperties.lookup("nucleation"));//switch between true and false.
    Switch growth(transportProperties.lookup("growth"));
    Switch coagulation(transportProperties.lookup("coagulation"));
    Switch monomerBal(transportProperties.lookup("monomerBal"));

	scalar node = nodes.value(); // number of size nodes
	int node1 = nodes.value();
	int node2 = nodes.value();
	int node3 = nodes.value();
    const dimensionedScalar epsilon0 = physicalConstant::epsilon0;
    const dimensionedScalar kb = physicalConstant::kb;
    const dimensionedScalar NA = physicalConstant::NA;
	const dimensionedScalar R = physicalConstant::R;
    scalar pi = constant::mathematical::pi;

//calcultation of first node mass
	dimensionedScalar v1 = MWDisp/(rhoDisp*NA);
	dimensionedScalar m1 = MWDisp/NA;
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
					

// read and create a an array for particle size nodes where size has dp[i]
    std::ifstream read("size.txt");
    if (read.is_open()) {
        for (int i=0; i<node; i++) {
            read >> dp[i].value();
			vp[i] = pi/6*pow(dp[i],3);
            Info << "diameter for node " << i+1 << " " << dp[i].dimensions() << " = " << dp[i].value() << " " << vp[i].value() << " " <<endl;
        }
    }
    else {
        Info << "cannot find the file" <<endl;
    }
    read.close();
    Info << nl;

//calculate coagulation splitting operators
	PtrList<PtrList<PtrList<dimensionedScalar> > >Xijk(node1);

	forAll(Xijk, k) {
        Xijk.set (
            k, 
            new PtrList<PtrList<dimensionedScalar> >(node2)
       		 );
			forAll(Xijk[k], j) {
		  	  Xijk[k].set (
		       			 j, 
            				new PtrList<dimensionedScalar>(node3)
		    );
			forAll(Xijk[k][j], i) {
				Xijk[k][j].set (
				    i, 
				    new dimensionedScalar("X", dimensionSet(0,0,0,0,0,0,0),0)

				);			
			}
		}
	}

						Info << "nodes1 " << node1 << endl;
						Info << "nodes2 " << node2 << endl;
						Info << "nodes3 " << node3 << endl;

	for (label k = 1; k < node; k++) {
		for (label i = 0; i < node; i++) {
			for (label j = 0; j < node; j++) {
				if (k == node-1)
				{
					if((vp[k].value() <= (vp[i].value()+vp[j].value())))
					{
						Xijk[i][j][k].value()=1;
					}
					else if ((vp[k-1].value()<= (vp[i].value()+vp[j].value())) && ((vp[i].value()+vp[j].value())<vp[k].value()))
					{
						Xijk[i][j][k].value()=(vp[i].value()+vp[j].value()-vp[k-1].value())/(vp[k].value()-vp[k-1].value());
					}
					else
					{
					Xijk[i][j][k].value()=0;
					}
				}
				else
				{
					if((vp[k].value() <= (vp[i].value()+vp[j].value())) && ((vp[i].value()+vp[j].value())<vp[k+1].value()))
					{
						Xijk[i][j][k].value()=(vp[k+1].value()-vp[i].value()-vp[j].value())/(vp[k+1].value()-vp[k].value());
					}
					else if ((vp[k-1].value()<= (vp[i].value()+vp[j].value())) && ((vp[i].value()+vp[j].value())<vp[k].value()))
					{
						Xijk[i][j][k].value()=(vp[i].value()+vp[j].value()-vp[k-1].value())/(vp[k].value()-vp[k-1].value());
					}
					else
					{
					Xijk[i][j][k].value()=0;
					}
				}
				if (j == 20 && i == 20)
		       		Info << i << " " << j << " " << k << " " << Xijk[i][j][k].value() <<endl;
			}
		}
	}




    Info<< "Reading field nPart\n" << endl;
    multivariateSurfaceInterpolationScheme<scalar>::fieldTable airNs;
    PtrList<volScalarField> Nair(nodes.value());

    for (label i = 0; i < nodes.value(); i++)
    {
	char rnPart[16];
	label indx = i+1;
	sprintf (rnPart, "nPartNode%d", indx);
        Nair.set
        (
            i,
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
        airNs.add(Nair[i]);//airNs calculate total number of particles, Nair[i] calculate number of particles at node i.
    }

    multivariateSurfaceInterpolationScheme<scalar>::fieldTable NSatks;
    PtrList<volScalarField> NSatk(nodes.value());

    for (label i = 0; i < nodes.value(); i++)
    {
	char rNSatk[16];
	label indx = i+1;
    sprintf (rNSatk, "NSatk%d", indx);
        NSatk.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    rNSatk,
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            )
        );
        NSatks.add(NSatk[i]);
    }


    multivariateSurfaceInterpolationScheme<scalar>::fieldTable gammaNs;
    PtrList<PtrList<volScalarField> > Ngamma(nodes.value());

    forAll(Ngamma, k) {
	    Ngamma.set (
	    	k, 
	    	new PtrList<volScalarField>(nodes.value())
	    );
		forAll(Ngamma[k], i) {
			char rgamma[16];
            label indxk = k+1;
            label indxi = i+1;
			sprintf (rgamma, "gamma%d_%d", indxk, indxi);
		    Ngamma[k].set
		    (
		        i,
		        new volScalarField
		        (
		            IOobject
		            (
		                rgamma,
		                runTime.timeName(),
		                mesh,
		                IOobject::MUST_READ,
		                IOobject::NO_WRITE
		            ),
		            mesh
		        )
		    );
		    gammaNs.add(Ngamma[k][i]);
		}
	}



    multivariateSurfaceInterpolationScheme<scalar>::fieldTable zetaNs;
    PtrList<volScalarField> Nzeta(nodes.value());

    for (label i = 0; i < nodes.value(); i++)
    {
	char rzeta[16];
	label indx = i+1;
	sprintf (rzeta, "zeta%d", indx);
        Nzeta.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    rzeta,
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            )
        );
        zetaNs.add(Nzeta[i]);
    }




    Info<< "Creating Fields for Flow Calculation\n" << endl;
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
            "mu",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        mu,
        zeroGradientFvPatchScalarField::typeName
    );
    muAir.correctBoundaryConditions();

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        rhoAir,
        zeroGradientFvPatchScalarField::typeName
    );
    muAir.correctBoundaryConditions();


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
        linearInterpolate(rhoAir*Uair) & mesh.Sf()
    );

    // inverse matrix diagonal field rUA
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

    // Pressure reference point in air
    label pAirRefCell = 0;
    scalar pAirRefValue = 0.0;

    setRefCell
    (
        pAir,
        mesh.solutionDict().subDict("PISO"),
        pAirRefCell,
        pAirRefValue
    );

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

    volScalarField Cp
    (
        IOobject
        (
            "Cp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        CpAir,
        zeroGradientFvPatchScalarField::typeName
    );

    surfaceScalarField CpPhiAir
    (
        IOobject
        (
            "CpPhiAir",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "CpPhiAir",
            dimVelocity*dimDensity*dimSpecificHeatCapacity*dimArea,
            0
        )
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

    surfaceScalarField rhoTFlux
    (
        IOobject
        (
            "rhoTFlux",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh.magSf()*fvc::snGrad(T)
    );


//    Info<< "Reading field nPart\n" << endl;
//    volScalarField nPart
//    (
//        IOobject
//        (
//            "nPart",
//            runTime.timeName(),
//            mesh,
//            IOobject::MUST_READ,
//            IOobject::AUTO_WRITE
//        ),
//        mesh
//    );

//new fields for calculation of evaporation/condensation kernels
//saturation ratio
    volScalarField S
    (
        IOobject
        (
            "S",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
//saturation vapor pressure
    volScalarField pSat
    (
        IOobject
        (
            "pSat",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
//nondimensional surface tension
    volScalarField theta
    (
        IOobject
        (
            "theta",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

//nucleation rate
    volScalarField Jk
    (
        IOobject
        (
            "Jk",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField Jk2
    (
        IOobject
        (
            "Jk2",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

//surface tension field
    volScalarField sigma
    (
        IOobject
        (
            "sigma",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
//critical cluster size
    volScalarField dpstar
    (
        IOobject
        (
            "dpstar",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
//critical cluster size (volume)
    volScalarField vpstar
    (
        IOobject
        (
            "vpstar",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
//critical cluster size node
    volScalarField kstar
    (
        IOobject
        (
            "kstar",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

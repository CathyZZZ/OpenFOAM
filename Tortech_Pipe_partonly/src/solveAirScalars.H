{

        Info << nl << "Solving air scalars" << endl;
//Info<<nl<<"XXXXX"<< Nair[0][1] << Nair[0][5000]<<endl;
   
        forAll(dp, i) {
        Info << "diameter for node " << i << " " << Nair[i][0] << " = " << dp[i].value() << " " <<endl;
}
/*forAll(Nair[0], cellI) {
    Info << Nair[0][cellI]<< nl<<endl;
}
*/
 //       rhoTFlux = mesh.magSf()*fvc::snGrad(T);//-3*mu*rhoAir*rhoTFlux/4/T
volScalarField& N1 = Nair[1];
Info << nl << "N1[400]="<< N1[400] << nl << "N1[800]=" << N1[800] << nl << "N1[5000]"<< N1[5000] << endl; 

Info << nl << "T[400]="<< T[400] << nl << "T[800]=" << T[800] << nl << "T[5000]"<< T[5000] << endl; 
Info << nl << "S" << S[400] << nl << "S" << S[800] << nl << "S" << S[5000] << endl;

volScalarField term1 = S*0.0;
volScalarField term2 = S*0.0; //Foam::pow((2*sigma)/(pi*m1),0.5);
volScalarField term3 = S*0.0; //Foam::log(S);
volScalarField term4 = S*0.0;
volScalarField term5 = S*0.0; //Foam::exp(term4);
volScalarField term6 = S*0.0; //2*theta/(3*term3);
volScalarField term7 = S*0.0;
volScalarField SN1 = S*0.0; //kb*T/pSat;
volScalarField JkN1 = S*0.0; //SN1*v1*term1*term2*term5;
forAll(N1,cellI) {
    sigma[cellI] = (ACoeff.value()-BCoeff.value()*T[cellI])/1000;
    pSat[cellI] = pAir[cellI]*Foam::exp(CCoeff.value()-DCoeff.value()/T[cellI]);
    S[cellI] = N1[cellI]* kb * T[cellI]/pSat[cellI];
    theta[cellI] = s1.value()*sigma[cellI]/(T[cellI] * kb );
    term1[cellI] = Foam::pow((pSat[cellI]/( kb * T[cellI])),2);
    term2[cellI] = Foam::pow((2*sigma[cellI])/(pi*m1.value()),0.5);
    term3[cellI] = Foam::log(S[cellI]+0.01);
    term4[cellI] = theta[cellI] - (4*Foam::pow(theta[cellI],3))/(27*Foam::pow(term3[cellI],2));
    term5[cellI] = Foam::exp(term4[cellI]);
    term6[cellI] = 2*theta[cellI]/(3*term3[cellI]);
    SN1[cellI] = kb * T[cellI]/pSat[cellI];//N1 at saturated condition
    JkN1[cellI] = SN1[cellI]*v1.value()*term1[cellI]*term2[cellI]*term5[cellI];
    Jk[cellI] = S[cellI]*v1.value()*term1[cellI]*term2[cellI]*term5[cellI];
    kstar[cellI] = Foam::pow(term6[cellI],3);
    dpstar[cellI] = 4*sigma[cellI]*v1.value()/(kb * T[cellI] * term3[cellI]);
    vpstar[cellI] = pi/6*Foam::pow(dpstar[cellI],3);
 //   Info << "RRRRRRRRRRRR"<< endl;
//    Info << sigma[cellI]<<", "<< term6[cellI]<<", "<< vpstar[cellI]<<endl;
}
Info << "dpstar" << dpstar[5000]<< endl;
/*
Info << nl << "pAir[400]="<< pAir[400] << nl << "pAir[800]=" << pAir[800] << nl << "pAir[5000]"<< pAir[5000] << endl;
Info << nl << "sigma[400]="<< sigma[400] << nl << "sigma[800]=" << sigma[800] << nl << "sigma[5000]"<< sigma[5000] << endl;
Info << nl << "pSat[400]="<< pSat[400] << nl << "pSat[800]=" << pSat[800] << nl << "pSat[1500]"<< pSat[1500] <<nl<< "pSat[5000]="<< pSat[5000]<< endl;
Info <<" s1"<< s1<<nl<<"kb"<<kb <<"m1"<<m1 << endl;
Info << nl << "S[5000]="<< S[5000] << nl << "theta[5000]=" << theta[5000] << nl << "term1[5000]"<< term1[5000] <<nl<< "term2[5000]="<< term2[5000]<< endl;

Info << nl << "term3[5000]="<< term3[5000] << nl << "term4[5000]=" << term4[5000] << nl << "term5[5000]"<< term6[5000] <<nl<< "SN1[5000]="<< SN1[5000]<< endl;

Info << nl << "JkN1[5000]="<< JkN1[5000] << nl << "Jk[5000]=" << Jk[5000] << nl << "kstar[5000]"<< kstar[5000] <<nl<< "dpstar[5000]="<< dpstar[5000]<< endl;

Info << nl << "Jk[400]="<< Jk[400] << nl << "Jk[800]=" << Jk[800] << nl << "Jk[5000]"<< Jk[5000] << endl;
Info << nl << "kstar[400]="<< kstar[400] << nl << "kstar[800]=" << kstar[800] << nl << "kstar[5000]"<< kstar[5000] << endl;

Info << nl << "dpstar[400]="<< dpstar[400] << nl << "dpstar[800]=" << dpstar[800] << nl << "dpstar[5000]"<< dpstar[5000] << endl;
Info << nl << "vpstar[400]="<< vpstar[400] << nl << "vpstar[800]=" << vpstar[800] << nl << "vpstar[5000]"<< vpstar[5000] << endl;
*/
//collision frequency beta calculation
forAll(term7,cellI) {
    term7[cellI] = Foam::pow((3.0/( pi * 4.0)),0.166666) * Foam::pow(( 6* kb * T[cellI] /rhoDisp.value()), 0.5);
}
//Info << nl << "term7[400]="<< term7[400] << nl << "term7[800]=" << term7[800] << nl << "term7[1500]"<< term7[1500] << endl;
//    Info << "MW"<< MWDisp << endl;
//saturation concentration over particles of size i
    //Info << T<< endal;
 
for (int k =0; k<node; k++) {
    forAll( NSatk[k], cellI) {       
         NSatk[k][cellI] =  pSat[cellI]/ (kb *T[cellI]) *Foam::exp(4000 * sigma[cellI] * MWDisp.value()/(R.value() * T[cellI] * rhoDisp.value() * dp[k].value()));        
}

//    Info << "MW"<< MWDisp << endl;
  //       Info <<" k=" <<k << endl;

//Info <<"NSatl["<<k<<"]="<< NSatk[k][5000] << endl;
/*
double aa= pSat[5000] / (kb * T[5000]) ;
double bb = 4 * sigma[5000] * MWDisp.value() / (R.value() * T[5000] * rhoDisp.value() * dp[k].value());
double cc = Foam::exp(bb);
double dd = aa * cc;
*/
/*
Info << nl << "NSatk["<<k<<"][5000]="<< NSatk[k][5000] << endl;
//Info << aa << nl << bb << nl << cc << dd <<endl;
Info <<"NSatk_total"<< pSat[5000] / (kb * T[5000]) * Foam::exp(4 * sigma[5000] * MWDisp.value() / (R.value() * T[5000] * rhoDisp.value() * dp[k].value())) <<endl;
//Info << "1="<< dp[k].value() << endl;*/
/*Info << "pSat=" <<pSat[5000] <<nl<<"kb="<< kb <<nl<<"T="<< T[5000] <<nl<<"sigma="<< sigma[5000] <<"MWDisp="<< MWDisp.value() <<nl<<"R="<< R.value() <<nl<<"rhoDisp="<< rhoDisp.value() <<nl<<"dp[k]"<< dp[k].value() << endl;
*/
}


//scalar beta = 0;
 //   Info << "XXXXXXXx"<< endl;
//collision frequency kernels used for coagulation and growth
for (label j = 1; j < node; j++) {
 //   Info<<"j="<< j<< endl;
    for (label i = 1; i < node; i++) {
        forAll(Beta[j][i], cellI) {
             scalar tem1 = 1/vp[i].value() + 1/vp[j].value();
             scalar beta = Foam::pow(vp[j].value(), 0.33333)+Foam::pow(vp[i].value(),0.33333);
//            Info<<"qqqq"<< endl;
            Beta[j][i][cellI] = dimCorr.value() *term7[cellI] *Foam::pow(tem1,0.5)*Foam::pow(beta,2.0);
        }
//Info<<"1111"<< endl;
        //        double tem2= Foam::pow(1/vp[i].value()+1/vp[j].value(),0.5);
//        Info <<" temporary_term["<<j<<"]["<<i<< "]="<< beta << endl;
//        Info << "dimCorr="<<dimCorr<<endl;
//        Info <<" temporary_term2[" << j << "][" <<i << "]="<< tem2 << endl;
//        Info << "Beta["<<j<<"]["<<i<<"]="<< Beta[j][i][5000] << endl;
                                                                                                                    
//        Info << "Beta[" << i << "][ " << j <<"]="<< Beta[i][j][5000] << endl;
       // Info <<"size_ofBeta"<< Beta.size() << endl;
       // Info <<"size_ofBeta[0]"<< Beta[0].size() << endl;
    }
}
//Info<<"1111111111"<< endl;

//Info << N1<< endl;
volScalarField srcMonoNuc = 0.0 * Jk/N1 ;
volScalarField srcMonoPlus = 0.0*Jk/N1;
//volScalarField srcMonoMinus = 0.0 *Jk;
//Info << Jk<< endl;
//Info << N1<< endl;
//Info<<srcMonoNuc<< endl;
volScalarField srcNucleation = 0.0 * Jk;
 //   Info<<"22222222"<< endl;

volScalarField tem_growth = Nair[1];
//    Info<<"22222222"<< endl;
volScalarField t1Plus = 0.0 * Jk;
volScalarField t2Plus = 0.0 * Jk;
volScalarField t3Plus = 0.0 * Jk;
volScalarField t4Plus = 0.0 * Jk;


for (label k = 1; k<node; k++) {
    forAll(Nair[k], cellI) {
        Nairold[k][cellI] = Nair[k][cellI];
    }
}
//SOURCE TERM CALCULATION

for (label k = 1; k<node; k++) {
    volScalarField& Nk = Nair[k];
    volScalarField srcCoagulation = 0.0 *Jk;
    volScalarField srcgrowth = 0.0*Jk;
    volScalarField addterm = 0.0 *Jk;
    volScalarField subterm = 0.0 *Jk2;
    volScalarField coagterm1 = 0.0 *Jk;//addterm
   // Info<<"N1"<< N1<< nl<<endl;
    volScalarField coagterm2 = 0.0 *Jk/N1;//sub term
    volScalarField coagterm3 = 0.0 *Jk;
//nucleation source terms
//Info<<"33333"<<endl;
    if (nucleation) {
        forAll (Nzeta[k], cellI) {
            if (vpstar[cellI] <= v1.value()) {
              Nzeta[1][cellI] = vpstar[cellI] / vp[1].value();
            }
            else if (( vp[k-1].value() <= vpstar[cellI] && (vpstar[cellI] <= vp[k].value()))) {
            Nzeta[k][cellI] = vpstar[cellI] / vp[k].value();
            }
            else 
                Nzeta[k][cellI] = 0;
            
            srcNucleation[cellI] = Jk[cellI] * Nzeta[k][cellI];
            srcMonoNuc[cellI] = -Jk[cellI] * kstar[k];//different from Robert's
       }

        Info << nl << "srcNucleation[400]="<< srcNucleation[400] << "srcNucleation[800]=" << srcNucleation[800]  << "srcNucleation[1500]"<< srcNucleation[1500] << endl;
    }

    if (coagulation) {
        forAll (Nk, cellI) {
            for ( label i=1; i < node; i++) {
                coagterm2[cellI] += Beta[k][i][cellI] * Nairold[i][cellI];
                for ( label j=1; j <=k; j++) {
                coagterm1[cellI] += Xijk[k][j][i].value() * Beta[j][i][cellI] * Nairold[i][cellI] * Nairold[j][cellI];
               
                }
            }
        }
 //       Info << "Beta[1][1][5000]"<<Beta[1][1][5000]<<endl;
 //       Info << "Beta[1][40][5000]"<<Beta[1][40][5000]<<endl;
 //       Info << "Nair[1][5000]" << Nair[1][5000]<< endl;
        forAll(srcCoagulation, cellI) {

            srcCoagulation[cellI] = 0.5 * coagterm1[cellI];
        }
        Info<<"coagterm1[5000]"<<coagterm1[5000]<< endl;
        Info<<"coagterm2[5000]"<<coagterm2[5000]*Nk[5000]<< endl;
        Info << nl << "srcCoagulation[400]="<< srcCoagulation[400] << "srcCoagulation[800]=" << srcCoagulation[800] << "srcCoagulation[1500]"<< srcCoagulation[1500] << endl;
    }

    if (growth) {
        if (k < Nair.size()-1) {
            //double oem_growth[cellI] = Nair[1][cellI];

            forAll ( N1, cellI) {
                if (N1[cellI] > NSatk[k-1][cellI]) {
                    if (k==1) {
                        addterm[cellI] = 0;
                    }
                    else
                     
                        addterm[cellI] = v1.value() /(vp[k].value() - vp[k-1].value()) * Beta[0][k-1][cellI]* (N1[cellI] -NSatk[k-1][cellI]) * tem_growth[cellI];
                        t1Plus[cellI] = t1Plus[cellI] -Beta[0][k-1][cellI] *(N1[cellI] - NSatk[k-1][cellI]) * Nair[k-1][cellI];
                        //addterm[cellI] = v1.value()/(vp[k].value() - vp[k-1].value()) * Beta[0][k-1][cellI] (N1[cellI] -NSatk[k-1][cellI]) * Nair[k-1][cellI];
                }
                if (N1[cellI]<NSatk[k+1][cellI]) {
                    addterm[cellI] = -(v1.value() /( vp[k+1].value() - vp[k].value())) * Beta[0][k+1][cellI] *(N1[cellI] - NSatk[k+1][cellI]) * Nair[k+1][cellI];
                    t2Plus[cellI] = t2Plus[cellI] + Beta[0][k+1][cellI] *(N1[cellI] - NSatk[k+1][cellI]) * Nair[k+1][cellI];
                }
                if (N1[cellI] < NSatk[k][cellI]) {
                    subterm[cellI] = -(v1.value() / (vp[k].value() - vp[k-1].value())) * Beta[0][k][cellI] * (N1[cellI] - NSatk[k][cellI]) * Nair[k][cellI];
                    t3Plus[cellI] =t3Plus[cellI] + Beta[0][k][cellI] *(N1[cellI] - NSatk[k][cellI]) * Nair[k][cellI];
                }

                if (N1[cellI] > NSatk[k][cellI]) {
                    subterm[cellI] = (v1.value() / (vp[k+1].value() - vp[k].value())) * Beta[0][k][cellI] * (N1[cellI] - NSatk[k][cellI]) * Nair[k][cellI];
                    t4Plus[cellI] = t4Plus[cellI] - Beta[0][k][cellI] *(N1[cellI] - NSatk[k][cellI]) * Nair[k][cellI];
                
                }
                if (k != 1) {
                    tem_growth[cellI] = Nair[k-1][cellI]; 
                }         
                srcgrowth[cellI] = addterm[cellI] + subterm[cellI];
            }
        }
    
        Info << nl << "srcgrowth[400]="<< srcgrowth[400] << "srcgrowth[800]=" << srcgrowth[800] << "srcgrowth[1500]"<< srcgrowth[1500] << endl;
    }
/*
Info<<"Dimension1"<< srcNucleation.dimensions() << endl;
Info<<"Dimension2"<< srcCoagulation.dimensions() << endl;
Info<<"Dimension3"<< srcgrowth.dimensions() << endl;
Info<<"Dimension4"<< phiAir.dimensions() << endl;
Info<<"Dimension5"<< Nk.dimensions() << endl;
*/



    solve
        (
         fvm::div(phiAir, Nk, "div(phi,n)")
         + fvm::Sp(rho * coagterm2, Nk)
         - rho * srcCoagulation
         - rho * srcNucleation
         - rho * srcgrowth,
         mesh.solver("nPart")
        );
    Nk.max(Foam::doubleScalarSMALL);
}


 
    forAll(srcMonoPlus, cellI) {
        srcMonoPlus[cellI] = t1Plus[cellI] + t2Plus[cellI] + t3Plus[cellI] + t4Plus[cellI];
    }
    
    if (monomerBal)
    {
        volScalarField N1old = N1;
        fvScalarMatrix monoBalEqn(
                fvm::div(phiAir,N1,"div(phi,n)")
                - fvm::Sp(rho * srcMonoNuc, N1)
 //               - rho * srcMonoPlus
                );
        monoBalEqn.relax(0.99);
        monoBalEqn.solve(mesh.solver("nIonsN"));

    }
    
    Info<<"N1"<< N1<< nl<<endl;
    scalar sumOfVolumes = 0;
    scalar sumOfVolumes2 = 0;
    for (label i=0; i<node; i++) {
        sumOfVolumes += vp[i].value() * Nair[i][0];
        sumOfVolumes2 += vp[i].value() * Nair[i][5000];
    }
    Info <<nl<<"sumOfVolumes=" << sumOfVolumes << endl;
    Info <<nl<<"sumOfVolumes2="<< sumOfVolumes2 << endl;
    }

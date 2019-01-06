#include <boost/math/distributions/normal.hpp>
#include <iostream> // input/output functions
#include <string>   // strings
#include <cmath>	// includes exp and log
#include <fstream>
#include <sstream>

using namespace std; 
 
double* linspace(double, double, int); 
double* concentration(double, double, int, double*, long);
 
int main() 
{
	string line;
	size_t pos;
	string numStr;
	int p1;
	int p2;
	double CMDValue = 0;
	double stdevValue = 0;
	double NtotValue = 0;
	int NbinValue = 0;

  	ifstream myfile ("./run/constant/transportProperties"); //stream transportProperties File
  	if (myfile.is_open())
  	{
    	while ( getline (myfile,line) )
    	{
			pos = line.find("Nbin"); //search for keyword "rho" and save the respective line in string "line"
			if(pos!=string::npos) // string::npos is returned if string is not found
			{
				p1 = line.find("]"); //position of "]"-character
				p2 = line.find(";"); //position of the ";"-character
				numStr = line.substr(p1+1,p2-p1-1); //extract string between ] and ;
				stringstream streamValue(numStr); //get rid of space before number
				streamValue >> NbinValue;
				
//				cout << "value of Nbin: " << NbinValue << '\n';
			}
			pos = line.find("CMD"); 
			if(pos!=string::npos)
			{
				p1 = line.find("]");
				p2 = line.find(";");
				numStr = line.substr(p1+1,p2-p1-1);
				stringstream streamValue(numStr);
				streamValue >> CMDValue;
				
//				cout << "value of CMD: " << CMDValue << '\n';
			}
			pos = line.find("stdev"); 
			if(pos!=string::npos)
			{
				p1 = line.find("]");
				p2 = line.find(";");
				numStr = line.substr(p1+1,p2-p1-1);
				stringstream streamValue(numStr);
				streamValue >> stdevValue;
				
//				cout << "value of stdev: " << stdevValue << '\n';
			}
			pos = line.find("Ntot"); 
			if(pos!=string::npos)
			{
				p1 = line.find("]");
				p2 = line.find(";");
				numStr = line.substr(p1+1,p2-p1-1);
				stringstream streamValue(numStr);
				streamValue >> NtotValue;
				
//				cout << "value of Ntot: " << NtotValue << '\n';
			}
    	}
    myfile.close();
  	}

  	else cout << "Unable to open file"; 


   double m = log(CMDValue); //Mean, mediam, mode of the normal [m]
   double stdv = log(stdevValue);	//Standard deviation of the normal[-]
   long Ntot = NtotValue;	//Particles/m^3
   int Nbin = NbinValue;	//Number of bins

   if (Nbin == 1){stdv = log(1.00000001);}

   boost::math::normal_distribution<double> normal(m, stdv);	//Creates the distribution
   double a = quantile(normal, 0.0001);		//Lower bound (0.01 percentile)
   double b = quantile(normal, 0.9999);	//Upper bound (99.99 percentile)
   double *d2 = linspace (a, b, Nbin);		//FIRST FUNCTION
   double *conc= concentration(m, stdv, Nbin, d2, Ntot);	//SECOND FUNCTION 
   double d[Nbin+1];
   
   for (int i=0; i<=Nbin; i++){
	   d[i] = exp(d2[i]);
   }

   for (int i=0; i<Nbin; i++){
	   cout << round((d[i]+d[i+1])/4e-9)*1e-9  << "	" << conc[i] <<endl;
   }
   
   return 0;
}

double* linspace (double a, double b, int Nbin){
	
	double c = (b-a)/Nbin;
	double* dtem = new double[Nbin+1];	//d:(Nbin+1) bounds.  d[0]=a and d[N]=b
	for (int i=0; i<=Nbin; i++) {
		dtem[i]=a+i*c;
	}

	for  (int i=2; i<=Nbin-2; i++) {
		if (dtem[i+1] == dtem[i]) {
			printf("The bin width collapsed to 0");
			exit (EXIT_FAILURE);
			}
	}
	
	return dtem;
	delete[] dtem;
}


double* concentration (double m, double stdv, int Nbin, double* d2, long Ntot) {

    boost::math::normal_distribution<double> normal(m, stdv); //Distribution
	double* conctem = new double[Nbin];
	
	for (int i=0; i<Nbin; i++){
		conctem[i]=Ntot*(cdf(normal, d2[i+1])-cdf(normal,d2[i]));	// Concentration of each bin
	}
	
	return conctem;
	delete[] conctem;
}

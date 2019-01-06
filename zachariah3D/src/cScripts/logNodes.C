//#include <boost/math/distributions/normal.hpp>
#include <iostream> // input/output functions
#include <string>   // strings
#include <cmath>	// includes exp and log
#include <fstream>
#include <sstream>

using namespace std; 
 
double* logspace(double, double, int); 
double* concentration(double, double*, double, double, int);
 
int main() 
{
	string line;
	size_t pos;
	string numStr;
	int p1;
	int p2;
	double dpMaxValue = 0;
	int nodesValue = 0;
	double rhoDispValue = 0;
	double MWDispValue = 0;
	double CCoeffValue = 0;
	double DCoeffValue = 0;
	double NpDispValue = 0;
	double dpDispValue = 0;
	double TinitValue = 0;

  	ifstream myfile ("./constant/transportProperties"); //stream transportProperties File
  	if (myfile.is_open())
  	{
    	while ( getline (myfile,line) )
    	{
			pos = line.find("rhoDisp"); //search for keyword "rho" and save the respective line in string "line"
			if(pos!=string::npos) // string::npos is returned if string is not found
			{
				p1 = line.find("]"); //position of "]"-character
				p2 = line.find(";"); //position of the ";"-character
				numStr = line.substr(p1+1,p2-p1-1); //extract string between ] and ;
				stringstream streamValue(numStr); //get rid of space before number
				streamValue >> rhoDispValue;
			}
			pos = line.find("dpMax"); 
			if(pos!=string::npos)
			{
				p1 = line.find("]");
				p2 = line.find(";");
				numStr = line.substr(p1+1,p2-p1-1);
				stringstream streamValue(numStr);
				streamValue >> dpMaxValue;
			}
			pos = line.find("nodes"); 
			if(pos!=string::npos)
			{
				p1 = line.find("]");
				p2 = line.find(";");
				numStr = line.substr(p1+1,p2-p1-1);
				stringstream streamValue(numStr);
				streamValue >> nodesValue;
			}
			pos = line.find("MWDisp"); 
			if(pos!=string::npos)
			{
				p1 = line.find("]");
				p2 = line.find(";");
				numStr = line.substr(p1+1,p2-p1-1);
				stringstream streamValue(numStr);
				streamValue >> MWDispValue;
			}
			pos = line.find("CCoeff"); 
			if(pos!=string::npos)
			{
				p1 = line.find("]");
				p2 = line.find(";");
				numStr = line.substr(p1+1,p2-p1-1);
				stringstream streamValue(numStr);
				streamValue >> CCoeffValue;
			}
			pos = line.find("DCoeff"); 
			if(pos!=string::npos)
			{
				p1 = line.find("]");
				p2 = line.find(";");
				numStr = line.substr(p1+1,p2-p1-1);
				stringstream streamValue(numStr);
				streamValue >> DCoeffValue;
			}
			pos = line.find("NpDisp"); 
			if(pos!=string::npos)
			{
				p1 = line.find("]");
				p2 = line.find(";");
				numStr = line.substr(p1+1,p2-p1-1);
				stringstream streamValue(numStr);
				streamValue >> NpDispValue;
			}
			pos = line.find("dpDisp"); 
			if(pos!=string::npos)
			{
				p1 = line.find("]");
				p2 = line.find(";");
				numStr = line.substr(p1+1,p2-p1-1);
				stringstream streamValue(numStr);
				streamValue >> dpDispValue;
			}
			pos = line.find("Tinit"); 
			if(pos!=string::npos)
			{
				p1 = line.find("]");
				p2 = line.find(";");
				numStr = line.substr(p1+1,p2-p1-1);
				stringstream streamValue(numStr);
				streamValue >> TinitValue;
			}
    	}
    myfile.close();
  	}

  	else cout << "Unable to open file"; 
	

	double v1 = MWDispValue/rhoDispValue/6.022e23;
	double d1 = pow(6/3.14159*v1,0.33333333333);
	double d2 = dpMaxValue;
	double psat = 1.01325e5*exp(CCoeffValue-DCoeffValue/TinitValue);
	double N1 = 1.001*psat/(1.38e-23*TinitValue);
   	double *dVec = logspace (d1, d2, nodesValue);	//node spacing
	double *conc= concentration(N1, dVec, NpDispValue, dpDispValue, nodesValue);

   	for (int i=0; i<nodesValue; i++){
//	   	cout << "node" << i+1 << ": size " << dVec[i] << ", concentration " << conc[i] <<endl;
	   	cout << dVec[i] << " " <<conc[i] <<endl;
   	}

	return 0;
}

double* logspace (double dmin, double dmax, int nodes){
	
	double logMin = log10(dmin);
	double logMax = log10(dmax);
	double delta = (logMax - logMin)/(nodes-1);
	double* v = new double[nodes];
	
	double exponent = 0;

	for (int i=0; i< nodes; i++)
	{
		exponent = logMin+i*delta;
		v[i] = pow(10,exponent);
	}	
	return v;
	delete[] v;
}




double* concentration (double N1, double* dVec, double Np, double dp, int nodes ) {

	int minPos=0;
	double* conctem = new double[nodes];	
	double* abstem = new double[nodes];
	for (int i=0; i< nodes; i++)
	{
		conctem[i]=0;
		abstem[i]=abs(dVec[i]-dp);
		if(abstem[i]<abstem[minPos]) {
			minPos = i;
		}
	}
	
	conctem[0] = N1;
	conctem[minPos-1] = Np;
	
	return conctem;
	delete[] conctem;
}





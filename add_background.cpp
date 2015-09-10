//This program will accept an XPS and spectral function file and use it to generate an XPS with background added
//Jonathan Curtis
//09/07/15

#include <iostream>
#include <fstream>
#include <sstream>
#include <complex>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <ctime>
using namespace std;

//First we define some constants 
const double pi = 3.1415926535897931;	//pi to high precision
const double hartree = 27.21138505;		//1 Hartree, in eV
const double AUtime	= .02418884326505; 	//The atomic unit of time, in fs, to high precision
const double hbar = 4.135667516;		//hbar, in eV*fs, to high precision
const complex<double> I(0.0,1.0);		//The imaginary unit, to high precision

int main(int argc, char * argv[]){
	//First we check the implementation
	//This implementation will be an input of the XPS file and the spectral function file and then the output file name
	if(argc!=4){
		cout<<"Incorrect implementation. Usage is: 'program_name' 'XPS.dat' 'spec.dat' 'out.dat'"<<endl;
		return 0;
	}

	//We now parse the input for the appropriate files
	const char * XPSIN = argv[1];
	const char * SPECIN = argv[2];
	const char * OUTFILE = argv[3];

	//Now we open the files and read in the XPS. We then center it. 
	//These variables will hold the frequencies and the XPS function
	int num_xps_w = 0;	//The number of frequency steps in the XPS file
	//We auto-detect this by checking the size of the XPS file 

	string line;	//We will read each line into this variable 
	ifstream XPScount(XPSIN);	//We open an ifstream for counting the number of lines
	while(getline(XPScount,line)){
		//We incremement the counter for each line 
		num_xps_w++;	//Increment
	}
	XPScount.close();	//Close the file 

	cout<<"num_xps_w = "<<num_xps_w<<endl;

	//Now we allocate the frequencies and XPS arrays
	double * xps_ws = new double[num_xps_w];	//Allocate the frequencies
	double * xps = new double[num_xps_w];		//Allocate the xps

	//And now we go through the file and fill the arrays

	ifstream XPSfile(XPSIN);
	int line_count = 0;	//This will count the line so we can access the appropriate array element to fill
	while(getline(XPSfile,line)){
		//We need three tokens
		//The frequencies are the first token per line and the XPS is the third token per line 
		string token1,token2,token3;	//Our tokens
		istringstream iss(line);	//This tokenizes each line
		iss>>token1>>token2>>token3;	//This extracts each token

		//The frequencies will be in eV and the XPS will be in 1/eV units 
		xps_ws[line_count] = atof( token1.c_str() );	//We convert the first token to a double and fill the array
		xps[line_count] = atof( token3.c_str() );	//We convert the third token and fill the array

		line_count++;	//Increment the line counter 
	}
	XPSfile.close();

	double * cum_xps = new double[num_xps_w];	//This is the cumulative distribution for the spectral function
	//It has as many entries as xps does
	//It is related to the XPS by
	//F(w) = 1 - integral_{-w}^{infty} XPS(w')dw'
	//We need B(w) = integral_{w}^{infty} XPS(w')dw'
	//Thus we have B(w) = 1-F(-w)
	//Because XPS and cumulative function are output on the same array (up to a minus sign)
	//We read it in and fill it with F[i] in reverse to get 1-C(w) = F(-w). 
	//We use as a background term 1-F(-w) = integral_{w}^{infty}A(w')dw'

	//And now we read in the spectral function array 
	ifstream SPECfile(SPECIN);
	line_count = 0;	//Re-zero the line-counter
	while(getline(SPECfile,line)){
		//We need three tokens
		//We just want the third token
		string token1,token2,token3;
		istringstream iss(line);	//Tokenize the line 
		iss>>token1>>token2>>token3;
		//The cumulative function has units that are unitless
		cum_xps[num_xps_w-1-line_count] = atof( token3.c_str() );	//We convert the token and fill the array IN REVERSE ORDER
		line_count++;	//Increment the counter
	}

	const double background_ratio = .1;	//This is the relative ratio of the background term
	cout<<"background ratio: "<<background_ratio<<endl;

	//We use a signal of 
	//xps(w) + background_ratio*B(w)

	//Now we print out the data with the background term added on 
	//Open an ofstream
	ofstream OUTfile(OUTFILE);
	for(int i = 0; i < num_xps_w; i++){
		OUTfile<<xps_ws[i]<<" "<<xps[i]+background_ratio*(cum_xps[i])<<" "<<xps[i]<<" "<<background_ratio*(cum_xps[i])<<endl;
	}
	OUTfile.close();

	delete [] xps_ws;
	delete [] xps;
	delete [] cum_xps;
		
	return 0;
}
//This is a simple script designed to compute the spectral function from a core response
//Jonathan Curtis
//08/05/15
//This is a rather inelegant and hacked together program, but it is hopefully quicker to debug and write

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


//This method will check to see if label==reference, and if it does will assign variable = int(value)
void parse_int(string label, string value, const string reference, int& variable){
	if(label==reference){
		variable = int(atof(value.c_str()));
	}
}

//This method will check to see if label==reference, and if it does will assign variable = double(value)
void parse_doub(string label, string value, const string reference, double& variable){
	if(label==reference){
		variable = double(atof(value.c_str()));
	}
}

//This method will check to see if label==reference, and if it does will assign variable = value
void parse_str(string label, string value, const string reference, string& variable){
	if(label==reference){
		variable = value;
	}
}

//This function will snap a given point onto a grid 
//It will return the largest array index that has a value less than the passed value
//If the element is out of bounds, it will return 0 or arrysize-1 depending on which end of the array it crosses
//It assumes that the array is uniformly spaced
int snapToGrid(double x, int size, double * arry){
  //First we deal with the case that x<values[0]
  if(x<arry[0]){
    return 0;
  }
  //Next, if it is in the array range
  else if(x<arry[size-1]){
    return floor( (x-arry[0])/(arry[1]-arry[0]) ); //We assume the array is uniform and has spacing values[1]-values[0]
  }
  else{
    return size-1;  //It is larger than the largest entry so we return the largest index
  }
}

int main(int argc, char * argv[]){
	//Analysis variables and parameters
	int num_core_steps = 0;	//The number of times steps in our core response
	int num_beta_steps = 0;	//The number of frequency steps we use to compute the beta function
	int num_clnt_steps = 0;	//The number of time steps we use to compute the cumulant
	int num_xps_steps = 0;	//The number of frequency steps we use to compute the XPS

	double * core_times = NULL;	//The actual times in the core response
	double * beta_freqs = NULL;	//The frequencies array for beta
	double * clnt_times = NULL;	//The times for the cumulant 
	double * xps_freqs = NULL;		//The frequencies for the xps used
	double * spec_freqs = NULL;		//The array of frequencies for the spectral function

	double * core = NULL;		//The core response values
	double * beta = NULL;			//The beta values
	complex<double> * clnt = NULL;	//The cumulant values
	complex<double> * xps = NULL;	//The xps spectrum
	double * spec_den = NULL;	//The spectral function density
	double * spec_cum = NULL;	//The spectral function cumulative distribution

	double core_uniform_width = 0.0;	//The width of the uniform widening used on the core data. Default of 0.0
	double core_uniform_order = 2.0;	//The order of the uniform widening used on the core data. Default of 2.0
	double core_high_w_width = 0.0;		//The width of the high frequency broadening applied to the data. Default of 0.0
	double core_high_w_order = 2.0;		//The order of the high frequency broadening applied to the data. Default of 2.0

	double max_t_core = 0.0;	//The last value of the time at which the core data is sampled.
	double max_w_beta = 0.0;	//The maximum frequency we compute beta to. Default of zero will be changed later to the cutoff
	double max_t_clnt = 0.0;	//The maximum time we compute the cumulant to. Default of zero will be changed later to the cutoff
	double max_w_xps = 0.0;		//The maximum frequency we compute the XPS for. We compute to this value symmetrically about zero. Default of zero will be changed later to the cutoff

	double xps_width = 0.0;	//The size of Gaussian widening applied to the xps spectrum
	double xps_order = 2.0;		//The order of the damping used on the XPS integral. Default of 2.0

	double delta_t_core = 0.0;	//The spacing ofthe time steps for the core data
	double cutoff_w_beta = 0.0;	//The cutoff frequency for the beta frequency, based on .666 Nyquist frequency
	double delta_w_beta = 0.0;	//The frequency spacing for beta
	double cutoff_t_clnt = 0.0;	//The cutoff time for the cumulant, based on the the sampling interval for frequency as .666 Nyquist time
	double delta_t_clnt = 0.0;	//The time spacing for the cumulant times
	double cutoff_w_xps = 0.0;	//The cutoff frequency for the xps, based on .666 the Nyquist frequency from sampling the cumulant
	double delta_w_xps = 0.0;		//The frequency spacing for the xps

	double coreav1 = 0.0;	//This will be used to compute the first average of the core response
	double coreav2 = 0.0;	//The average response after the weigthing is done

	double spec_norm = 0.0;	//This computes the integral of the xps and checks to see if it is close to one
	double quasi_peak_weight = 0.0;	//This is the weight of the first quasi particle excitation peak. It is the integral of beta(w)/w^2
	double quasi_peak_shift = 0.0;	//This is the quasi-particle peak shift, given by the integral of beta(w)/w
	double sat_peak_weight = 0.0;	//This computes the weight of the satellite  peaks
	double spec_peak_value = 0.0; //This is the global maximum of the spectral function
	double spec_peak_point = 0.0;	//This will be used to record the location of the maxima of the spectral function

	//File names
	string COREIN;
	string COREOUT;
	string BETAOUT;
	string CLNTOUT;
	string XPSOUT;
	string SPECOUT;
	
	cout<<"--------------------------------------------------"<<endl;
	//First we load in the data from a parameter file passed by command line arg
	if(argc!=2){
		cout<<"Implementation error: program accepts one argument, the parameter file name."<<endl;
		cout<<"--------------------------------------------------"<<endl;

		return 0;
	}

	//We proceed by opening the passed name
	//Then we read it line by line and parse for the parameters
	cout<<"Loading parameters..."<<endl;
	ifstream paramfile(argv[1]);
	string paramline;
	while(getline(paramfile,paramline)){
		string label, value;	//These will be used to parse for a label and extract value
		istringstream iss(paramline);	//This is used because it is good at parsing strings for values

		iss>>label>>value;	//Extract the first entry on each line into label and the next entry into value
		parse_str(label,value,"COREIN",COREIN);
		parse_str(label,value,"COREOUT",COREOUT);
		parse_str(label,value,"BETAOUT",BETAOUT);
		parse_str(label,value,"CLNTOUT",CLNTOUT);
		parse_str(label,value,"XPSOUT",XPSOUT);
		parse_str(label,value,"SPECOUT",SPECOUT);

		parse_int(label,value,"num_beta_steps",num_beta_steps);
		parse_int(label,value,"num_clnt_steps",num_clnt_steps);
		parse_int(label,value,"num_xps_steps",num_xps_steps);

		parse_doub(label,value,"core_uniform_width",core_uniform_width);
		parse_doub(label,value,"core_uniform_order",core_uniform_order);
		parse_doub(label,value,"core_high_w_width",core_high_w_width);
		parse_doub(label,value,"core_high_w_order",core_high_w_order);

		parse_doub(label,value,"max_w_beta",max_w_beta);
		parse_doub(label,value,"max_t_clnt",max_t_clnt);
		parse_doub(label,value,"max_w_xps",max_w_xps);

		parse_doub(label,value,"xps_width",xps_width);
		parse_doub(label,value,"xps_order",xps_order);
	}
	paramfile.close();
	cout<<"Done"<<endl;

	//We force these parameters to be positive.
	core_uniform_width = abs(core_uniform_width);
	core_uniform_order = abs(core_uniform_order);
	core_high_w_width = abs(core_high_w_width);
	core_high_w_order = abs(core_high_w_order);
	xps_width = abs(xps_width);
	xps_order = abs(xps_order);

	//Now we extract the core data from COREIN
	cout<<"Loading data..."<<endl;
	ifstream corefile(COREIN.c_str());
	string coreline;
	int corecounter=1;
	while(getline(corefile,coreline)){
		string token1, token2, token3, token4;	//This will parse for up to four tokens per line 
		istringstream iss(coreline);
		//We extract the number of time points from the fourth token of the third line 

		if(corecounter==3){
			iss>>token1>>token2>>token3>>token4;
			num_core_steps = int( atof( token4.c_str() ) ); //We extract and convert it to an integer
			//And we allocate the data arrays now that we know their size
			core_times = new double[num_core_steps];
			core = new double[num_core_steps];
			
		}

		//We now extract the time, data pairs, 
		//These are on the 5th lines to the end and are the first and second tokens respectively
		//We automatically convert to atomic units
		if(corecounter>4){
			iss>>token1>>token2;
			core_times[corecounter-5] = atof( token1.c_str() )/AUtime;	//We make sure to convert from fs to AU
			core[corecounter-5] = atof( token2.c_str() )/2.0;		//It is in Rydberg and we divide by 2 to get Hartrees
		}
		corecounter++;
	}
	corefile.close();

	cout<<"Done"<<endl;
	cout<<"Checking parameters and allocating..."<<endl;

	//We implement more value guards and defaults. Also compute deltas and cutoffs
	//From this we compuet the Nyquist frequency, the time spacing, ...
	delta_t_core = core_times[1]-core_times[0];
	max_t_core = core_times[num_core_steps-1];	//The last time point is the max time and the total elapsed time 
	cutoff_w_beta = pi/(delta_t_core*1.5);
	//Now we check to make sure that the other parameters are given in acceptable ranges
	//We also enforce defaults. 
	if(num_beta_steps<2){
		num_beta_steps = num_core_steps;	//The default value for the number of beta steps is the number of core steps 
	}
	if(num_clnt_steps<2){
		num_clnt_steps = num_core_steps;	//The default value for the number of cumulant steps is the number of core steps 
	}
	if(num_xps_steps<2){
		num_xps_steps = 2 * num_clnt_steps;	//The default for this is twice the number of beta steps 
	}
	
	//We make sure that the maximum frequency is not above the cutoff. If it is, we set it to the cutoff
	if(max_w_beta > cutoff_w_beta || max_w_beta <= 0.0){
		max_w_beta = cutoff_w_beta;
	}
	delta_w_beta = max_w_beta/double(num_beta_steps-1);
	//From this we compute the cutoff time
	cutoff_t_clnt = pi/(1.5*delta_w_beta);
	if(max_t_clnt > cutoff_t_clnt || max_t_clnt <= 0.0){
		max_t_clnt = cutoff_t_clnt;
	}
	delta_t_clnt = max_t_clnt/double(num_clnt_steps-1);
	//And now we compute the cutoff frequency for the XPS
	cutoff_w_xps = pi/(1.5*delta_t_clnt);
	if(max_w_xps > cutoff_w_xps || max_w_xps <= 0.0){
		max_w_xps = cutoff_w_xps;
	}
	delta_w_xps = 2.0*max_w_xps/double(num_xps_steps-1);	//We use a symmetric interval so delta must be twice as large

	//Now we can allocate all the remaining arrays 
	beta = new double[num_beta_steps];
	clnt = new complex<double>[num_clnt_steps];
	xps = new complex<double>[num_xps_steps];

	beta_freqs = new double[num_beta_steps];
	clnt_times = new double[num_clnt_steps];
	xps_freqs = new double[num_xps_steps];

	//We generate the frequencies and times 
	for(int i = 0; i < num_beta_steps; i++){
		beta_freqs[i] = double(i)*delta_w_beta;
	}
	for(int i = 0; i < num_clnt_steps; i++){
		clnt_times[i] = double(i)*delta_t_clnt;
	}
	for(int i = 0; i < num_xps_steps; i++){
		xps_freqs[i] = -max_w_xps + double(i)*delta_w_xps;
	}

	cout<<"Done"<<endl;
	cout<<"Processing core data..."<<endl;

	//Now we process the core data
	//Compute the average
	for(int i = 1; i < num_core_steps; i++){
		coreav1 += 0.5 * delta_t_core * (core[i] + core[i-1])/max_t_core;
	}

	//Subtract the average 
	for(int i = 0; i < num_core_steps; i++){
		core[i] -= coreav1;	
	}

	//Now we weight the core data and compute the average as we go
	for(int i = 0; i < num_core_steps; i++){
		core[i] *= exp(-pow(core_uniform_width*core_times[i],core_uniform_order));
		if(i>0){
			coreav2 += 0.5 * delta_t_core * (core[i] + core[i-1])/max_t_core;
		}
	}

	//Subtract average again
	for(int i = 0; i < num_core_steps; i++){
		core[i] -= coreav2;	
	}

	cout<<"Done"<<endl;
	cout<<"Computing beta..."<<endl;

	//We compute it using a double loop to do the Fourier transform using trapezoidal rule
	//beta(w) = w*re(integral_0^tmax dt e^{-iwt}*core(t))
	//We also include frequency dependent broadening of the form
	//e^(-(width*w*t)^order)
	for(int i = 0; i < num_beta_steps; i++){
		//initialize to zero
		beta[i] = 0.0;
		for(int j = 1; j < num_core_steps; j++){
			beta[i] += 0.5 * delta_t_core * beta_freqs[i] * real(
					core[j  ]*exp(-I*beta_freqs[i]*core_times[j  ] - pow(core_high_w_width*core_times[j  ]*beta_freqs[i],core_high_w_order)) +
					core[j-1]*exp(-I*beta_freqs[i]*core_times[j-1] - pow(core_high_w_width*core_times[j-1]*beta_freqs[i],core_high_w_order)) );
		}
		//Josh divides his beta by pi/2 so we do this as well
		//beta[i] /= pi/2.0;
	}

	cout<<"Done"<<endl;
	cout<<"Computing cumulant..."<<endl;

	//Now we must compute the cumulant
	//To do this we use the formula
	//C(t) = int_0^{wmax} dw beta(w)/w^2 (e^(iwt)-iwt -1)
	//Using trapezoidal rule
	//To avoid division by zero we using a small pole shift in the denominator 
	//We also compute the quasiparticle peak weight, the main peak weight, and the quasiparticle peak shift given by
	//quasi_peak_weight = int_0^wmax beta(w)/w^2
	//quasi_peak_shift = int_0^wmax beta(w)/w
	//main_peak_weight = e^(-quasi_peak_weight)
	complex<double> pole_shift = (1.0e-30)*I;	//We make it imaginary and small

	for(int i = 0; i < num_beta_steps; i++){
		sat_peak_weight +=  0.5 * delta_w_beta * real(
					beta[i  ]/pow(beta_freqs[i  ]-pole_shift,2)+
					beta[i-1]/pow(beta_freqs[i-1]-pole_shift,2));

		quasi_peak_shift +=  0.5 * delta_w_beta * real(
					beta[i  ]/pow(beta_freqs[i  ]-pole_shift,1)+
					beta[i-1]/pow(beta_freqs[i-1]-pole_shift,1));
	}

	quasi_peak_weight = exp(-sat_peak_weight);


	for(int i = 0; i < num_clnt_steps; i++){
		clnt[i] = 0.0;
		for(int j = 1; j< num_beta_steps; j++){
			clnt[i] += 0.5 * delta_w_beta * (
					beta[j  ]/pow(beta_freqs[j  ]-pole_shift,2)*(
					exp(I*beta_freqs[j  ]*clnt_times[i]) - I*beta_freqs[j  ]*clnt_times[i] - 1.0 )+
					beta[j-1]/pow(beta_freqs[j-1]-pole_shift,2)*(
					exp(I*beta_freqs[j-1]*clnt_times[i]) - I*beta_freqs[j-1]*clnt_times[i] - 1.0 ));
		}
	}

	cout<<"Done"<<endl;
	cout<<"Computing XPS..."<<endl;

	//Now we compute the XPS 
	//We have 
	//XPS(w) = I/pi int_{0}^tmax dt e^{-iwt +C(t) - (width*t)^2}
	//We use trapezoidal rule
	for(int i = 0; i < num_xps_steps; i++){
		xps[i] = 0.0;
		for(int j = 1; j < num_clnt_steps; j++){
			xps[i] += (I)/pi * 0.5 * delta_t_clnt * (
				exp(+I*xps_freqs[i]*clnt_times[j  ] - pow(xps_width*clnt_times[j  ],xps_order) + clnt[j  ])+
				exp(+I*xps_freqs[i]*clnt_times[j-1] - pow(xps_width*clnt_times[j-1],xps_order) + clnt[j-1]));
		}
	}

	//And we compute the normalization of the spectral function
	//We also use this loop to find the location and value of the global max 
	spec_peak_value = xps_freqs[0];

	for(int i =1; i < num_xps_steps; i++){
		spec_norm += 0.5*delta_w_xps * imag(xps[i]+xps[i-1]);
		if(imag(xps[i])>spec_peak_value){
			spec_peak_value = imag(xps[i]);	//The new largest value
			spec_peak_point = xps_freqs[i];	//And its location
		}
	}

	cout<<"Done"<<endl;
	cout<<"Computing Cumulative Spectral Function..."<<endl;
	//We compute the spectral function and its cumulative distribution and print them to a file 
	//We also flip the spectral function around and center the main peak
	//First we compute the spectral function
	//We do this by taking the imag(XPS) and getting each entry in reverse order, then putting that into the spectral function entry
	//We allocate the arrays
	spec_freqs = new double[num_xps_steps];
	spec_den = new double[num_xps_steps];
	spec_cum = new double[num_xps_steps];

	//Next we fill the spectral function array
	//We also fill the frequencies
	//we fill in reverse order and also flip the sign and subtract off spec_peak_point
	for(int i = 0; i < num_xps_steps; i++){
		spec_den[i] = imag(xps[num_xps_steps-1-i]);	//This flips the order
		spec_freqs[i] = -(xps_freqs[num_xps_steps-1-i] - spec_peak_point);
	}

	//Now we compute the cumulative distribution 
	//We use trapezoidal rule 

	for(int i = 0; i < num_xps_steps; i++){
		spec_cum[i] = 0.0;	//Zero out the value
		for(int j = 1; j < i; j++){
			//We only integrate up to i as this is the indefinite integral we are computing 
			spec_cum[i] += 0.5 * delta_w_xps * (spec_den[j-1]+spec_den[j]);
		}
	}
	cout<<"Done"<<endl;
	cout<<"Saving calculations..."<<endl;

	//We write the outputs to files
	ofstream coreout(COREOUT.c_str());
	ofstream betaout(BETAOUT.c_str());
	ofstream clntout(CLNTOUT.c_str());
	ofstream xpsout(XPSOUT.c_str());
	ofstream specout(SPECOUT.c_str());

	for(int i = 0; i < num_core_steps; i++){
		//We convert the time to fs and the core to Hartrees
		coreout<<AUtime*core_times[i]<<" "<<hartree*core[i]<<endl;
	}
	for(int i = 0; i < num_beta_steps; i++){
		//we convert the frequencies to eV and beta to eV
		betaout<<beta_freqs[i]*hartree<<" "<<hartree*beta[i]<<endl;
	}
	for(int i = 0; i < num_clnt_steps; i++){
		//Convert the time to fs but leave the cumulant unitless
		clntout<<AUtime*clnt_times[i]<<" "<<real(clnt[i])<<" "<<imag(clnt[i])<<endl;
	}
	for(int i = 0; i < num_xps_steps; i++){
		//We convert frequency to eV and the XPS we convert to fs
		xpsout<<hartree*xps_freqs[i]<<" "<<real(xps[i])/hartree<<" "<<imag(xps[i])/hartree<<endl;
	}
	for(int i = 0; i < num_xps_steps; i++){
		//We convert the frequency to eV and the spectral function to 1/eV. The cumulative is unitless
		specout<<hartree*spec_freqs[i]<<" "<<spec_den[i]/hartree<<" "<<spec_cum[i]<<endl;
	}

	coreout.close();
	betaout.close();
	clntout.close();
	xpsout.close();
	specout.close();

	cout<<"Done"<<endl;
	cout<<"Dumping parameters..."<<endl;
	
	//Now we compute the sumulative weight below a certain value
	//We want to know the integral of the XPS up to a given cutoff
	//If the cutoff is x and C is the inetgral of the spectral function,
	//Where the spectral function is the XPS flipped around the central peak,
	//Then we want int_{-infty}^x XPS(t)dt = 1-cumulative(-x)
	double weight_below_E = -5.0/hartree;	//We check the spectral weight up to this value (in eV). We convert to hartree by dividing by hartree
	int weight_below_E_index = snapToGrid(-weight_below_E,num_xps_steps,spec_freqs);	//This is the index of the weight for Cumlative(-weight_below_E) 
	double weight_below = 1.0 - spec_cum[weight_below_E_index]; //The cumulative weight of the XPS below the given cutoff. Note it is 1-cumulative because that is reveresed about the origin	

	//Now we write all the parameters to standard out 
	//The parameters are all written in atomic units 
	cout<<"--------------------------------------------------"<<endl;
	cout<<"INFILE "<<argv[1]<<endl;
	cout<<"COREIN "<<COREIN<<endl;
	cout<<"COREOUT "<<COREOUT<<endl;
	cout<<"BETAOUT "<<BETAOUT<<endl;
	cout<<"CLNTOUT "<<CLNTOUT<<endl;
	cout<<"XPSOUT "<<XPSOUT<<endl;
	cout<<"SPECOUT "<<SPECOUT<<endl;
	cout<<endl;
	cout<<"num_core_steps "<<num_core_steps<<endl;
	cout<<"num_beta_steps "<<num_beta_steps<<endl;
	cout<<"num_clnt_steps "<<num_clnt_steps<<endl;
	cout<<"num_xps_steps "<<num_xps_steps<<endl;
	cout<<endl;
	cout<<"delta_t_core "<<delta_t_core<<" AUtime"<<endl;
	cout<<"delta_w_beta "<<delta_w_beta<<" Hartree"<<endl;
	cout<<"delta_t_clnt "<<delta_t_clnt<<" AUtime"<<endl;
	cout<<"delta_w_xps "<<delta_w_xps<<" Hartree"<<endl;
	cout<<endl;
	cout<<"max_t_core "<<max_t_core<<" AUtime"<<endl;
	cout<<"max_w_beta "<<max_w_beta<<" Hartree"<<endl;
	cout<<"max_t_clnt "<<max_t_clnt<<" AUtime"<<endl;
	cout<<"max_w_xps "<<max_w_xps<<" Hartree"<<endl;
	cout<<endl;
	cout<<"cutoff_w_beta "<<cutoff_w_beta<<" Hartree"<<endl;
	cout<<"cutoff_t_clnt "<<cutoff_t_clnt<<" AUtime"<<endl;
	cout<<"cutoff_w_xps "<<cutoff_w_xps<<" Hartree"<<endl;
	cout<<endl;
	cout<<"core_uniform_width "<<core_uniform_width<<" Hartree"<<endl;
	cout<<"core_uniform_order "<<core_uniform_order<<endl;
	cout<<endl;
	cout<<"core_high_w_width "<<core_high_w_width<<endl;
	cout<<"core_high_w_order "<<core_high_w_order<<endl;
	cout<<endl;
	cout<<"xps_width "<<xps_width<<" Hartree"<<endl;
	cout<<"xps_order "<<xps_order<<endl;
	cout<<endl;
	cout<<"coreav1 "<<coreav1<<" Hartree"<<endl;
	cout<<"coreav2 "<<coreav2<<" Hartree"<<endl;
	cout<<"quasi_peak_shift "<<quasi_peak_shift*hartree<<" eV"<<endl;
	cout<<"sat_peak_weight "<<sat_peak_weight<<endl;
	cout<<"quasi_peak_weight "<<quasi_peak_weight<<endl;
	cout<<"spec_norm "<<spec_norm<<endl;
	cout<<"xps_peak_value "<<spec_peak_value/hartree<<" 1/eV"<<endl;
	cout<<"xps_peak_point "<<spec_peak_point*hartree<<" eV"<<endl;
	cout<<"XPS weight below "<<weight_below_E*hartree<<" eV: "<<weight_below<<endl;
	cout<<"--------------------------------------------------"<<endl;

	//Clean up allocations
	delete [] core_times;
	delete [] core;
	delete [] beta_freqs;
	delete [] beta;
	delete [] clnt;
	delete [] xps_freqs;
	delete [] xps;
	delete [] spec_freqs;
	delete [] spec_den;
	delete [] spec_cum;

	return 0;
}











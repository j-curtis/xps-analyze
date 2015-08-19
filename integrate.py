#This is a quick numpy script designed to compute the cumulative distribution of the spectral function
#Jonathan Curtis
#08/19/15

import numpy as np
import sys 

def main():
	#First we parse the command line arguments
	#We are expecting " program_name 'input_xps.dat' 'output.dat'"
	#We will load in the data from an xps data file which as the first column as frequency, the second as Real part, and the third column as the Imag part (spectral function).
	#We will print the file as column 1: frequency, column 2: spectral function, column 3: cumulative spectral function (integral)

	argc = len(sys.argv)
	
	if argc!=3:
		print "Usage error. Usage is: 'program name' 'input_xps.dat' 'output_spec.dat' "
		return 0
	
	INFILE = sys.argv[1]
	OUTFILE = sys.argv[2]

	#Next we load the data in. We use the numpy loadtxt option to get an array with the correct sizing
	
	in_data = np.loadtxt(INFILE,usecols=(0,2))	#We only need the  first and third columns
	
	print in_data.shape

if __name__ == "__main__":
	main()

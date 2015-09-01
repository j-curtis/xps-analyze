#This is a simple python script desigend to generate a test file to perform calculations with 
#Jonathan Curtis
#08/24/15

import numpy as np

#We choose a signal that is analytically solvable

def main():
	#PARAMETERS
	
	FILEOUT = "testcoreresponse.dat"

	num_t_steps = 10000

	freq = 5

	tot_t = 50

	#GENERATE THE SIGNAL
	
	t = np.linspace(0.0,tot_t,num_t_steps)

	f = np.sin(freq*t)

	#PRINT TO FILE, USING THE CORRECT FORMAT
	#The first 2 lines will just begin with #
	#The thrid line must read # token2 token3 num_t_steps
	#The fourth line must also be comments beginning with #
	
	with open(FILEOUT,"w") as fileout:
		fileout.write("#\n")
		fileout.write("#\n")
		fileout.write("# # # "+str(num_t_steps)+"\n")
		fileout.write("# \n")
	
		for i in xrange(num_t_steps):
			fileout.write(str(t[i])+" "+str(f[i])+"\n")

if __name__ == "__main__":
	main() 


# xps-analyze
C++ code that generates an XPS spectrum from raw core response data

To use the quick_analyze program, run it on a parameter file that is formatted as follows:

Each line will be parsed.
The first string of a line is treated as a variable flag. 
If the flag is recognized by the parser, it will try interpret the second string as the value for that parameter.
If it is not recognized, it will ignore that line.
If the parser fails to interpret the value for a flag, undefined behavior could result.
Any strings after the value on a given line will be ignored (and can be used for comments).
The lines can be parsed in any order.
Most variables have default values specified or enforced. However, it is encouraged that each variable be given a value.
The recognized flags are as follows:


COREIN = The file where the core response data will be read from.

COREOUT = The file where the processed core response data will be written to.

BETAOUT = The file where the function beta will be written to.

CLNTOUT = The file where the cumulant will be written to.

XPSOUT = The file where the XPS will be written to.

SPECOUT = The file where the spectral function density and cumulative distribution will be written to.

core_uniform_width = The width of the uniform spectral broadening applied to the core response data. Default of 0.0. Value is in Hartrees.

core_uniform_order = The order of the uniform spectral broadening applied to the core response data. Default of 2.0. Value is unitless.

core_high_w_width = The width of the high frequency broadening applied to the core response data. Default is 0.0. Value is unitless.

core_high_w_order = The order of the high frequency broadening applied to the core response data. Default is 2.0. Value is unitless.

xps_width = The width of the uniform spectral broadening applied to the XPS spectrum. Default is 0.0. Value is in Hartrees. 

xps_order = The order of the uniform spectral broadening applied to the XPS spectrum. Default is 2.0. Value is unitless.

num_beta_steps = The number of frequency steps used in computing beta. Default is num_core_steps.

num_clnt_steps = The number of time steps used in computing the cumulant. Default is num_core_steps.

num_xps_steps = The number of frequency steps used in calculation the XPS spectrum. Default is 2*num_clnt_steps.

max_w_beta = The highest frequency we compute beta to. Default is the Nyqusit frequency/1.5. Value is in Hartree.

max_t_clnt = The largest time we compute the cumulant to. Default is the Nyquist time/1.5 (analogous to the Nyquist frequency for an inverse FT). Value is in AUtime.

max_w_xps = The largest time we compute the xps to (symmetrically about 0). The default is the Nyquist frequency from the cumulant time step/1.5. Value is in Hartree.


The output of the program contains a detailed list of all the variables and results (except for the functions, which are printed to the given files). 
It is, in theory, possible to parse the output of the program into the program to run the same exact program, though this is not recommended or tested substantially.

The program add_background.cpp/bin accepts as an argument (in this order) an XPS file, then a SPEC file, and then the name of the output file. It uses the cumulative spectral function from SPEC to generate a background term and add it to the XPS from the XPS file. 
This is then printed to the output file with the columns: frequency, XPS+background, XPS, background.
It requires the XPS and SPEC files to use the same frequency grids. It also uses a preset background scaling constant which is a parameter in the program itself. This is called 'background_ratio' and is used to determine the relative weight of the background term (and has units of eV, so that the background term units work out).


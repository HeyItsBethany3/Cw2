This is a class system for solving the Elliptic PDE in Q1.

The 'Function1' and 'Function2' classes derive from 'AbstractFunction'.
They contain all the functions necessary to specify the elliptic problem; the
f(x) function and the exact u(x), so the problem setup can be easily changed.

The 'Elliptic' class is the main class for solving the PDE. It contains all the
relevant algorithms for solving the PDE and is constructed using an
'AbstractFunction' object.

'Driver.cpp' is the main cpp file that we are running our tests in.

Compile and run all the files by using the makefile.
Navigate to this directory using the terminal and then type 'make'.
Use 'make clean' to restore the directory back to its original state.

When creating plots, use the relevant function in driver. The code will
automatically copy the csv files created into a MATLAB directory but you have
to change the path in the code, since this will be different on each system.
Copy in the matlab files 'MatlabScripts' permanently into this your MATLAB
working directory. Now tests can be run quickly, by running 'Driver.cpp' and
immediately running the MATLAB files. The MATLAB files to run for each test
corresponds to the name of the csv file created.

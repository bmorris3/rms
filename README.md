# rms
Rotational Modulation Simulator

Requires STSP to be built on your machine, and that you create an environment
variable called `STSP_EXECUTABLE` which contains a path to the STS executable
file. 

### Installation

To build from source, first clone this git repository, `cd` into it, and run:
```
python setup.py install
```
When the `rms` installation is complete, cd into the `STSP` directory, and run
the following command to compile STSP: 
```
gcc -lm stsp.c -o stsp_rms
```
Lastly, create an environment variable in your path `.cshrc` or `.bashrc` file
called `STSP_EXECUTABLE` which contains the absolute path to the `stsp_rms` 
executable file you just compiled. Now you're ready to go!

### Examples

Check out the examples notebook in the `notebooks` directory.

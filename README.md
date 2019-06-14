# Compiling the Sources

There is a [makefile](src/makefile) file to compile the sources. This [makefile](src/makefile) mut be edited prior to the compilation to set the correct location of the CPLEX path and the name of the CPLEX libraries, which depend on the version you have installed. The current [makefile](src/makefile) is prepared for CPLEX 12.6.2 running on a Linux distribution.

Once you edited the [makefile](src/makefile), you can run it with the make command to generate an executable linked with the static library of CPLEX. This corresponds to the first goal of the [makefile](src/makefile), named NRP-static.

You also have the option to prepare an executable dynamically linked to CPLEX. To prepare this executable use the goal NRP-dynamic. In the latter case, remember to set/modify the LD_LIBRARY_PATH environment variable to contain the path to the dynamic library before running the executable. The output of the make command should report the exact command you need to do this.

# Running the Executable

The archive "parameters_NRP" must be in the same folder as the executable

```
./NRP (Input) (Max time execution) (Algorithm) (Option_1) [Option_2] [Option_3]
```

**Input**: Name of the archive to execute

**Max time execution**: Measured in seconds. When no limit time, type '0'

**Algorithm**: Choose between the following algorithms
 * 'econst1': Epsilon constraint with 1 ILP per iteration
 * 'econst2': Epsilon constraint with 2 ILPs per iteration
 * 'augmecon': Improved Epsilon constraint, AUGMEGON
 * 'spf': Supported Pareto Front
 * 'hybrid': Hybrid (parametrization + epsilon constraint)
 * 'tchebycheff': Augmented Thebycheff algorithm
 * 'mixed': Mixed two types of algorithms in {hybrid, tchebycheff}

# Algorithms and Options
	ALGORITHM  OPTION_1  [OPTION_2]  [OPTION_3]
		econst1 {f1,f2}
		econst2 {f1,f2}
		augmecon exact {f1,f2} {'fix','variable'}
		augmecon anytime {f1 , f2} {'fix','variable'}
		spf {exact,approximate}	
		hybrid exact
		hybrid anytime {normal,quadrants} {area,unexplored}
		tchebycheff {exact,anytime}	
		mixed {hybrid tchebycheff} {tchebycheff hybrid} [spf]


# Examples

* ./NRP nrp1 0 econst1 f1
* ./NRP nrp5 0 augmecon exact f1 variable
* ./NRP nrp-g2 0 augmecon exact f1 0.00001
* ./NRP nrp-e1 300 mixed hybrid tchebycheff
* ./NRP nrp-m4 60 mixed hybrid tchebycheff spf



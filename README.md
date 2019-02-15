# Short Description

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
	ALGORITHM  OPTION_1  OPTION_2  OPTION_3
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



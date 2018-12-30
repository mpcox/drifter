# drifter

This Perl program simulates genetic drift and sampling error in population genetic samples.  It is primarily intended to simulate drift in haploid populations of varying size given known, or predicted, starting conditions (i.e., allele frequencies and times of population origin).

The program was reported in:

Cox MP. 2006. [Extreme patterns of variance in small populations: Placing limits on human Y-chromosome diversity through time in the Vanuatu Archipelago](https://doi.org/10.1111/j.1469-1809.2006.00327.x). *Annals of Human Genetics* 71: 390-406.

INSTALLATION

*drifter* requires a standard working Perl installation and has been confirmed to work with Perl versions up to 5.18.2.

USAGE

For full information, see the complete [documentation file](documentation.pdf).

The following usage assumes a standard installation (*i.e.,* with *drifter.pl* aliased to *drifter*).

Basic command line information can be obtained by typing:

```
drifter
```

*drifter* offers a number of command line options:

```
-b         Verbose – outputs verbose comments
-h         Horizontal – outputs a horizontal format output file
-v         Vertical – outputs a vertical format output file
-t         Timefile – outputs a log file containing run time information
-s         Sampling – creates a sampling profile only for the input populations
-d NUMBER  Seed – allows the user to enter a random number seed
-g NUMBER  Generations – the number of generations for input populations to experience drift
-i NUMBER  Iterations – defines the number of independent replicates
```

*drifter* reads information from an input file.  The following example, ['solomons 1'](examples/solomons1.input), lists the populations to be simulated (here, Malaita and Rendova), their effective population sizes and the sample sizes of the current study, and the names and frequencies of the observed alleles (here, P, K, M and O).

```
#drifter
poplabel = [Malaita]
Np = [10000]
Ns = [12]
A1 = [K][0.667]
A2 = [S][0.083]
A3 = [M][0.250]
;
poplabel = [Rendova]
Np = [1000]
Ns = [20]
A1 = [K][0.550]
A2 = [O][0.450]
;
```

At its simplest, *drifer* is run by invoking the input file:

```
drifter input.file
```

EXAMPLE

For details, see the worked examples [here](examples).

To use *drifter* on the solomons 1 dataset for 10 generations and 2 iterations, run:

```
drifter -g 10 -i 2 solomons1.input
```

*drifter* generates random outcomes, but the output formatting looks like the following, with each block containing the final frequencies of the input alleles after 10 generations of drift:

```
Number of iterations: 2
Number of generations: 10

K    M    O    S
Malaita    0.58333    0.25000    0    0.16667    
Rendova    0.45000    0    0.55000    0    
;
Malaita    0.75000    0.08333    0    0.16667    
Rendova    0.65000    0    0.35000    0    
;
```


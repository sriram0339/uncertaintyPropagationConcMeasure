This project has the source code for the paper

Olivier Bouissou, Eric Goubault, Sylvie Putot, Aleksandar Chakarov, and Sriram Sankaranarayanan, Uncertainty Propagation using Probabilistic Affine Forms and Concentration of Measure Inequalities In Tools and Algorithms for Construction and Analysis of Systems (TACAS), Volume 9636 of Lecture Notes in Computer Science pp. 225-243 (2016). 

http://www.cs.colorado.edu/~srirams/papers/affine-forms-conc-measure-tacas-16.pdf

== Compilation Instructions

=== Prerequisites

You will need the following libraries installed:
  *  MPFR: Multiprecision floating point library
  * MPFI: An interval extension to MPFR
  * BOOST
  * BOOST TIMER

###  Compile

Simply go into the `src/` directory and run `Make`. You will need to edite `Makefile` to make sure that the compilation flags `-I` and `-L` point to the
right paths for the installed libraries. Currently, they are configured to a MAC OSX with homebrew setup.


This will produce many binaries called `ssExamples`, `eulerMaruyama`,  `filter` and `fersonExample`.


##  Running various benchmarks

Each of the benchmarks are described in the appendix of the paper referenced above.

| Benchmark Name |  Command | File where benchmark is coded |
| -------------  | ---------| -------------------------|
| Simple 2D arm motion example | `./ssExamples 1` | `ssExamples.cpp` |
| Cart steering example | `./ssExamples 2` |     |
| Tank filling | `./ssExamples 3 ` |     |
| Anesthesia model simulation | `./ssExamples 4` |   |
| Tumor Model | `./ssExamples 5` |   |
| Double Well Model | `./ssExamples 6` |   | 
| Rimless wheel model | `./ssExamples 7` |    |
| Cartpole balance model | `./ssExamples 8` |   | 
| Polynomial Cartpole balance model |  `./ssExamples 9 ` |   |
| Ferson et al Example | `./ferson` | ferson.cpp |
| Euler Maruyama SDE example | `./eulerMaruyama` | eulerMaruyma.cpp |
| Filter Example | `./filter` | filter.cpp |


## Questions/Comments

Sriram Sankaranarayanan <srirams@NAMEOFMYSTATE.EDU>

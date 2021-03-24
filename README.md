# Torquato-Jiao Sequential Linear Program Solver

## Preliminaries

In order to compile this program, you will need to have GSL and GLPK installed.
This version was tested with GSL version 1.15 and GLPK version 4.55.  You may 
also need to make sure that the #include statements in the .cpp and .h files are
correct for your GLPK installation.

## Usage

The program call for the TJ program is:

```bash
$ ./tj_3d iconfig.dat parameters.txt output
```

* Change "3d" to the dimension that you'd like to use.  While the program was
written with a lot of flexibility available at runtime, the vectors and
matrices used are hard-coded according to dimension due to their ubiquity,
so you need a different executable for each dimension.  Fortunately, most 
people will only need two different version at most!

* The first argument is the initial configuration that will be packed using TJ.
* The second argument tells TJ where to read parameter values such as the
influence sphere radius, compression and translation limits, etc.
* The third argument tells TJ where to write the final configuration as well as
any intermediate configurations (see the print_every parameter).  IMPORTANT: It 
will automatically append a ".dat" at the end of the name you supply.  This is
because intermediate configurations may be written as so many iterations; these
are named "output_xx.dat", where "xx" is the LP iteration after which the
packing was written.

## Packing file format

The TJ algorithm is designed to generate packings in any dimension.  The way 
that the files are written is as follows:

[dim]

[number of spheres]

[sphere diameter]

[Lattice (dim rows)]

[sphere centers (in global coordinates)]

Note that the lattice is written so that its ROWS span the fundamental cell.
The sphere centers' locations are written in GLOBAL coordinates, meaning that
these are their positions, INDEPENDENT OF THE FUNDAMENTAL CELL.  Their LOCAL
coordinates (i.e. in terms of the fundamental cell's lattice matrix) may be
obtained by multiplying the global coordinates by the lattice matrix's inverse,
and vice versa.

## Potential Issues

One potential issue that may arise while running TJ is that GLPK will print
messages saying:

    "Warning: numerical instability (primal simplex, phase I)"
    
or

    "Warning: numerical instability (primal simplex, phase II)"

The cause of this is normally that the primal basis has become infeasible due 
to roundoff errors.  GLPK will attempt to fix the problem and continue solving,
but in the event that this doesn't work, try lowering "feasible_tol" in the
parameters file.  The included parameters file has some more suggestions
regarding the proper selection of parameter values.

Another alternative is to implement a more advanced LP solver.  While GLPK is a
decent free option, better options exist that are free to those who are 
affiliated with an academic institution.  My personal favorite is Gurobi. The 
base class lp_class has been designed so that other solvers can inherit from it 
and be implemented seamlessly.  It should be straightforward to follow the 
example set by the GLPK implementation for your solver of choice. 

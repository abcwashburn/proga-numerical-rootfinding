# proga-numerical-rootfinding
This repository was made to showcase my work in Dr. Proga's physics lab at UNLV, where I worked during the Summer of 2015.
Contains programs that use numerical methods to find the roots of the Bernoulli function. These calculations aided the professor's research in black hole accretion.

Explanation of files:
* nr.h, nrutil_nr.h, nrtypes_nr.h
  * Header files used for the NR class.
  * Source: Numerical Recipes: The Art of Scientific Computing, 3rd edition CD.
* roots_bondi.cpp
  * Uses the bracketing and bisection numerical method to determine the roots of the Bernoulli function.  Output is printed to the terminal.
  * To run this program:
    * Ensure that nr.h, nrutil_nr.h, and nrtypes_nr.h are in the same folder as this file.
    * Change X1, X2, lambda, and position as needed.
    * Compile and execute.
* roots_bondi_table.cpp
  * Similar to roots_bondi.cpp, but output is printed to a text file.  This text file can then be used with a plotting program (such as      gnuplot) to plot the data onto a graph.
  * To run this program:
    * Ensure that nr.h, nrutil_nr.h, and nrtypes_nr.h are in the same folder as this file.
    * Change X1, X2, lambda, position, and filename as needed.
    * Compile and execute.
* roots_bondi_bulge_table.cpp
  * Similar to roots_bondi_table.cpp, but altered to simulate the accretion rate of a galaxy's bulge.  This text file can then be used       with a plotting program (such as gnuplot) to plot the data onto a graph.
  * To run this program:
    * Ensure that nr.h, nrutil_nr.h, and nrtypes_nr.h are in the same folder as this file.
    * Change X1, X2, lambda, position, filename, and the constraints in the k for loop as needed.
    * Compile and execute.
* Graphs (folder)
  * Some graph samples in postscript format. The data in these graphs are from the output of these programs.

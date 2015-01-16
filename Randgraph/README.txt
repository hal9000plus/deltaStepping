This directory contains a C++ command line tool for generating various families of random graphs (e.g, bullseye, hierarchical) in a simple text format.

Running the program with no arguments prints out a list of the options (graph type + edge length distribution).  The format for the graphs is: first line gives n, the number of vertices, numbered 0..n-1.  Subsequent lines give one edge each.  The undirected edge (u,v) with weight w appears as u v w<NEWLINE>   The last line is "-1".

Seth Pettie <pettie@mpi-sb.mpg.de> and Vijaya Ramachandran <vlr@cs.utexas.edu>.

February 6, 2006.


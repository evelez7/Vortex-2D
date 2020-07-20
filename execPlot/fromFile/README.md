This is a modified version of execPlot that reads a matrix from a file instead of doing the math.

The values of the matrix must be delimited by whitespaces. 
The matrix is expected to be NxN (not tested as NxM).
The values will be read as doubles.

The names of the .vtk files are derived from the name of the files being parsed.
Must use C++17 for std::filesystem.

Depending on the grid spacing used in the interpolation, `h_g` may need to be changed on line 71 of main.cpp
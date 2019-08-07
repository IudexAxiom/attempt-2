# attempt-2
My attempt-2.cc code, which is a modified step-26 from the tutorials.
This code is designed to run the Heat Equation, in 3D, using imported meshes, with Neumann Boundary Conditions(changeable).
Post processing performed on the solution include the gradient and the heat flux.

To start, first navigate yourself to deal.ii's "examples" folder, create a new folder in there called "attempt-2", and then place all three of my files into there(CMakeText, attempt-2, and TestCube).
You will have to make sure you change the name of the imported mesh files and change the location of the output.
I frequently use Gmsh for mesh creation..
http://gmsh.info/

The equation solved in attempt-2 is the standard heat equation, 
u_t - D*(u_xx + u_yy + u_zz) = f

Within the code, we have f = 0, and 3 boundary conditions applied to our domain.
Boundaries with an ID of: 0 are given 0 Dirichlet, 1 are given 0 Neumann, and 2 are given a constant 1 Neumann 
This is easily changed within the code.

The files are output in vtk format and easily read in VisIt, which is my main visulization tool.
Here is a link..
https://wci.llnl.gov/simulation/computer-codes/visit/features

Now, I know there are some warnings that pop up, and I will try to iron those things out eventually. 
Deepest apolgies for not having those iron out now.

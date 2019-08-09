// Gmsh project created on Fri Aug 09 10:32:55 2019
Point(1) = {0,0,0,.1};
Point(2) = {1,0,0,.1};
Point(3) = {1,1,0,.1};
Point(4) = {0,1,0,.1};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(1) = {1,2,3,4};
Transfinite Line{1:4} = 10;

Plane Surface(2) ={1};
Transfinite Surface{2} = {1,2,3,4};
Recombine Surface{2};

Extrude {0,0,1} {
Surface{2};
Layers{10};
Recombine;
}

Physical Surface(0) = {2};
Physical Surface(1) = {21,17,25,13};
Physical Surface(2) = {26};
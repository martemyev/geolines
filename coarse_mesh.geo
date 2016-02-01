LX = 500.0;
LY = 500.0;
p = LX/30;

// domain
Point(1)  = { 0,   0, 0, p};
Point(2)  = { 0,  LY, 0, p};
Point(3)  = { LX, LY, 0, p};
Point(4)  = { LX,  0, 0, p};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};

Transfinite Surface{6}; // Alternate; 

//mark
Physical Surface(1) = {6};
Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4};

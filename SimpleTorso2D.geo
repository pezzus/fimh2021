Hx = -4.0; Hy = 2.0;
R1 = 3.0;  R2 = 2.0;
Tx = 10.0; Ty = 15.0;

DefineConstant[sh = {0.04, Min 0, Max 1.0, Step 0.02,
                     Name "Heart mesh size"} ];
DefineConstant[st = {0.5, Min 0, Max 10.0, Step 1.0,
                     Name "Torso mesh size"} ];

sb = st;

SetFactory("OpenCASCADE");

// Split the chest with electrodes
Point(1) = {   0, 0, 0 };
Point(2) = { +Tx, 0, 0 };
Point(3) = { -Tx, 0, 0 };
Point(4) = { 0, +Ty, 0 };
Point(5) = { 0, -Ty, 0 };
Point(6) = { Tx*Cos(3.0/4*Pi), Ty*Sin(3.0/4*Pi), 0 };
Point(7) = { Tx*Cos(1.0/4*Pi), Ty*Sin(1.0/4*Pi), 0 };

// chest
Ellipse(1) = {5, 1, 5, 2};
Ellipse(2) = {2, 1, 4, 7};
Ellipse(3) = {7, 1, 4, 4};
Ellipse(4) = {4, 1, 4, 6};
Ellipse(5) = {6, 1, 4, 3};
Ellipse(6) = {3, 1, 5, 5};

// heart
Circle(11)  = {Hx, Hy, 0, R1, 0, 2*Pi};
Circle(12)  = {Hx, Hy, 0, R2, 0, 2*Pi};
Point(100)  = {Hx, Hy, 0.0, sb};

Curve Loop(1) = {1,2,3,4,5,6};  // chest
Curve Loop(2) = {11};           // epi
Curve Loop(3) = {12};           // endo

Plane Surface(1) = {1, 2};  // torso
Plane Surface(2) = {3, 2};  // heart
Plane Surface(3) = {3};     // blood

Characteristic Length{ PointsOf{ Curve{1,2,3,4,5}; } } = st;
Characteristic Length{ PointsOf{ Curve{11,12}; } } = sh;
Point{100} In Surface{3};

Mesh.Algorithm = 5;
Physical Surface("torso", 1) = {1};
Physical Surface("heart", 2) = {2};
Physical Surface("blood", 3) = {3};
Physical Curve("chest", 100) = {1,2,3,4,5,6};
Physical Curve("epi",   110) = {11};
Physical Curve("endo",  120) = {12};
//Physical Point("centre", 1000) = {100};

Physical Point("elecVF", 1001) = {5};
Physical Point("elecVR", 1002) = {7};
Physical Point("elecVL", 1003) = {6};
Physical Point("elecV1", 1004) = {3};
Physical Point("s0", 1005) = {2};


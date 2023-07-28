SetFactory("OpenCASCADE");
Merge "ShroomLock-V7-Lock.step";
//+
//+
Physical Surface("screws", 1) = {27, 25, 31, 29};
//+
//+
//+
Physical Surface("inner-radius", 2) = {33, 32, 39, 38, 37, 36, 35, 34};
//+
Physical Volume("volume", 3) = {1};
//+
Field[1] = Cylinder;
//+
Field[1].Radius = 9;
//+
Field[1].VIn = 0.3; // 0.4; // 0.6;
//+
Field[1].VOut = 3;
//+
Field[1].ZAxis = 3;
//+
Background Field = 1;
//
// export in meters
Mesh.ScalingFactor = 0.001;//+


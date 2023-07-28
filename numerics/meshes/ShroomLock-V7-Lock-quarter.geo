SetFactory("OpenCASCADE");
Merge "ShroomLock-V7-Lock-quarter.step";
//+
//+
Physical Surface("inner-radius", 2) = {17, 16};
//+
Physical Surface("screws", 1) = {15};
//+
Physical Surface("sym-x", 101) = {31};
//+
Physical Surface("sym-y", 102) = {32};
//+
Physical Volume("volume", 3) = {1};
//+
Field[1] = Cylinder;
//+
Field[1].Radius = 8;
//+
Field[1].VIn = 0.3;
//+
Field[1].VOut = 3;
//+
Field[1].ZAxis = 3;
//+
Background Field = 1;
//
// export in meters
Mesh.ScalingFactor = 0.001;//+

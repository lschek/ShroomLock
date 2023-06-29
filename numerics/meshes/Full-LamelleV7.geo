SetFactory("OpenCASCADE");
Merge "Full_LamelleV7.step";
//+
// specify a global mesh size
//Mesh.MeshSizeMin = 0;
//Mesh.MeshSizeMax = 0.7;
//+
Physical Surface("screws", 1) = {43, 52, 49, 46};
//+
Physical Surface("inner-radius", 2) = {21, 17, 13, 9, 37, 33, 25, 29};
//+
Physical Volume("volume", 3) = {1};
//+
Field[1] = Cylinder;
//+
Field[1].Radius = 10;
//+
Field[1].VIn = 0.6; // 0.3;
//+
Field[1].VOut = 4;
//+
Field[1].ZAxis = 1.5;
//+
Background Field = 1;
//
// export in meters
Mesh.ScalingFactor = 0.001;
//+



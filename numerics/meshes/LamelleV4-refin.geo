SetFactory("OpenCASCADE");
Merge "LamelleV4.step";
//+
// specify a global mesh size
//Mesh.MeshSizeMin = 0;
//Mesh.MeshSizeMax = 0.7;
//+
Physical Surface("screws", 1) = {50, 47, 44, 41};
//+
Physical Surface("inner-radius", 2) = {35, 31, 7, 11, 23, 27, 19, 15};
//+
Physical Volume("volume", 3) = {1};
// specify a global mesh size
//Mesh.MeshSizeMin = 0;
//Mesh.MeshSizeMax = 0.7;
//+
//+
Field[1] = Cylinder;
//+
Field[1].Radius = 10;
//+
Field[1].VIn = 0.3; // 0.3;
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


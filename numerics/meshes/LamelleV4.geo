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
Mesh.MeshSizeMin = 0;
Mesh.MeshSizeMax = 0.7;

// export in meters
Mesh.ScalingFactor = 0.001;
//+


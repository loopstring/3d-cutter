# 3d-cutter
Simplified implementation of the paper "An Adaptive Virtual Node Algorithm with Robust Mesh Cutting" http://www.seas.upenn.edu/~cffjiang/research/cut/paper.pdf. For now it doesn't support robust intersection and cut through node or edge or face.

Written in C++11, the header only cutter only depend on standard library with a simply interface:
TetMesh Cutter3D::run(const TetMesh& tetMesh, const TriMesh& triMesh).

Easy construction of the input meshes:
TetMesh<T>(vector<array<T,3>>& nodes, vector<array<int,4>>& mesh)
TriMesh<T>(vector<array<T,3>>& nodes, vector<array<int,3>>& mesh)

main.cpp is an interactive cutting interface that depends on OpenGL and GLUT.

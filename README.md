# 3d-cutter

## Description
If you want to cut a tetrahedral mesh with an arbitrary 3d surface, this is the code to use!

It's a simplified implementation of the paper "An Adaptive Virtual Node Algorithm with Robust Mesh Cutting" http://www.seas.upenn.edu/~cffjiang/research/cut/paper.pdf. The resolution of the cutting details is constrained by the tetrahedral mesh, but as the mesh is adaptively refined based on the cuts, so you can sequentially apply arbitrary number of cuts. For now it doesn't support robust intersection and cut through node or edge or face, so you will need to perturb your meshes to avoid these degenerate cases if they don't split.

## Usage
Written in C++11, the header only cutter only depend on standard library, and it has a simple interface:

TetMesh Cutter3D::run(const TetMesh& tetMesh, const TriMesh& triMesh).

You simply feed the cutter a tetrahedral mesh and a cutting surface mesh, and it will return a new tetrahedral mesh that is cut.

The meshes can be easily constructed from std::vector and std::array, so just convert your meshes into these std data structures and start running the cutter!

TetMesh<T>(vector<array<T,3>>& nodes, vector<array<int,4>>& mesh)

TriMesh<T>(vector<array<T,3>>& nodes, vector<array<int,3>>& mesh)


## Demo
main.cpp is an interactive cutting interface that depends on OpenGL and GLUT. It supports drawing a curve that will be extended to a cutting surface, dragging pieces around, rotation and shifting, for more details, take a look at the comments in the keyboard callback function "void key(unsigned char key, int x, int y)".

The video demo.mov is added to showcase the cutting effects.

# 3d-cutter

## Description
If you want to cut a tetrahedral mesh with an arbitrary surface in 3d, this is the code to use!

It's a simplified implementation of the paper "An Adaptive Virtual Node Algorithm with Robust Mesh Cutting" http://www.seas.upenn.edu/~cffjiang/research/cut/paper.pdf. The resolution of the cutting details is constrained by the tetrahedral mesh, which is adaptively refined after each cut. For now it doesn't support robust intersection computations or cutting through nodes, edges or faces, so you will need to perturb your meshes to avoid these degenerate cases if they happen.

## Usage
Written in C++11, the header only cutter depends on the standard library alone, and it has a simple interface:

TetMesh Cutter3D::run(const TetMesh& tetMesh, const TriMesh& triMesh).

You simply feed the cutter a tetrahedral mesh and a cutting surface mesh, and it will return a new tetrahedral mesh that is cut.

The meshes can be easily constructed from std::vector and std::array, so just convert your meshes into these std data structures and start running the cutter!

TetMesh<T>(vector<array<T,3>>& nodes, vector<array<int,4>>& mesh)

TriMesh<T>(vector<array<T,3>>& nodes, vector<array<int,3>>& mesh)


## Demo
main.cpp is an interactive cutting interface that depends on OpenGL and GLUT. It supports drawing a curve that will be extended along the z direction into a cutting surface, dragging pieces around, rotation and shifting. For more details, please see the comments in the keyboard callback function "void key(unsigned char key, int x, int y)".

The video demo.mov is added to showcase the cutting effects.

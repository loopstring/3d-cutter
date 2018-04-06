//
//  Copyright (c) 2018 Yuting Wang. All rights reserved.
//

#ifdef __APPLE__//Mac OS
//#  include <GL/glew.h>
#  include <OpenGL/OpenGL.h>
#  include <GLUT/glut.h>
#else
//#  include <GL/glew.h>
#  include <GL/freeglut.h>
#  include <GL/freeglut_ext.h>
#endif  // __APPLE__

#include <iostream>
#include <vector>
#include <array>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include "Cutting.h"
#include <algorithm>
#include <map>
#include <fstream>
#include <sstream>
using namespace std;

#define d 3
#define T float
#define TV array<T, d>
#define T2 array<T, 2>
#define T4 array<T, 4>
#define Idp1 array<int, d+1>
#define Id array<int, d>
#define I3 array<int, 3>
#define I2 array<int, 2>
#define I4 array<int, 4>

#ifdef far
#undef far
#endif

#ifdef near
#undef near
#endif

TetMesh<T> tetMesh;
TriMesh<T> triMesh;

void zoom(T scaleSpeed) {
    for (auto& v : tetMesh.nodes_) {
        v = divide<T, d>(v, 1/scaleSpeed);
    }
    for (auto& v : triMesh.nodes_) {
        v = divide<T, d>(v, 1/scaleSpeed);
    }
}

void shift(TV shiftSpeed) {
    for (auto& v : tetMesh.nodes_) {
        v = add<T, d>(v, shiftSpeed);
    }
    for (auto& v : triMesh.nodes_) {
        v = add<T, d>(v, shiftSpeed);
    }
}

void simpleTet1() {
    tetMesh.nodes_.clear();
    tetMesh.mesh_.clear();
    tetMesh.nodes_.push_back(TV{0,0,0});
    tetMesh.nodes_.push_back(TV{1,0,0});
    tetMesh.nodes_.push_back(TV{0,1,0});
    tetMesh.nodes_.push_back(TV{0,0,1});
    tetMesh.nodes_.push_back(TV{1,1,1});
    tetMesh.mesh_.push_back(Idp1{0,1,2,3});
    tetMesh.mesh_.push_back(Idp1{4,2,1,3});
//    array<T,4> w;
//    Cutter3D<T>::computeIntersection(array<TV,4>{tetMesh.nodes_[0],tetMesh.nodes_[1],tetMesh.nodes_[2],tetMesh.nodes_[3]}, array<TV,1>{TV{0.1,0.1,0.1}}, w);
//    cout << w[0] << ", " << w[1] << ", " << w[2] << ", " << w[3] << endl;
    tetMesh.initializeSurfaceMesh();
    tetMesh.computeConnectedComponents();
    
    triMesh.nodes_.push_back(TV{0.3f,-0.1f,1.1f});
    triMesh.nodes_.push_back(TV{0.3f,-0.1f,-1});
    triMesh.nodes_.push_back(TV{0.3f,1.1f,1.1f});
    triMesh.nodes_.push_back(TV{0.3f,1.1f,-1});
    triMesh.mesh_.push_back(I3{0,1,2});
    triMesh.mesh_.push_back(I3{2,3,1});
    
    Cutter3D<T> cutter;
    tetMesh = cutter.run(tetMesh, triMesh);
}

void loadTriMesh(const string& filename) {
    ifstream fs(filename);
    triMesh.nodes_.clear();
    triMesh.mesh_.clear();
    string line;
    int l = 0;
    while(getline(fs, line)) {
        if (line.find("end_header") != std::string::npos) {
            while(getline(fs, line)) {
                if (line.substr(0,2) != "3 ") {
                    stringstream ss(line);
                    T t1,t2,t3;
                    ss >> t1 >> t2 >> t3;
                    triMesh.nodes_.push_back(divide<T,3>(TV{t1,t2,t3},13));
                } else {
                    stringstream ss(line);
                    int t1,t2,t3;
                    ss >> t1 >> t1 >> t2 >> t3;
                    triMesh.mesh_.push_back(I3{t1,t2,t3});
                }
            }
            break;
        }
    }
    cout << triMesh.nodes_.size() << " " << triMesh.mesh_.size() << endl;
}

void meshCutGrid() {
    tetMesh.nodes_.clear();
    tetMesh.mesh_.clear();
    int width = 23;
    int height = 23;
    int depth = 23;
    T low = -0.5;
    T high = 0.5;
    T left = -0.5;
    T right = 0.5;
    T far = 0.6;
    T near = -0.5;
    T lw = (right-left)/width;
    T lh = (high-low)/height;
    T ld = (far-near)/depth;
    int offset = (width+1)*(depth+1);

    for (int i = 0; i < height+1; i++) {
        //cout << "height: " << low+i*lh << endl;
        for (int j = 0; j < width+1; j++) {
            for (int k = 0; k < depth+1; k++) {
                tetMesh.nodes_.push_back(TV{left+j*lw, low+i*lh,near+k*ld});
            }
        }
    }
    
    int tetIndex = 0;
    for (int k = 0; k < height; k++) {
        for (int i = 0; i < width; i++) {
            for (int j = 0; j < depth; j++) {
                tetMesh.mesh_.push_back(I4{tetIndex+offset+1, tetIndex+1, tetIndex, tetIndex+depth+1});
                tetMesh.mesh_.push_back(I4{tetIndex, tetIndex+offset, tetIndex+offset+1, tetIndex+depth+1});
                tetMesh.mesh_.push_back(I4{tetIndex+depth+1, tetIndex+offset, tetIndex+offset+1, tetIndex+offset+depth+1});
                tetMesh.mesh_.push_back(I4{tetIndex+1, tetIndex+depth+1, tetIndex+offset+1, tetIndex+depth+2});
                tetMesh.mesh_.push_back(I4{tetIndex+depth+1, tetIndex+depth+2, tetIndex+offset+depth+2, tetIndex+offset+1});
                tetMesh.mesh_.push_back(I4{tetIndex+depth+1, tetIndex+depth+2+offset, tetIndex+offset+depth+1, tetIndex+offset+1});
                tetIndex++;
            }
            tetIndex++;
        }
        tetIndex += (depth+1);
    }

    tetMesh.initializeSurfaceMesh();
    tetMesh.computeConnectedComponents();
    
//    triMesh.nodes_.push_back(TV{0.3,-0.1,1.1});
//    triMesh.nodes_.push_back(TV{0.3,-0.1,-1});
//    triMesh.nodes_.push_back(TV{0.3,1.1,1.1});
//    triMesh.nodes_.push_back(TV{0.3,1.1,-1});
//    triMesh.mesh_.push_back(I3{0,1,2});
//    triMesh.mesh_.push_back(I3{2,3,1});

//    loadTriMesh("/Users/yutingwang/Downloads/beethoven.ply");
//    Cutter3D<T> cutter;
//    tetMesh = cutter.run(tetMesh, triMesh);
}


T scaleSpeed=1.1;
T shiftSpeed=0.02;
bool cutting = false;
void key(unsigned char key, int x, int y) {
    switch( key ) {
        case 033: // Escape Key
            exit(EXIT_SUCCESS);
            break;
        case 'c': // enter/exit cutting mode. In cutting mode, drag the mouse to draw a cutting curve. Otherwise, drag the red tet mesh to move it, or drag the mouse in the background to rotate the scene.
            cutting = !cutting;
            break;
        case 'r': // remove cutting triangle mesh
            triMesh.nodes_.clear();
            triMesh.mesh_.clear();
            break;
        case 'i': // zoom in
            zoom(scaleSpeed);
            break;
        case 'o': // zoom out
            zoom(1/scaleSpeed);
            break;
        case 'a': // shift left
            shift(TV{-shiftSpeed,0,0});
            break;
        case 'd': // shift right
            shift(TV{shiftSpeed,0,0});
            break;
        case 'w': // shift up
            shift(TV{0,shiftSpeed,0});
            break;
        case 's': // shift down
            shift(TV{0,-shiftSpeed,0});
            break;
        case 'f': // reload initial meshes
            meshCutGrid();
            break;
    }
    glutPostRedisplay();
}

int windowWidth = 600;
int windowHeight = 600;
T2 startingPosition;
T transSpeed=0.1;
set<int> draggingParticles;
int triMeshRes = 30;
T triMeshNear = -1.1;
T triMeshFar = 1.1;
T triMeshEdgeLen = (triMeshFar-triMeshNear) / triMeshRes;
void mouse(int button, int state, int x, int y) {
    array<T, 2> location{2*x/T(windowWidth)-1,1-2*y/T(windowHeight)};
    if(button==4){
        zoom(scaleSpeed);
        glutPostRedisplay();
        return;
    } else if(button==3){
        zoom(1/scaleSpeed);
        glutPostRedisplay();
        return;
    }
    if(state==GLUT_DOWN){
        if(button==GLUT_LEFT_BUTTON){
            startingPosition = location;
            if (cutting) {
                triMesh.clear();
                for (int i = 0; i <= triMeshRes; ++i) {
                    triMesh.nodes_.push_back(TV{location[0], location[1], (T)-1.1+i*triMeshEdgeLen});
                }
            } else {
                glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
                glBegin(GL_TRIANGLES);
                for(int i = 0; i < tetMesh.mesh_.size(); ++i) {
                    glColor4d(tetMesh.connectedComponents_[i]/255.,0,0,0);
                    const auto& tet = tetMesh.mesh_[i];
                    for (const auto& fi: FaceIndexes) {
                        for (int k = 0; k < 3; ++k) {
                            auto p = tetMesh.nodes_[tet[fi[k]]];
                            glVertex3d(p[0], p[1], p[2]);
                        }
                    }
                }
                glEnd();
                unsigned char val;
                glReadPixels(x, windowHeight - y, 1, 1, GL_RED, GL_UNSIGNED_BYTE, &val);
                draggingParticles.clear();
                for (int i = 0; i < tetMesh.connectedComponents_.size(); ++i) {
                    if (tetMesh.connectedComponents_[i] == (int)val) {
                        for (auto j : tetMesh.mesh_[i]) {
                            draggingParticles.insert(j);
                        }
                    }
                }
            }
        }
    }
    else if(state==GLUT_UP){
        if (cutting) {
            if (triMesh.mesh_.size()==0) {
                return;
            }
            tetMesh = Cutter3D<T>::run(tetMesh, triMesh);
            glutPostRedisplay();
        }
    }
}

void motion(int x, int y) {
    T2 location{2*x/T(windowWidth)-1,1-2*y/T(windowHeight)};
    auto shift = substract<T, 2>(location, startingPosition);
    T l = sqrt(pow(shift[0],2)+pow(shift[1],2));
    if (fabs(shift[0]) < 0.002 || fabs(shift[1]) < 0.002) {
        return;
    }
    if (cutting){
        //cout << location[0] << ", " << location[1] << endl;
        for (int i = 0; i <= triMeshRes; ++i) {
            if (i > 0) {
                int n = triMesh.nodes_.size();
                triMesh.mesh_.push_back(I3{n, n-1, n-triMeshRes-1});
                triMesh.mesh_.push_back(I3{n-triMeshRes-1, n-1, n-triMeshRes-2});
            }
            triMesh.nodes_.push_back(TV{location[0], location[1], (T)-1.1+i*triMeshEdgeLen});
        }
    } else {
        if (draggingParticles.size()) {
            for (auto i : draggingParticles) {
                tetMesh.nodes_[i] = add<T, d>(tetMesh.nodes_[i], TV{shift[0], shift[1], 0});
            }
        } else {
            T dx = -shift[1];
            T dy = shift[0];
            T costheta = cos(-l);
            T sintheta = sin(-l);
            dx /= l;
            dy /= l;
            array<array<T,d>,d> rotationMatrix{
                TV{costheta+dx*dx*(1-costheta), dx*dy*(1-costheta), dy*sintheta},
                TV{dx*dy*(1-costheta), dy*dy*(1-costheta)+costheta, -dx*sintheta},
                TV{-dy*sintheta, dx*sintheta, costheta}
            };
            for (auto& p: tetMesh.nodes_) {
                p = TV{dot<T,d>(rotationMatrix[0], p),dot<T,d>(rotationMatrix[1], p),dot<T,d>(rotationMatrix[2], p)};
            }
            //print<T,d>(tetMesh.nodes_);
            for (auto& p: triMesh.nodes_) {
                p = TV{dot<T,d>(rotationMatrix[0], p),dot<T,d>(rotationMatrix[1], p),dot<T,d>(rotationMatrix[2], p)};
            }
        }
        startingPosition = location;
    }
    glutPostRedisplay();
}

void render(const vector<TV>& nodes, const vector<I3>& faces, const T4& triColor, const T4& segColor) {
    vector<TV> tris, segs;
    for (const auto& face: faces) {
        tris.push_back(nodes[face[0]]);
        tris.push_back(nodes[face[1]]);
        tris.push_back(nodes[face[2]]);
        for (int i = 0; i < 3; ++i) {
            segs.push_back(nodes[face[i]]);
            segs.push_back(nodes[face[(i+1)%3]]);
        }
    }
    glVertexPointer(3, GL_FLOAT, 0, segs.data());
    glColor4f(segColor[0], segColor[1], segColor[2], segColor[3]);
    glDrawArrays(GL_LINES, 0, segs.size());
    
    glVertexPointer(3, GL_FLOAT, 0, tris.data());
    glColor4f(triColor[0], triColor[1], triColor[2], triColor[3]);
    glDrawArrays(GL_TRIANGLES, 0, tris.size());
}

void render() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnableClientState(GL_VERTEX_ARRAY);
    
    render(tetMesh.nodes_, tetMesh.surfaceMesh_, T4{1,0,0,1}, T4{0,1,1,1});
    render(triMesh.nodes_, triMesh.mesh_, T4{0,0,1,1}, T4{0,1,1,1});
    
    glutSwapBuffers();
    glDisableClientState(GL_VERTEX_ARRAY);
}

void reshape(GLint newWidth,GLint newHeight) {
    glViewport(0, 0, newWidth, newHeight);
    windowWidth = newWidth;
    windowHeight = newHeight;
    glutPostRedisplay();
}

int main(int argc, char **argv) {
    meshCutGrid();
    
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA | GLUT_ALPHA);
    glutInitWindowSize(windowWidth, windowHeight);
    string window_name="cutting";
    glutCreateWindow(window_name.c_str());
    glClearColor(0.0,0.0,0.0,1.0);
    glEnable(GL_DEPTH_TEST);
    
    //glutSpecialFunc(Special_Key);
    glutKeyboardFunc(key);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    glutDisplayFunc(render);
    glutReshapeFunc(reshape);
    
    glutMainLoop();
    
    return 0;
}

//
//  Copyright (c) 2018 Yuting Wang. All rights reserved.
//

#ifndef DataStructures_h
#define DataStructures_h

#include <vector>
#include <array>
#include <limits>
#include <algorithm>
#include <iostream>
#ifdef WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif
#include <stack>
#include <map>
#include <numeric>

using namespace std;

template<typename T, int d>
void print(const array<T, d>& v) {
    cout << "{";
    for (int i = 0; i < d; ++i) {
        cout << v[i] << ", ";
    }
    cout << "}\n";
}

template<typename T, int d>
void print(const vector<array<T, d>>& v) {
    cout << "{\n";
    for (const auto& a: v) {
        print<T,d>(a);
    }
    cout << "}\n";
}

template<typename T, int d1, int d2>
void print(const array<array<T, d1>, d2>& v) {
    cout << "{\n";
    for (const auto& a: v) {
        print<T,d1>(a);
    }
    cout << "}\n";
}

template<typename T>
T sorted(const T& v) {
    auto sv = v;
    sort(sv.begin(), sv.end());
    return sv;
}

class UnionFind {
    int *id, cnt, *sz;
public:
    // Create an empty union find data structure with N isolated sets.
    UnionFind(int N)   {
        cnt = N;
        id = new int[N];
        sz = new int[N];
        for(int i=0; i<N; i++)	{
            id[i] = i;
            sz[i] = 1;
        }
    }
    ~UnionFind()	{
        delete [] id;
        delete [] sz;
    }
    // Return the id of component corresponding to object p.
    int find(int p)	{
        int root = p;
        while (root != id[root])
            root = id[root];
        while (p != root) {
            int newp = id[p];
            id[p] = root;
            p = newp;
        }
        return root;
    }
    // Replace sets containing x and y with their union.
    void merge(int x, int y)	{
        int i = find(x);
        int j = find(y);
        if (i == j) return;
        
        // make smaller root point to larger one
        if   (sz[i] < sz[j])	{
            id[i] = j;
            sz[j] += sz[i];
        } else	{
            id[j] = i;
            sz[i] += sz[j];
        }
        cnt--;
    }
    // Are objects x and y in the same set?
    bool connected(int x, int y)    {
        return find(x) == find(y);
    }
    // Return the number of disjoint sets.
    int count() {
        return cnt;
    }
};

template<typename T, int d>
array<T,4> toI4(const array<T,d>& Id, T fill = -1) {
    array<T,4> a;
    a.fill(fill);
    for (int i = 0; i < min(d, 4); ++i) {
        a[i] = Id[i];
    }
    return a;
}
template<typename T, int d>
array<T, d> add(const array<T, d>& a1, const array<T, d>& a2) {
    array<T, d> s;
    for (int j = 0; j < d; ++j) {
        s[j] = a1[j] + a2[j];
    }
    return s;
}

template<typename T, int d>
array<T, d> substract(const array<T, d>& a1, const array<T, d>& a2) {
    array<T, d> s;
    for (int j = 0; j < d; ++j) {
        s[j] = a1[j] - a2[j];
    }
    return s;
}

template<typename T, int d>
array<T, d> divide(const array<T, d>& a, T t) {
    array<T, d> s;
    for (int j = 0; j < d; ++j) {
        s[j] = a[j] / t;
    }
    return s;
}

template<typename T, int d, int n>
array<T, d> center(const array<array<T, d>, n>& nodes, const array<T, n> weights) {
    array<T, d> c;
    for (size_t i = 0; i < d; ++i) {
        c[i] = 0;
        for (size_t j = 0; j < n; ++j) {
            c[i] += nodes[j][i] * weights[j];
        }
    }
    return c;
}

template<typename T, int d>
array<T, d> elementCenter(const vector<array<T,d>>& nodes, const array<int, 4>& element, const array<T, 4>& weights) {
    array<T, d> c;
    c.fill(0);
    for (size_t i = 0; i < 4; ++i) {
        if (element[i] >= 0) {
            const auto& node = nodes[element[i]];
            auto w = weights[i];
            for (int j = 0; j < d; ++j) {
                c[j] += (node[j] * w);
            }
        } else {
            break;
        }
    }
    return c;
}

template<typename T, int d>
T dot(const array<T, d>& a1, const array<T, d>& a2) {
    T s = 0;
    for (int j = 0; j < d; ++j) {
        s += a1[j] * a2[j];
    }
    return s;
}

template<typename T>
T cross(const array<T, 2>& a1, const array<T, 2>& a2) {
    return a1[0] * a2[1] - a1[1] * a2[0];
}

template<typename T>
array<T, 3> cross(const array<T, 3>& a1, const array<T, 3>& a2) {
    array<T, 3> s;
    for (int j = 0; j < 3; ++j) {
        s[j] = a1[(j+1)%3] * a2[(j+2)%3] - a2[(j+1)%3] * a1[(j+2)%3];
    }
    return s;
}

template<typename T>
T volume(const array<T,3>& node1, const array<T,3>& node2, const array<T,3>& node3, const array<T,3>& node4) {
    return dot<T,3>(cross<T>(substract<T,3>(node2,node1), substract<T,3>(node3,node1)), substract<T,3>(node4,node1));
}

template<typename T, int d>
T norm(const array<T, d>& a1) {
    return sqrt(dot(a1, a1));
}

template<typename T, int d>
T pointEdgeWeight(const array<T, d>& e1, const array<T, d>& e2, const array<T, d>& p) {
    array<T, d> v1 = substract<T, d>(p, e1);
    array<T, d> v2 = substract<T, d>(e2, e1);
    return dot<T, d>(v1, v2) / dot<T, d>(v2, v2);
}

template<typename T, int d>
bool pointOnEdge(const array<T, d>& e1, const array<T, d>& e2, const array<T, d>& p, T& w) {
    array<T, d> v1 = substract<T, d>(p, e1);
    array<T, d> v2 = substract<T, d>(e1, e2);
    if (cross<T>(v1, v2) == 0) {
        w = pointEdgeWeight<T, d>(e1, e2, p);
        return true;
    }
    return false;
}

template<typename T>
bool edgeEdgeIntersect(const array<T, 2>& p1, const array<T, 2>& p2, const array<T, 2>& q1, const array<T, 2>& q2, T& w) {
    array<T, 2> e1 = substract<T, 2>(p2, p1);
    T a1 = cross<T>(substract<T, 2>(q1, p1), e1);
    T a2 = cross<T>(substract<T, 2>(q2, p1), e1);
    if ((a1 < 0 && a2 > 0) || (a1 > 0 && a2 < 0)) {
        array<T, 2> e2 = substract<T, 2>(q2, q1);
        T a3 = cross<T>(substract<T, 2>(p1, q1), e2);
        T a4 = cross<T>(substract<T, 2>(p2, q1), e2);
        if ((a3 < 0 && a4 > 0) || (a3 > 0 && a4 < 0)) {
            w = a3 / (a3 - a4);
            return true;
        }
    }
    return false;
}

template<typename T>
bool pointInTriangle(const array<array<T, 2>, 3>& triangle, const array<T, 2>& point, array<T, 3>& w) {
    array<T, 3> areas;
    for (int i = 0; i < 3; ++i) {
        areas[(i+2)%3] = cross<T>(substract<T, 2>(point, triangle[i]), substract<T, 2>(triangle[(i+1)%3], triangle[i]));
    }
    if ((areas[0] < 0 && areas[1] < 0 && areas[2] < 0) || (areas[0] > 0 && areas[1] > 0 && areas[2] > 0)) {
        w = divide<T, 3>(areas, areas[0]+areas[1]+areas[2]);
        return true;
    }
    return false;
}

template<typename T, int d, int d2>
array<T, d> elementCenter(const vector<array<T, d>>& vertices, const array<int, d2>& element) {
    array<T, d> center;
    for (auto i : element) {
        for (int j = 0; j < d; ++j) {
            center[j] += vertices[i][j];
        }
    }
    for (int i = 0; i < d; ++i) {
        center[i] /= (T)element.size();
    }
    return center;
}

template<typename T, int d>
struct Box {
    array<T, d> lowerLeft_, upperRight_;
    
    Box() {}
    
    Box(const Box& b) {
        for (int i = 0; i < d; ++i) {
            lowerLeft_[i] = b.lowerLeft_[i];
            upperRight_[i] = b.upperRight_[i];
        }
    }
    
    Box(const Box& b1, const Box& b2) {
        for (int i = 0; i < d; ++i) {
            lowerLeft_[i] = min(b1.lowerLeft_[i], b2.lowerLeft_[i]);
            upperRight_[i] = max(b1.upperRight_[i], b2.upperRight_[i]);
        }
    }
    
    bool intersects(const Box& b) {
        for (int i = 0; i < d; ++i) {
            if (lowerLeft_[i] > b.upperRight_[i] || upperRight_[i] < b.lowerLeft_[i]) {
                return false;
            }
        }
        return true;
    }
};

template<typename T, int d, int d1>
Box<T,d> buildBox(const vector<array<T, d>>& vertices, const array<int, d1>& element) {
    Box<T,d> b;
//    print<int,d1>(element);
    for (int i = 0; i < d; ++i) {
        b.lowerLeft_[i] = numeric_limits<T>::max();
        b.upperRight_[i] = numeric_limits<T>::lowest();
    }
//    print<T,d>(b.lowerLeft_);
//    print<T,d>(b.upperRight_);
    for (size_t i = 0; i < d1; ++i) {
//        print<T,d>(vertices[element[i]]);
        for (int j = 0; j < d; ++j) {
            b.lowerLeft_[j] = min(b.lowerLeft_[j], vertices[element[i]][j]);
            b.upperRight_[j] = max(b.upperRight_[j], vertices[element[i]][j]);
        }
    }
//    print<T,d>(b.lowerLeft_);
//    print<T,d>(b.upperRight_);
    return b;
}

template<typename T, int d>
struct BoxNode {
    int n_;
    BoxNode *left_, *right_;
    Box<T, d> box_;
    
    BoxNode(int n, const vector<Box<T, d>>& boxes) : n_(n), left_(nullptr), right_(nullptr), box_(boxes[n]) {}
    
    BoxNode(BoxNode* left, BoxNode* right) : n_(-1), left_(left), right_(right), box_(left->box_, right->box_) {}
    
    ~BoxNode() {
        if (left_ != nullptr) {
            delete left_;
        }
        if (right_ != nullptr) {
            delete right_;
        }
    }
};

template<typename T, int d>
BoxNode<T, d>* buildBoxHierarchy(const vector<Box<T, d>>& boxes, const vector<array<T, d>>& centers, vector<size_t>& elementIndexes, int begin, int end, int level) {
    BoxNode<T, d>* root = nullptr;
    if (elementIndexes.size() == 0) {
        return nullptr;
    }
    if (begin == end) {
        root = new BoxNode<T, d>(elementIndexes[begin], boxes);
    } else {
        nth_element(elementIndexes.begin()+begin, elementIndexes.begin()+(begin+end)/2, elementIndexes.begin()+end, [&](std::size_t i, std::size_t j){ return centers[i][level%d] < centers[j][level%d]; });
        BoxNode<T, d> *left = buildBoxHierarchy<T, d>(boxes, centers, elementIndexes, begin, (begin+end)/2, ++level);
        BoxNode<T, d> *right = buildBoxHierarchy<T, d>(boxes, centers, elementIndexes, (begin+end)/2+1, end, level);
        root = new BoxNode<T, d>(left, right);
    }
    return root;
}

template<typename T, int d>
class BoxHierarchy {
    int n_; //number of boxes
    BoxNode<T, d>* root_;
public:
    BoxHierarchy(const vector<Box<T, d>>& boxes, const vector<array<T, d>>& centers) : n_(centers.size()) {
        if (boxes.size()) {
            vector<size_t> elementIndexes(boxes.size());
            iota(elementIndexes.begin(), elementIndexes.end(), 0);
            root_ = buildBoxHierarchy<T, d>(boxes, centers, elementIndexes, 0, boxes.size()-1, 0);
        }
    }
    
    void intersect(const BoxHierarchy<T, d>& bh, vector<vector<int>>& intersectingElements) const {
        intersectingElements.clear();
        intersectingElements.resize(n_);
        if (!root_ || !bh.root_) {
            return;
        }
        stack<pair<BoxNode<T, d>*, BoxNode<T, d>*>> s;
        s.push(pair<BoxNode<T, d>*, BoxNode<T, d>*>(root_, bh.root_));
        while (s.size()) {
            pair<BoxNode<T, d>*, BoxNode<T, d>*> top = s.top();
            s.pop();
            if (top.first->box_.intersects(top.second->box_)) {
                if (top.first->n_ != -1 && top.second->n_ != -1) {
                    intersectingElements[top.first->n_].push_back(top.second->n_);
                } else if (top.second->n_ != -1) {
                    s.push(pair<BoxNode<T, d>*, BoxNode<T, d>*>(top.first->left_, top.second));
                    s.push(pair<BoxNode<T, d>*, BoxNode<T, d>*>(top.first->right_, top.second));
                } else if (top.first->n_ != -1) {
                    s.push(pair<BoxNode<T, d>*, BoxNode<T, d>*>(top.first, top.second->left_));
                    s.push(pair<BoxNode<T, d>*, BoxNode<T, d>*>(top.first, top.second->right_));
                } else {
                    s.push(pair<BoxNode<T, d>*, BoxNode<T, d>*>(top.first->left_, top.second->left_));
                    s.push(pair<BoxNode<T, d>*, BoxNode<T, d>*>(top.first->left_, top.second->right_));
                    s.push(pair<BoxNode<T, d>*, BoxNode<T, d>*>(top.first->right_, top.second->left_));
                    s.push(pair<BoxNode<T, d>*, BoxNode<T, d>*>(top.first->right_, top.second->right_));
                }
            }
        }
    }
    
    ~BoxHierarchy() {
        delete root_;
    }
};

template<typename T, int d1, int d2>
BoxHierarchy<T, d1> buildBoxHierarchy(const vector<array<T,d1>>& nodes, const vector<array<int, d2>>& elements) {
    vector<Box<T, d1>> boxes;
    vector<array<T, d1>> centers;
    for (const auto& e: elements) {
        boxes.push_back(buildBox<T,d1,d2>(nodes, e));
        centers.push_back(elementCenter<T,d1,d2>(nodes, e));
    }
    return BoxHierarchy<T, d1>(boxes, centers);
}

template<typename T>
class TriMesh {
    typedef std::array<int, 2> I2;
    typedef std::array<int, 3> I3;
    typedef std::array<int, 4> I4;
    typedef std::array<T, 3> TV;
    
public:
    vector<TV> nodes_;
    vector<I3> mesh_;
    
    void clear() {
        nodes_.clear();
        mesh_.clear();
    }
};

const array<array<int, 3>, 4> FaceIndexes = {
    array<int, 3>{0,1,2},
    array<int, 3>{0,2,3},
    array<int, 3>{1,2,3},
    array<int, 3>{0,3,1}
};

array<array<int,3>,4> tetFaces(const array<int,4>& tet) {
    array<array<int,3>,4> faces;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 3; ++j) {
            faces[i][j] = tet[FaceIndexes[i][j]];
        }
    }
    return faces;
}

array<int,3> tetFace(const array<int,4>& tet, int i) {
    array<int,3> face;
    for (int j = 0; j < 3; ++j) {
        face[j] = tet[FaceIndexes[i][j]];
    }
    return face;
}

const array<array<int, 2>, 6> EdgeIndexes = {array<int,2>{0,1}, array<int,2>{0,2}, array<int,2>{0,3}, array<int,2>{1,2}, array<int,2>{1,3}, array<int,2>{2,3}};

array<array<int,2>,6> tetEdges(const array<int,4>& tet) {
    array<array<int,2>,6> faces;
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 2; ++j) {
            faces[i][j] = tet[EdgeIndexes[i][j]];
        }
    }
    return faces;
}

array<array<int,2>,3> faceEdges(const array<int,3>& face) {
    array<array<int,2>,3> edges;
    for (int i = 0; i < 3; ++i) {
        edges[i] = array<int,2>{face[i], face[(i+1)%3]};
    }
    return edges;
}

template<typename T, int d, int d1>
array<array<T,d>,d1> elementNodes(const vector<array<T,d>>& nodes, const array<int,d1>& element) {
    array<array<T,d>,d1> ps;
    for (int i = 0; i < d1; ++i) {
        ps[i] = nodes[element[i]];
    }
    return ps;
}
template<typename T>
class TetMesh {
    typedef array<int, 2> I2;
    typedef array<int, 3> I3;
    typedef array<int, 4> I4;
    typedef array<T, 3> TV;

public:
    vector<TV> nodes_;
    vector<I4> mesh_;
    vector<I3> surfaceMesh_;
    vector<int> connectedComponents_; //connected component id of each element

    TetMesh() {}

    TetMesh(vector<TV>&& nodes, vector<I4>&& mesh): nodes_(nodes), mesh_(mesh) {
        initializeSurfaceMesh();
        computeConnectedComponents();
    }
    void initializeSurfaceMesh() {
        surfaceMesh_.clear();
        map<I3, I3> surfaceElements; //sorted to unsorted elements
        for (const auto& tet: mesh_) {
            for (const auto& fi: FaceIndexes) {
                auto face = I3{tet[fi[0]], tet[fi[1]], tet[fi[2]]};
                auto sortedFace = face;
                sort(sortedFace.begin(), sortedFace.end());
                if (surfaceElements.count(sortedFace)) {
                    surfaceElements.erase(sortedFace);
                } else {
                    surfaceElements[sortedFace] = face;
                }
            }
        }
        for (const auto& e: surfaceElements) {
            surfaceMesh_.push_back(e.second);
        }
    }
    
    void computeConnectedComponents() {
        connectedComponents_.resize(mesh_.size());
        UnionFind nodeClasses(nodes_.size());
        for (int i = 0; i < mesh_.size(); ++i) {
            for (int j = 0; j < 3; ++j) {
                nodeClasses.merge(mesh_[i][j], mesh_[i][j+1]);
            }
        }
        map<int, int> nodeClassToCC;
        int c = 1;
        for (int i = 0; i < mesh_.size(); ++i) {
            if (!nodeClassToCC.count(nodeClasses.find(mesh_[i][0]))) {
                connectedComponents_[i] = c;
                nodeClassToCC[nodeClasses.find(mesh_[i][0])] = c;
                ++c;
            } else {
                connectedComponents_[i] = nodeClassToCC[nodeClasses.find(mesh_[i][0])];
            }
        }
        cout << "found " << c-1 << " connected components\n";
    }
};

#endif /* DataStructures_hpp */

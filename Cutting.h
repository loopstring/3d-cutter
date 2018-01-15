//
//  Copyright (c) 2018 Yuting Wang. All rights reserved.
//
#ifndef Cutting_h
#define Cutting_h

#include "DataStructures.h"

template<typename T>
class Cutter3D {
    typedef std::array<int, 1> I1;
    typedef std::array<int, 2> I2;
    typedef std::array<int, 3> I3;
    typedef std::array<int, 4> I4;
    typedef std::array<int, 5> I5;
    typedef std::array<T, 2> T2;
    typedef std::array<T, 3> T3;
    typedef std::array<T, 4> T4;
    typedef map<I4, T4> Intersections;
    typedef map<I4, vector<int>> TetBoundary2TetIds;

    struct CutElement {
        int parentElementIndex;
        array<bool, 4> subElements; // in the same order as the tet nodes
        
        CutElement(int i, bool fill = true): parentElementIndex(i) {
            subElements.fill(fill);
        }

        int numPieces() const {
            return (int)subElements[0] + (int)subElements[1] + (int)subElements[2] + (int)subElements[3];
        }
    };
    
    static bool computeIntersection(const array<T3,2>& nodes1, const array<T3,3>& nodes2, array<T,2>& w1, array<T, 3>& w2) {
        T v1 = volume<T>(nodes1[0], nodes2[0], nodes2[1], nodes2[2]);
        T v2 = volume<T>(nodes1[1], nodes2[0], nodes2[1], nodes2[2]);
        T v3 = volume<T>(nodes1[0], nodes1[1], nodes2[0], nodes2[1]);
        T v4 = volume<T>(nodes1[0], nodes1[1], nodes2[1], nodes2[2]);
        T v5 = volume<T>(nodes1[0], nodes1[1], nodes2[2], nodes2[0]);
        if (v1*v2<0 && (v3>0)==(v4>0) && (v4>0)==(v5>0)) {
            w1[0] = fabs(v2) / (fabs(v1) + fabs(v2));
            w1[1] = 1 - w1[0];
            T v = fabs(v3) + fabs(v4) + fabs(v5);
            w2[0] = fabs(v4) / v;
            w2[1] = fabs(v5) / v;
            w2[2] = 1 - w2[0] - w2[1];
            //cout << w1[0] << endl;
            return true;
        } else {
            return false;
        }
    }

    static bool computeIntersection(const array<T3,2>& nodes1, const array<T3,3>& nodes2, array<T,2>& w) {
        array<T,3> w1;
        return computeIntersection(nodes1, nodes2, w, w1);
    }

    static bool computeIntersection(const array<T3,3>& nodes1, const array<T3,2>& nodes2, array<T,3>& w) {
        array<T,2> w1;
        return computeIntersection(nodes2, nodes1, w1, w);
    }

    static bool computeIntersection(const array<T3,4>& nodes1, const array<T3,1>& nodes2, array<T,4>& w) {
        T v1 = volume<T>(nodes1[0], nodes1[1], nodes1[2], nodes2[0]);
        T v2 = volume<T>(nodes1[0], nodes1[2], nodes1[3], nodes2[0]);
        T v3 = volume<T>(nodes1[0], nodes1[3], nodes1[1], nodes2[0]);
        T v4 = volume<T>(nodes2[0], nodes1[1], nodes1[2], nodes1[3]);
        if (v1 == 0 || v2 == 0 || v3 == 0 || v4 == 0) {
            cout << "point tet degenerate case" << endl;
        }
        // cout << v1 << ", " << v2 << ", " << v3 << ", " << v4 << endl;
        if ((v1>0) == (v2>0) && (v2>0) == (v3>0) && (v3>0) == (v4>0)) {
            T v = fabs(v1) + fabs(v2) + fabs(v3) + fabs(v4);
            w[0] = fabs(v4) / v;
            w[1] = fabs(v2) / v;
            w[2] = fabs(v3) / v;
            w[3] = 1 - w[0] - w[1] - w[2];
            return true;
        } else {
            return false;
        }
        return false;
    }

    template<int d1, int d2>
    static void computeIntersections(const vector<T3>& nodes1, const vector<T3>& nodes2, const vector<array<int,d1>>& e1, const vector<array<int,d2>>& e2, const BoxHierarchy<T,3>& b1, const BoxHierarchy<T,3>& b2, map<I4, T4>& intersections) {
        vector<vector<int>> intersectingBoxes; // intersecting boxes
        b1.intersect(b2, intersectingBoxes);
        for (size_t i = 0; i < intersectingBoxes.size(); ++i) {
            //cout << "e1 " << i << endl;
            //print<int,d1>(e1[i]);
            for (auto j : intersectingBoxes[i]) {
                //cout << j << ", " << endl;
                //print<int,d2>(e2[j]);
                auto tetNodes = elementNodes<T,3,d1>(nodes1, e1[i]);
                auto triNodes = elementNodes<T,3,d2>(nodes2, e2[j]);
                array<T,d1> w;
                if (computeIntersection(tetNodes, triNodes, w)) {
                    intersections[toI4<int,d1>(e1[i])] = toI4<T,d1>(w,0);
                }
            }
        }
    }

    static Intersections computeIntersections(const TetMesh<T>& tetMesh, const TriMesh<T>& triMesh, TetBoundary2TetIds& tetBoundary2TetIds) {
        map<I4, T4> intersections;
        
        // build box hierarchies for tetMesh
        set<I3> tetMeshFaces;
        set<I2> tetMeshEdges;
        for (int i = 0; i < tetMesh.mesh_.size(); ++i) {
            auto tet = tetMesh.mesh_[i];
            sort(tet.begin(), tet.end());
            tetBoundary2TetIds[tet].push_back(i);
            auto faces = tetFaces(tet);
            for (auto& face: faces) {
                sort(face.begin(), face.end());
                tetBoundary2TetIds[toI4<int,3>(face)].push_back(i);
                tetMeshFaces.insert(face);
            }
            auto edges = tetEdges(tet);
            for (auto& edge: edges) {
                sort(edge.begin(), edge.end());
                tetBoundary2TetIds[toI4<int,2>(edge)].push_back(i);
                tetMeshEdges.insert(edge);
            }
        }
        vector<I1> tetMeshNodeVec;
        for (int i = 0; i < tetMesh.nodes_.size(); ++i) {
            tetMeshNodeVec.push_back(I1{i});
        }
        vector<I3> tetMeshFaceVec(tetMeshFaces.begin(), tetMeshFaces.end());
        vector<I2> tetMeshEdgeVec(tetMeshEdges.begin(), tetMeshEdges.end());
        cout << "buliding tet mesh hierarchy" << endl;
        auto tetMeshHierarchy = buildBoxHierarchy<T,3,4>(tetMesh.nodes_, tetMesh.mesh_);
        auto tetMeshFaceHierarchy = buildBoxHierarchy<T,3,3>(tetMesh.nodes_, tetMeshFaceVec);
        auto tetMeshEdgeHierarchy = buildBoxHierarchy<T,3,2>(tetMesh.nodes_, tetMeshEdgeVec);
        auto tetMeshNodeHierarchy = buildBoxHierarchy<T,3,1>(tetMesh.nodes_, tetMeshNodeVec);
        cout << "tet mesh hierarchy built" << endl;

        // box hierarchy for triMesh
        set<I2> triMeshEdges;
        for (const auto& tri: triMesh.mesh_) {
            auto edges = faceEdges(tri);
            for (auto& edge: edges) {
                sort(edge.begin(), edge.end());
                triMeshEdges.insert(edge);
            }
        }
        vector<I2> triMeshEdgeVec(triMeshEdges.begin(), triMeshEdges.end());
        vector<I1> triMeshNodeVec;
        for (int i = 0; i < triMesh.nodes_.size(); ++i) {
            triMeshNodeVec.push_back(I1{i});
        }
//        cout << "trimeshhierarchy\n";
//        print<T,3>(triMesh.nodes_);
//        print<int,3>(triMesh.mesh_);
        auto triMeshHierarchy = buildBoxHierarchy<T,3,3>(triMesh.nodes_, triMesh.mesh_);
//        cout << "trimeshhierarchy\n";
        auto triMeshEdgeHierarchy = buildBoxHierarchy<T,3,2>(triMesh.nodes_, triMeshEdgeVec);
        auto triMeshNodeHierarchy = buildBoxHierarchy<T,3,1>(triMesh.nodes_, triMeshNodeVec);
        cout << "tri mesh hierarchy built" << endl;

        // compute intersections
        // v-v
        // v-e
        // v-f
        // e-v
        // e-e
        // e-f
        computeIntersections<2,3>(tetMesh.nodes_, triMesh.nodes_, tetMeshEdgeVec, triMesh.mesh_, tetMeshEdgeHierarchy, triMeshHierarchy, intersections);
        // f-v
        // f-e
        computeIntersections<3,2>(tetMesh.nodes_, triMesh.nodes_, tetMeshFaceVec, triMeshEdgeVec, tetMeshFaceHierarchy, triMeshEdgeHierarchy, intersections);
        // t-v
        computeIntersections<4,1>(tetMesh.nodes_, triMesh.nodes_, tetMesh.mesh_, triMeshNodeVec, tetMeshHierarchy, triMeshNodeHierarchy, intersections);

        return intersections;
    }

    static vector<CutElement> split(const TetMesh<T>& tetMesh, const Intersections& intersections, TetBoundary2TetIds& tetBoundary2TetIds, set<int>& cutTets) {
        cutTets.clear();
        for (const auto& t: tetBoundary2TetIds) {
            if (intersections.count(t.first)) {
                for (auto i: t.second) {
                    cutTets.insert(i);
                }
            }
        }
        cout << cutTets.size() << " tets cut\n";
        vector<CutElement> v;
        for (int i = 0; i < tetMesh.mesh_.size(); ++i) {
            if (cutTets.count(i)) {
                array<bool,4> added;
                added.fill(false);
                auto tet = tetMesh.mesh_[i];
                for (int j = 0; j < 4; ++j) {
                    if (!added[j]) {
                        // find all connected pieces
                        CutElement ce(i, false);
                        stack<int> s;
                        s.push(j);
                        while(s.size()) {
                            auto top = s.top();
                            ce.subElements[top] = true;
                            added[top] = true;
                            s.pop();
                            // add all the connected pieces that are not added yet
                            for (int k = 0; k < 4; ++k) {
                                if (!added[k]) {
                                    if (!intersections.count(toI4<int,2>(sorted(I2{tet[top],tet[k]})))) {
                                        s.push(k);
                                    }
                                }
                            }
                        }
                        v.push_back(ce);
                    }
                }
            }
        }
        return v;
    }

    void static newTet(int parentId, const I4& tet, const TetMesh<T>& tetMesh, vector<T3>& newNodes, vector<I4>& newMesh, map<int,int>& nodeMapping, UnionFind& uf) {
        I4 newTet;
        //cout << "parent id " << parentId << endl;
        for (int i = 0; i < 4; ++i) { // for each node
            int newId = uf.find(tet[i]);
            //cout << tet[i] << ", " << newId << endl;
            const auto& it = nodeMapping.find(newId);
            if (it != nodeMapping.end()) {
                newTet[i] = it->second;
            } else {
                newTet[i] = newNodes.size();
                nodeMapping[newId] = newNodes.size();
                newNodes.push_back(tetMesh.nodes_[tetMesh.mesh_[parentId][i]]);
            }
        }

        newMesh.push_back(newTet);
    }
    
    static void merge(const vector<CutElement>& cutElements, const TetMesh<T>& tetMesh, vector<T3>& newNodes, vector<I4>& newMesh, const Intersections& intersections) {
        newNodes.clear();
        newMesh.clear();
        UnionFind uf(tetMesh.nodes_.size() + 4 * cutElements.size());
        map<I5, int> faceNode2NewNode; // key = {face,materialNode,node}
        set<int> cutTets;
        int total = tetMesh.nodes_.size();
        for (const auto& ce: cutElements) { // need to do face-face merging even for tets that are touched by the cut but not split, so that if a neighbor splits they are all connected to it.
            cutTets.insert(ce.parentElementIndex);
            const auto& tet = tetMesh.mesh_[ce.parentElementIndex];
            for (int i = 0; i < 4; ++i) { // for each face
                auto face = tetFace(tet, i);
                sort(face.begin(), face.end());
                I5 key;
                for (int j = 0; j < 3; ++j) {
                    key[j] = face[j];
                }
                for (int j = 0; j < 3; ++j) { // for each node check for material
                    int fij = FaceIndexes[i][j];
                    if (ce.subElements[fij]) {
                        key[3] = tet[fij];
                        uf.merge(total+fij, key[3]);
                        for (int k = 0; k < 3; ++k) { // for each node, merge
                            int fik = FaceIndexes[i][k];
                            key[4] = tet[fik];
                            int newId = total+fik;
                            //print<int,5>(key);
                            const auto& it = faceNode2NewNode.find(key);
                            if (it != faceNode2NewNode.end()) {
                                //cout << "merging " << it->second << ", " << newId << endl;
                                uf.merge(it->second, newId);
                            } else {
                                faceNode2NewNode[key] = newId;
                            }
                        }
                    }
                }
            }
            total += 4;
        }
        total = tetMesh.nodes_.size();
        map<int,int> nodeMapping;
        for (const auto& ce: cutElements) {
            newTet(ce.parentElementIndex, I4{total, total+1, total+2, total+3}, tetMesh, newNodes, newMesh, nodeMapping, uf);
            total += 4;
        }
        for (int i = 0; i < tetMesh.mesh_.size(); ++i) {
            if (!cutTets.count(i)) {
                newTet(i, tetMesh.mesh_[i], tetMesh, newNodes, newMesh, nodeMapping, uf);
            }
        }

//        cout << "merged mesh \n";
//        print<T,3>(newNodes);
//        print<int,4>(newMesh);
    }
    
    static TetMesh<T> subdivide(const vector<CutElement>& cutElements, const TetMesh<T>& tetMesh, vector<T3>& newNodes, vector<I4>& newMesh, Intersections& intersections) {
        // add a new node inside the tet, connect with cuts on each face to subdivide the tet
        map<I4, int> newNodeMapping;
        for (int i = 0; i < cutElements.size(); ++i) {
            const auto& ce = cutElements[i];
            const auto& originalTet = tetMesh.mesh_[ce.parentElementIndex];
            const auto sortedOriginalTet = sorted(originalTet);
            const auto& tet = newMesh[i];
            
            // get all edge cuts and add them as new nodes
            const auto originalEdges = tetEdges(originalTet);
            const auto edges = tetEdges(tet);
            int cutEdges = 0;
            T4 averageEdgeWeight{0,0,0,0};
            map<int, T> originalNodeId2Weight;
            for (int k = 0; k < originalEdges.size(); ++k) {
                auto sortedOriginalEdge = toI4<int,2>(sorted(originalEdges[k]));
                auto sortedEdge = toI4<int,2>(sorted(edges[k]));
                const auto& it = intersections.find(sortedOriginalEdge);
                if (it != intersections.end()) {
                    ++cutEdges;
                    for (int j = 0; j < 2; ++j) {
                        originalNodeId2Weight[sortedOriginalEdge[j]] += it->second[j];
                    }
                    const auto& idIt = newNodeMapping.find(sortedEdge);
                    if (idIt == newNodeMapping.end()) {
                        newNodeMapping[sortedEdge] = newNodes.size();
                        newNodes.push_back(elementCenter<T,3>(tetMesh.nodes_, sortedOriginalEdge, it->second));
//                        cout << "edge node ";
//                        print<T,3>(elementCenter<T,3>(tetMesh.nodes_, sortedOriginalEdge, it->second));
                    }
                }
            }
            for (int j = 0; j < 4; ++j) {
                averageEdgeWeight[j] = originalNodeId2Weight[sortedOriginalTet[j]];
            }
            //cout << "cutEdges " << cutEdges << endl;

            // face cuts
            const auto originalFaces = tetFaces(originalTet);
            const auto faces = tetFaces(tet);
            for (int k = 0; k < faces.size(); ++k) {
                auto sortedOriginalFace = toI4<int,3>(sorted(originalFaces[k]));
                auto sortedFace = toI4<int,3>(sorted(faces[k]));
                const auto& it = intersections.find(sortedOriginalFace);
                if (it != intersections.end()) { // face center already computed
                    const auto& idIt = newNodeMapping.find(sortedFace);
                    if (idIt == newNodeMapping.end()) {
                        newNodeMapping[sortedFace] = newNodes.size();
                        newNodes.push_back(elementCenter<T,3>(tetMesh.nodes_, sortedOriginalFace, it->second));
                    }
//                    cout << "face center ";
//                    print<T,3>(elementCenter<T,3>(tetMesh.nodes_, sortedOriginalFace, it->second));
                } else { // use average of edge cuts if not
                    int numEdges = 0;
                    T4 faceWeights{0,0,0,0};
                    map<int, T> node2weight;
                    for (int j = 0; j < 3; ++j) {
                        auto sortedOriginalEdge = toI4<int,2>(sorted(array<int,2>{sortedOriginalFace[j], sortedOriginalFace[(j+1)%3]}));
                        const auto& edgeIt = intersections.find(sortedOriginalEdge);
                        if (edgeIt != intersections.end()) {
                            ++numEdges;
                            for (int e = 0; e < 2; ++e) {
                                node2weight[sortedOriginalEdge[e]] += edgeIt->second[e];
                            }
                        }
                    }
                    if (numEdges > 1) { // otherwise don't add new face center
                        newNodeMapping[sortedFace] = newNodes.size();
                        for (int j = 0; j < 3; ++j) {
                            faceWeights[j] = node2weight[sortedOriginalFace[j]] / numEdges;
                        }
//                        cout << "face weight ";
//                        print<T,4>(faceWeights);
//                        cout << "face center ";
//                        print<T,3>(elementCenter<T,3>(tetMesh.nodes_, sortedOriginalFace, faceWeights));
                        newNodes.push_back(elementCenter<T,3>(tetMesh.nodes_, sortedOriginalFace, faceWeights));
                        intersections[sortedOriginalFace] = faceWeights;
                    }
                }
            }
            
            // tet center
            int tetCenterId = newNodes.size();
            const auto& tetCenterIt = intersections.find(sortedOriginalTet);
            if (tetCenterIt != intersections.end()) {
                newNodes.push_back(elementCenter<T,3>(tetMesh.nodes_, sortedOriginalTet, tetCenterIt->second));
//                cout << "tet center ";
//                print<T,3>(elementCenter<T,3>(tetMesh.nodes_, sortedOriginalTet, tetCenterIt->second));
            } else { // if doesn't exist, use average of edge cuts or the center
                if (ce.numPieces() == 4) {
                    averageEdgeWeight.fill(0.25);
                } else {
                    averageEdgeWeight = divide<T,4>(averageEdgeWeight, cutEdges);
//                    print<T,4>(averageEdgeWeight);
                }
                newNodes.push_back(elementCenter<T,3>(tetMesh.nodes_, sortedOriginalTet, averageEdgeWeight));
//                cout << "tet center ";
//                print<T,3>(elementCenter<T,3>(tetMesh.nodes_, sortedOriginalTet, averageEdgeWeight));
                intersections[sortedOriginalTet] = averageEdgeWeight;
            }

            // add elements that are created by the new nodes added above
            vector<I4> newTets;
            for (int f = 0; f < faces.size(); ++f) {
                const auto& face = faces[f];
                const auto sortedFace = toI4<int,3>(sorted(face));
                const auto& newFaceCenterIt = newNodeMapping.find(sortedFace);
                if (newFaceCenterIt != newNodeMapping.end()) {
                    for (int j = 0; j < 3; ++j) {
                        auto sortedEdge = toI4<int,2>(sorted(array<int,2>{face[j], face[(j+1)%3]}));
                        const auto& newEdgeCenterIt = newNodeMapping.find(sortedEdge);
                        if (newEdgeCenterIt != newNodeMapping.end()) {
                            if (ce.subElements[FaceIndexes[f][j]]) {
                                newTets.push_back(I4{tetCenterId, newFaceCenterIt->second, face[j], newEdgeCenterIt->second});
                            }
                            if (ce.subElements[FaceIndexes[f][(j+1)%3]]) {
                                newTets.push_back(I4{tetCenterId, newFaceCenterIt->second, newEdgeCenterIt->second, face[(j+1)%3]});
                            }
                        } else if (ce.subElements[FaceIndexes[f][j]]) {
                            newTets.push_back(I4{tetCenterId, newFaceCenterIt->second, face[j], face[(j+1)%3]});
                        }
                    }
                } else if (ce.subElements[FaceIndexes[f][0]]) { // no face intersection, might have 0 or 1 edge cut
                    bool isSplit = false;
                    for (int j = 0; j < 3; ++j) {
                        auto sortedEdge = toI4<int,2>(sorted(array<int,2>{face[j], face[(j+1)%3]}));
                        const auto& newEdgeCenterIt = newNodeMapping.find(sortedEdge);
                        if (newEdgeCenterIt != newNodeMapping.end()) {
                            newTets.push_back(I4{tetCenterId, face[(j+2)%3], face[j], newEdgeCenterIt->second});
                            newTets.push_back(I4{tetCenterId, face[(j+2)%3], newEdgeCenterIt->second, face[(j+1)%3]});
                            isSplit = true;
                            break;
                        }
                    }
                    if (!isSplit) {
                        newTets.push_back(I4{tetCenterId, face[0], face[1], face[2]});
                    }
                }
            }
            newMesh[i] = newTets[0];
            for (int j = 1; j < newTets.size(); ++j) {
                newMesh.push_back(newTets[j]);
            }
        }
        return TetMesh<T>(move(newNodes), move(newMesh));
    }

public:
    static TetMesh<T> run(const TetMesh<T>& tetMesh, const TriMesh<T>& triMesh) {
        TetBoundary2TetIds tetBoundary2TetIds;
        auto intersections = computeIntersections(tetMesh, triMesh, tetBoundary2TetIds);
        cout << "finished computing " << intersections.size() << " intersections\n";
//        for (auto& a: intersections) {
//            print<int,4>(a.first);
//            print<T,4>(a.second);
//        }
        set<int> cutTets;
        vector<CutElement> cutElements = split(tetMesh, intersections, tetBoundary2TetIds, cutTets);
//        for (auto& ce: cutElements) {
//            cout << ce.parentElementIndex << endl;
//            print<bool,4>(ce.subElements);
//        }
        vector<T3> newNodes;
        vector<I4> newMesh;
        merge(cutElements, tetMesh, newNodes, newMesh, intersections);
        cout << "finished split-merge\n";
        return subdivide(cutElements, tetMesh, newNodes, newMesh, intersections);
    }
};

#endif /* Cutting_h */

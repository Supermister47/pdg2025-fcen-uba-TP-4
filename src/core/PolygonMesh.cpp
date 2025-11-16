//------------------------------------------------------------------------
//  Copyright (C) Gabriel Taubin
//  Time-stamp: <2025-08-05 23:13:04 taubin>
//------------------------------------------------------------------------
//
// PolygonMesh.cpp
//
// Software developed for the course
// Digital Geometry Processing
// Copyright (c) 2025, Gabriel Taubin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
//     * Redistributions of source code must retain the above
//       copyright notice, this list of conditions and the following
//       disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials
//       provided with the distribution.
//     * Neither the name of the Brown University nor the names of its
//       contributors may be used to endorse or promote products
//       derived from this software without specific prior written
//       permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL GABRIEL
// TAUBIN BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
// OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
// USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.

#include <iostream>
#include "PolygonMesh.hpp"
#include "Partition.hpp"
#include "Graph.hpp"

// TODO Mon Mar 6 2023
// - merge your code from Assignment 2

PolygonMesh::PolygonMesh(const int nVertices, const vector<int>& coordIndex) : HalfEdges(nVertices, coordIndex),
_nPartsVertex(),
_isBoundaryVertex(),
_dualGraph(getNumberOfFaces()) {
  int nV = getNumberOfVertices();
  int nE = getNumberOfEdges(); // Edges method
  int nC = getNumberOfCorners();

  // 1) classify the vertices as boundary or internal
  _createIsBoundaryVertex(nV, nE);
  // 2) create a partition of the corners in the stack
  Partition partition(nC);
  // 3) for each regular edge
  //    - get the two half edges incident to the edge
  //    - join the two pairs of corresponding corners accross the edge
  //    - you need to take into account the relative orientation of
  //      the two incident half edges
  for (int iE = 0; iE < nE; iE++) {
    int nFaces = getNumberOfEdgeHalfEdges(iE);

    // If there is exactly two faces incident to the edge
    // then the edge is regular and we can join the corners
    // accross the edge
    if (nFaces == 2) {

      int iC00 = getEdgeHalfEdge(iE, 0);
      int iC01 = getDst(iC00);

      int iC10 = getEdgeHalfEdge(iE, 1);
      int iC11 = getDst(iC10);

      if (HalfEdges::isOriented(iC00)) {
        // Consistently oriented
        partition.join(iC00, iC11);
        partition.join(iC01, iC10);
      }
      else {
        // Oposite orientation
        partition.join(iC00, iC10);
        partition.join(iC01, iC11);
      }
    }
  }

  // consistently oriented
  /* \                  / */
  /*  \ iC01 <-- iC00  /  */
  /*   X ---- iE ---- X   */
  /*  / iC10 --> iC11  \  */
  /* /                  \ */

  // oposite orientation
  /* \                  / */
  /*  \ iC01 --> iC00  /  */
  /*   X ---- iE ---- X   */
  /*  / iC10 --> iC11  \  */
  /* /                  \ */

  // a decision has to be made about inconsistently oriented faces
  // incident to the same edge, as well as how to deal with singular
  // edges; for the moment let's assume that the mesh does not have
  // singular edges, and that pairs of corners corresponding to the
  // same vertex across inconsistently oriented faces will be joined

  // note that the partition will end up with the corner separators as
  // singletons, but it doesn't matter for the last step, and
  // the partition will be deleteted upon return

  // 4) count number of parts per vertex
  //    - initialize _nPartsVertex array to 0's
  //    - for each corner iC which is a representative of its subset,
  //    - get the corresponding vertex index iV and increment _nPartsVertex[iV]
  //    - note that all the corners in each subset share a common
  //      vertex index, but multiple subsets may correspond to the
  //      same vertex index, indicating that the vertex is singular

  _create_nPartsVertex(nVertices, nC, partition);

  // Create the dual graph

  _createDualGraph();
}

void PolygonMesh::_createIsBoundaryVertex(int nV, int nE) {
  _isBoundaryVertex.resize(nV, false);

  // - for edge boundary iE label its two end vertices as boundary
  for (int iE = 0; iE < nE; iE++) {
    int nFaces = getNumberOfEdgeHalfEdges(iE);

    // If there is at least one edge that is regular or singular
    // then both vertices of the edge can not be boundary
    if (nFaces == 1) {
      _isBoundaryVertex[getVertex0(iE)] = true;
      _isBoundaryVertex[getVertex1(iE)] = true;
    }
  }
}

void PolygonMesh::_create_nPartsVertex(const int nVertices, int nC, Partition& partition) {
  _nPartsVertex.resize(nVertices, 0);

  for (int iC = 0; iC < nC; iC++) {
    if (partition.find(iC) == iC) {
      int iV = _coordIndex[iC];

      if (iV == -1)
        continue; // skip face separators

      _nPartsVertex[iV]++;
    }
  }
}

void PolygonMesh::_createDualGraph() {
  for (int iE = 0; iE < getNumberOfEdges(); iE++) {
    int nFaces = getNumberOfEdgeFaces(iE);

    // For each face adjacent to edge iE
    for (int j0 = 0; j0 < nFaces; j0++) {
      int iF0 = getEdgeFace(iE, j0);

      // Connect it to the others faces that are adjacent to edge iE
      for (int j1 = j0 + 1; j1 < nFaces; j1++) {
        int iF1 = getEdgeFace(iE, j1);
        _dualGraph.insertEdge(iF0, iF1);
      }

    }
  }
}

bool PolygonMesh::_isValidFace(const int iF) const {
  return 0 <= iF && iF < getNumberOfFaces();
}

int PolygonMesh::getNumberOfFaces() const {
  return _nF;
}

int PolygonMesh::getNumberOfEdgeFaces(const int iE) const {
  return getNumberOfEdgeHalfEdges(iE);
}

int PolygonMesh::getEdgeFace(const int iE, const int j) const {
  if (!(_isValidEdge(iE) && 0 <= j && j < getNumberOfEdgeFaces(iE))) {
    return -1;
  }

  int iHE = getEdgeHalfEdge(iE, j);
  return getFace(iHE);
}

bool PolygonMesh::isEdgeFace(const int iE, const int iF) const {
  if (!(_isValidEdge(iE) && _isValidFace(iF))) {
    return false;
  }

  int j = 0;
  bool res = false;
  int iF_j = getEdgeFace(iE, j);
  while (iF_j != -1) {
    res |= iF_j == iF;
    j++;
    iF_j = getEdgeFace(iE, j);
  }
  return res;
}



// classification of vertices

bool PolygonMesh::isBoundaryVertex(const int iV) const {
  int nV = getNumberOfVertices();
  return (0 <= iV && iV < nV) ? _isBoundaryVertex[iV] : false;
}

bool PolygonMesh::isSingularVertex(const int iV) const {
  int nV = getNumberOfVertices();
  return (0 <= iV && iV < nV && _nPartsVertex[iV] > 1);
}

// properties of the whole mesh

bool PolygonMesh::isRegular() const {
  int nV = getNumberOfVertices();
  int nE = getNumberOfEdges();

  bool res = true;

  for (int iE = 0; iE < nE; iE++) {
    res &= !isSingularEdge(iE);
  }

  for (int iV = 0; iV < nV; iV++) {
    res &= !isSingularVertex(iV);
  }
  return res;
}

bool PolygonMesh::hasBoundary() const {
  int nE = getNumberOfEdges();
  bool res = false;

  for (int iE = 0; iE < nE; iE++) {
    res |= isBoundaryEdge(iE);
  }

  return res;
}

//////////////////////////////////////////////////////////////////////
// CONNECTED COMPONENTS

// connected components of the primal graph
// - returns number of connected components nCC
// - fills the faceLabel array with connected component number iCC
// - for each face; 0<=iCC<nCC
int PolygonMesh::computeConnectedComponentsPrimal(vector<int>& faceLabel) const {
  int nCCprimal = 0;
  vector<int> vToCC(getNumberOfVertices(), -1);
  faceLabel.clear();
  faceLabel.resize(getNumberOfFaces(), -1);

  // Create a partition of the vertices using the edges of the primal graph
  Partition partition(getNumberOfVertices());

  for (int iE = 0; iE < getNumberOfEdges(); iE++) {
    int iV0 = getVertex0(iE);
    int iV1 = getVertex1(iE);
    partition.join(iV0, iV1);
  }

  // Fulfill vToCC

  for (int iV = 0; iV < getNumberOfVertices(); iV++) {
    // Get the vertex representative of the partition
    int root = partition.find(iV);
    // Assign a connected component number to the representative
    // if it doesn't have one already
    if (vToCC[root] == -1) {
      vToCC[root] = nCCprimal;
      nCCprimal++;
    }

    // Assign the cc number of the representative root corner
    // to the vertex
    vToCC[iV] = vToCC[root];
  }

  // Fulfill faceLabel
  for (int iC = 0; iC < getNumberOfCorners(); iC++) {
    if (_isValidCorner(iC)) {
      int iF = getFace(iC);
      int iV = _coordIndex[iC];
      faceLabel[iF] = vToCC[iV];
    }
  }


  return nCCprimal;
}

// connected components of the dual graph
// - returns number of connected components nCC
// - fills the faceLabel array with connected component number iCC
// - for each face; 0<=iCC<nCC
int PolygonMesh::computeConnectedComponentsDual(vector<int>& faceLabel) const {
  int nCCdual = 0;
  vector<int> vToCC(getNumberOfFaces(), -1);
  faceLabel.clear();
  faceLabel.resize(getNumberOfFaces(), -1);



  // Create a partition of the vertices using the edges of the dual graph
  Partition partition(getNumberOfFaces());

  // Vertices of the dual graph are faces of the primal graph
  for (int iE = 0; iE < _dualGraph.getNumberOfEdges(); iE++) {
    int iV0 = _dualGraph.getVertex0(iE);
    int iV1 = _dualGraph.getVertex1(iE);
    partition.join(iV0, iV1);
  }

  // Fulfill vToCC
  for (int iF = 0; iF < getNumberOfFaces(); iF++) {
    // Get the vertex representative of the partition
    int root = partition.find(iF);
    // Assign a connected component number to the representative
    // if it doesn't have one already
    if (faceLabel[root] == -1) {
      faceLabel[root] = nCCdual;
      nCCdual++;
    }

    // Assign the connected component number of the representative root corner
    // to the face
    faceLabel[iF] = faceLabel[root];
  }


  return nCCdual;
}

// ORIENTATION

// determines if the mesh is oriented
// - a mesh is oriented if all regular edges are consistently oriented
// - by definition, a mesh with singular edges is not oriented
// - note that isolated singular vertices play no role in this
//   definition (since cuting through them does not affect
//   orientation)
bool PolygonMesh::isOriented() const {
  // If there are singular edges, the mesh is not oriented
  if (hasSingularEdges()) return false;

  // Check all regular edges for consistent orientation
  for (int iE = 0; iE < getNumberOfEdges(); iE++) {

    if (isRegularEdge(iE)) {
      int iC0 = getEdgeHalfEdge(iE, 0);

      // If one regular edge is not consistently oriented,
      // then the mesh is not oriented
      if (HalfEdges::isOriented(iC0) == false) {
        return false;
      }
    }
  }

  return true;
}

bool PolygonMesh::_isOrientableInvertingFaces(vector<bool>& invertFace) const {
  if (hasSingularEdges()) return false;
  if (isOriented()) return true;
  // The mesh is not oriented but only has regular and boundary edges


  int nC = getNumberOfCorners();
  int nF = getNumberOfFaces();
  vector<bool> faceWasVisited(nF, false);

  invertFace.clear();
  invertFace.resize(nF, false);

  vector<int> cornerStack;
  int iC0 = 0;
  int faceRoot = _face[iC0];
  int nVisitedFaces = 0;

  // I don't know if this is the most efficient way to to get the
  // first corner of any face
  Faces faces(getNumberOfVertices(), _coordIndex);


  while (nVisitedFaces < nF) {

    // Find the next unvisited face
    // If there is one single connected component, this will
    // be executed only once, using the first face as root
    while (faceWasVisited[faceRoot]) {
      faceRoot++;
    }

    // Mark the root face as visited
    faceWasVisited[faceRoot] = true;
    nVisitedFaces++;

    int iC_j = faces.getFaceFirstCorner(faceRoot);

    // Traverse all the corners and add them to the stack if they are regular half edges
    _addAllHalfEdgesTo(iC_j, cornerStack);

    // Check all half edges in the stack of the current connected component
    while (!cornerStack.empty()) {
      int iC = cornerStack.back();
      cornerStack.pop_back();
      int iF = getFace(iC);
      int iC_twin = getTwin(iC);
      int iF_twin = getFace(iC_twin);


      if (HalfEdges::isOriented(iC)) {
        // If the half edge is properly oriented, then it depends on whether the current face iF
        // has to be inverted or not
        if (invertFace[iF]) {
          // I should invert the twin face if it is not oriented already
          if (faceWasVisited[iF_twin] && invertFace[iF_twin]) {
            // If the face was visited and has to be inverted, then it's okey
          }

          if (faceWasVisited[iF_twin] && !invertFace[iF_twin]) {
            // Otherwise, the mesh is not orientable
            return false;
          }

          else {
            // Otherwise, the twin face has to be inverted
            faceWasVisited[iF_twin] = true;
            nVisitedFaces++;
            invertFace[iF_twin] = true;

            _addAllHalfEdgesTo(iC_twin, cornerStack);
          }
        }
        else {
          // If the current face iF should not be inverted, then the twin face should not be inverted
          if (faceWasVisited[iF_twin]) {
            if (invertFace[iF_twin]) {
              // If it was visited and has to be inverted, then the mesh is not orientable
              return false;
            }
            // Otherwise, it's okey
          }

          else {
            // Otherwise, the twin face should not be inverted
            faceWasVisited[iF_twin] = true;
            nVisitedFaces++;
            invertFace[iF_twin] = false;

            _addAllHalfEdgesTo(iC_twin, cornerStack);
          }
        }
      }

      else {
        // If the half edge is not properly oriented, then it depends on whether the current face iF
        // has to be inverted or not
        if (invertFace[iF]) {
          // then, the twin face should not be inverted
          if (faceWasVisited[iF_twin]) {
            // If it was visited and has to be inverted, then the mesh is not orientable
            if (invertFace[iF_twin]) {
              return false;
            }
            // If it wasn't inverted, then it's okey
          }

          // Otherwise, the twin face should not be inverted
          else {
            faceWasVisited[iF_twin] = true;
            nVisitedFaces++;
            invertFace[iF_twin] = false;  // Just to be explicit

            _addAllHalfEdgesTo(iC_twin, cornerStack);
          }
        }

        // If the current face iF is not inverted, then the other face should be inverted
        else {
          if (faceWasVisited[iF_twin]) {
            // If it was visited and has to be inverted, then it's okey
            if (invertFace[iF_twin]) {
              // it's okey
            }
            // Otherwise, the mesh is not orientable
            else {
              return false;
            }
          }

          // Otherwise, the twin face should be inverted
          else {
            faceWasVisited[iF_twin] = true;
            nVisitedFaces++;
            invertFace[iF_twin] = true;

            _addAllHalfEdgesTo(iC_twin, cornerStack);
          }
        }
      }
    }
  }

  return true;
}

// Auxiliary method to add all half edges of a face to the corner stack
void PolygonMesh::_addAllHalfEdgesTo(int& iC_j, std::vector<int>& cornerStack) const {
  for (int j = 0; j < getFaceSize(iC_j) - 1; j++) {
    iC_j = getDst(iC_j);
    if (getTwin(iC_j) != -1) {
      cornerStack.push_back(iC_j);
    }
  }
}

// determines if the mesh is orientable
// - a mesh is orientable if a choice of orientation can be made for
//   each face so that, after the face orientation changes are made,
//   the resulting mesh is oriented
// - by definition, a mesh with singular edges is not orientable
// - note that isolated singular vertices play no role in this
//   definition (since cuting through them does not affect
//   orientation)
bool PolygonMesh::isOrientable() const {

  vector<bool> invertFace;

  return _isOrientableInvertingFaces(invertFace);
}

// orient
// - implementation requires a dual graph traversal algorithm
// - the number of connected components nCC of the dual graph are
//   determined as a by product
// - if multiple orientations are posible choose one
// - fills the ccIndex array, of size equal to the number of faces
//   nF, with the connected component number iCC assigned to each face
// - the values stored in the ccIndex shold be in the range 0<=iCC<nCC
// - fills the output invert_face, with values required to produce
//   an oriented mesh
// - the size of output invert_face array should equal to the number
//   of faces
// - returns the number of connected components nCC if successful,
//   and 0 if the mesh is not orientable
// - if not successful, the output arrays should be empty as well
int PolygonMesh::orient(vector<int>& ccIndex, vector<bool>& invertFace) {
  ccIndex.clear();

  int nCC = computeConnectedComponentsDual(ccIndex);
  int nC = getNumberOfCorners();
  int nF = getNumberOfFaces();

  invertFace.clear();
  invertFace.resize(nF, false);
  // If there are singular edges, the mesh is not orientable
  if (hasSingularEdges()) return 0;
  // If the mesh is already oriented, nothing to do
  if (isOriented()) return true;

  _isOrientableInvertingFaces(invertFace);

  // Invert invertFace, I do this because this orientation works better for the
  // he-test-xx.wrl models provided
  for (int i = 0; i < invertFace.size(); i++) {
    invertFace[i] = !invertFace[i];
  }


  return nCC;
}

//////////////////////////////////////////////////////////////////////
// MANIFOLD

// determine how many isolated vertices the mesh has
int PolygonMesh::numberOfIsolatedVertices() {
  vector<int> isolated_vertex;
  getIsolatedVertices(isolated_vertex);
  int nVIsolated = static_cast<int>(isolated_vertex.size());


  return nVIsolated;
}

// get array of isolated vertex indices
void PolygonMesh::getIsolatedVertices(vector<int>& isolated_vertex) {
  isolated_vertex.clear();

  int nV = getNumberOfVertices();

  for (int iV = 0; iV < nV; iV++) {

    // Even if this structure tells which vertices are singular, it can be
    // used to know which vertices are isolated too, since isolated vertices
    // have 0 parts
    if (_nPartsVertex[iV] == 0) {
      isolated_vertex.push_back(iV);
    }
  }

}



// remove isolated vertices
// - the new number of vertices nVout should be <= that the original
//   number of vertices nV
// - fills the coordMap array, of size nVout, with an input vertex
//   input in the range 0<=iV<nV for each output vertex index
//   0<=iVout<nVout
// - fills the output coordIndexOut array with vertex indices in the
//   output range 0<=iVout<nVout
// - the output coordIndexOut should be of the same size as the
//   input coordIndex array
// - returns true if one or more isolated vertices have been removed,
//   and false if no isolated vertices have been found
// - if no isolated vertices are found, the output arrays should be
//   empty as well
bool PolygonMesh::removeIsolatedVertices
(vector<int>& coordMap, vector<int>& coordIndexOut) {
  coordMap.clear();
  coordIndexOut.clear();

  int nV = getNumberOfVertices();

  vector<int> isolatedVertices;
  getIsolatedVertices(isolatedVertices);

  // If there are no isolated vertices, nothing to do
  if (isolatedVertices.empty()) {
    return false;
  }

  int nVOut = nV - isolatedVertices.size();
  int iVCoordIndex = 0;
  int iVIsolated = 0;
  coordMap.resize(nVOut, -1);
  // Build coordMap
  for (int iV = 0; iV < nV; iV++) {
    // Skip isolated vertices
    if (iV == isolatedVertices[iVIsolated]) {
      iVIsolated++;
      continue;
    }
    else {
      coordMap[iVCoordIndex] = iV;
      iVCoordIndex++;
    }
  }

  // Build an array to map old vertex indices to new vertex indices
  vector<int> oldToNewVertexMap(nV, -1);

  // Fill oldToNewVertexMap
  for (int iVOut = 0; iVOut < nVOut; iVOut++) {
    int iV = coordMap[iVOut];
    oldToNewVertexMap[iV] = iVOut;
  }

  // Build coordIndexOut
  int nC = getNumberOfCorners();
  coordIndexOut.resize(nC, -1);
  for (int iC = 0; iC < nC; iC++) {
    int iV = _coordIndex[iC];
    if (iV == -1) {
      // Assign face separators the same as input
      coordIndexOut[iC] = -1;
    }
    else {
      coordIndexOut[iC] = oldToNewVertexMap[iV];
    }
  }

  return true;
}


// cut through singular vertices
// - should only cut through singular vertices which belong to
//   different connected components of the dual graph
// - it should also remove isolated vertices
// - singular edges should not be modified
// - it should work on non-orientable meshes of any kind
// - the new number of vertices nVout should be => that the original
//   number of vertices nV
// - if nVout==nV, should return empty vIndexMap and coordIndexOut
// - otherwise
// - fills the vIndexMap array, of size nVout, with an input vertex
//   input in the range 0<=iV<nV for each output vertex index
//   0<=iVout<nVout
// - fills the output coordIndexOut array with vertex
//   indices in the output range 0<=iVout<nVout
// - the output coordIndexOut should be of the same size as the
//   input coordIndex array
void PolygonMesh::cutThroughSingularVertices
(vector<int>& vIndexMap, vector<int>& coordIndexOut) {
  vIndexMap.clear();
  coordIndexOut.clear();

  vector<int> coordMap;
  vector<int> coordIndexRemovingIsolatedVertices;
  bool removedIsolatedVertices = removeIsolatedVertices(coordMap, coordIndexRemovingIsolatedVertices);

  // If we did not remove isolated vertices, we need to create an identity coordMap
  if (!removedIsolatedVertices) {
    coordIndexRemovingIsolatedVertices = _coordIndex;

    // Make coordMap identity
    int nV = getNumberOfVertices();
    coordMap.resize(nV, -1);
    for (int iV = 0; iV < nV; iV++) {
      coordMap[iV] = iV;
    }
  }

  int nV = getNumberOfVertices();
  int nVOut = nV;
  int nC = coordIndexRemovingIsolatedVertices.size();

  coordIndexOut.resize(nC, -1);


  Partition partitionRegularSingular(nC);

  // Build a partition of the corners by joining pairs of
  // matching corners across regular and singular edges
  // Even if we removed isolated vertices, the number of edges is still the same
  // so we can use getNumberOfEdges() safely
  for (int iE = 0; iE < getNumberOfEdges(); iE++) {
    int nFaces = getNumberOfEdgeHalfEdges(iE);

    // We join corners across regular edges and singular edges
    if (nFaces >= 2) {
      for (int j = 0; j < nFaces; j++) {
        int iCj = getEdgeHalfEdge(iE, j);
        int iCj1 = getDst(iCj);

        int iVj = coordMap[coordIndexRemovingIsolatedVertices[iCj]];
        int iVj1 = coordMap[coordIndexRemovingIsolatedVertices[iCj1]];

        for (int k = j + 1; k < nFaces; k++) {
          int iCk = getEdgeHalfEdge(iE, k);
          int iCk1 = getDst(iCk);

          int iVk = coordMap[coordIndexRemovingIsolatedVertices[iCk]];
          int iVk1 = coordMap[coordIndexRemovingIsolatedVertices[iCk1]];


          // this case considers that the two half edges are oriented
          if (iVj == iVk1) {
            partitionRegularSingular.join(iCj, iCk1);
            partitionRegularSingular.join(iCj1, iCk);
          }
          // Otherwise, this case considers that the two half edges are opositely oriented
          else {
            partitionRegularSingular.join(iCj, iCk);
            partitionRegularSingular.join(iCj1, iCk1);
          }
        }
      }
    }
  }

  // For each vertex, count how many connected components it belongs to
  vector<int> _nCCVertex;
  _nCCVertex.resize(nV, 0);
  vector<int> rootCornerToNewVertex(nC, -1);

  // Create a partition of the corners by joining pairs of
  // matching corners across regular and singular edges
  // Even if we removed isolated vertices, the number of edges is still the same
  // so we can use getNumberOfEdges() safely

  vIndexMap.resize(nV, -1);
  bool hasSingularVertices = false;
  
  for (int iC = 0; iC < nC; iC++) {
    // For each root corner of the partition
    //  if iC is a root corner, process it
    if (partitionRegularSingular.find(iC) == iC) {
      // If the corner is a face separator, skip it
      if (coordIndexRemovingIsolatedVertices[iC] == -1)
        continue;

      // Get the vertex index of the corner
      int iV = coordMap[coordIndexRemovingIsolatedVertices[iC]];

      // Increment the number of connected components the vertex belongs to
      _nCCVertex[iV]++;
      // If the vertex is in more than one partition, then it means that it is singular
      // but belongs to different components in the dual graph
      // So, I assign a new vertex to the new CC and map it in vIndexMap
      if (_nCCVertex[iV] > 1) {
        rootCornerToNewVertex[iC] = nVOut;
        nVOut++;
        vIndexMap.push_back(iV);
        hasSingularVertices = true;
      }
      // Otherwise, it is not neccesary to duplicate the vertex, because for now it is in one single CC
      else {
        rootCornerToNewVertex[iC] = iV;
        vIndexMap[iV] = iV;
      }
    }
  }

  if (!hasSingularVertices && !removedIsolatedVertices) {
    vIndexMap.clear();
    coordIndexOut.clear();

    return;
  }
  else if (!hasSingularVertices && removedIsolatedVertices) {
    coordIndexOut = coordIndexRemovingIsolatedVertices;
    vIndexMap = coordMap;
    return;
  }


  int iCOut = 0;
  nVOut = coordMap.size();
  // Assign the new vertices into the coordIndexOut
  for (int iC = 0; iC < nC; iC++) {

    if (coordIndexRemovingIsolatedVertices[iC] == -1) {
      coordIndexOut[iCOut] = -1; // face separator
      iCOut++;
      continue;
    }

    int iCroot = partitionRegularSingular.find(iC);

    coordIndexOut[iCOut] = rootCornerToNewVertex[iCroot];

    iCOut++;
  }

  // It is not necessary to repeat the computation of dualGraph because the conections
  // between faces have not changed (since we only cut through singular vertices)


}

// convert to manifold
// - removes isolated vertices, cuts through singular vertices and
//   through singular edges
// - the new number of vertices nVout should be => that the original
//   number of vertices nV
// - if nVout==nV, should return empty vIndexMap and coordIndexOut
// - otherwise
// - fills the vIndexMap array, of size nVout, with an input vertex
//   input in the range 0<=iV<nV for each output vertex index
//   0<=iVout<nVout
// - fills the output coordIndexOut array with vertex
//   indices in the output range 0<=iVout<nVout
// - the output coordIndexOut should be of the same size as the
//   input coordIndex array

void PolygonMesh::convertToManifold
(vector<int>& vIndexMap, vector<int>& coordIndexOut) {
  bool success = false;
  vIndexMap.clear();
  coordIndexOut.clear();

  vector<int> coordMap;
  vector<int> coordIndexRemovingIsolatedVertices;
  bool removedIsolatedVertices = removeIsolatedVertices(coordMap, coordIndexRemovingIsolatedVertices);

  // If we did not remove isolated vertices, we need to create an identity coordMap
  if (!removedIsolatedVertices) {
    coordIndexRemovingIsolatedVertices = _coordIndex;

    // Make coordMap identity
    int nV = getNumberOfVertices();
    coordMap.resize(nV, -1);
    for (int iV = 0; iV < nV; iV++) {
      coordMap[iV] = iV;
    }
  }

  int nV = getNumberOfVertices();
  int nVOut = nV;
  int nC = coordIndexRemovingIsolatedVertices.size();

  coordIndexOut.resize(nC, -1);

  Partition partition(nC);

  int nE = getNumberOfEdges();
  for (int iE = 0; iE < nE; iE++) {
    int nFaces = getNumberOfEdgeHalfEdges(iE);

    // If there is exactly two faces incident to the edge
    // then the edge is regular and we can join the corners
    // accross the edge
    if (nFaces == 2) {

      int iC00 = getEdgeHalfEdge(iE, 0);
      int iC01 = getDst(iC00);

      int iC10 = getEdgeHalfEdge(iE, 1);
      int iC11 = getDst(iC10);

      if (HalfEdges::isOriented(iC00)) {
        // Consistently oriented
        partition.join(iC00, iC11);
        partition.join(iC01, iC10);
      }
      else {
        // Oposite orientation
        partition.join(iC00, iC10);
        partition.join(iC01, iC11);
      }
    }
  }



  vector<int> rootCornerToNewVertex(nC, -1);
  vector<int> local_nPartsVertex(nV, 0);

  // Create a partition of the corners by joining pairs of
  // matching corners across regular and singular edges
  // Even if we removed isolated vertices, the number of edges is still the same
  // so we can use getNumberOfEdges() safely

  vIndexMap.resize(nV, -1);

  bool hasSingularVertices = false;
  for (int iC = 0; iC < nC; iC++) {
    if (partition.find(iC) == iC) {
      if (coordIndexRemovingIsolatedVertices[iC] == -1)
        continue;
      int iV = coordMap[coordIndexRemovingIsolatedVertices[iC]];

      local_nPartsVertex[iV]++;
      // If the vertex is in more than one partition, then it mean that it is singular
      // but belongs to different components in the dual graph
      // So, I assign a new vertex to the new CC and map it in vIndexMap
      if (local_nPartsVertex[iV] > 1) {
        rootCornerToNewVertex[iC] = nVOut;
        nVOut++;
        vIndexMap.push_back(iV);
        hasSingularVertices = true;
      }
      // Otherwise, it is not neccesary to duplicate the vertex, because for now it is in one single CC
      else {
        rootCornerToNewVertex[iC] = iV;
        vIndexMap[iV] = iV;
      }
    }
  }

  if (!hasSingularVertices && !removedIsolatedVertices) {
    vIndexMap.clear();
    coordIndexOut.clear();

    return;
  }
  else if (!hasSingularVertices && removedIsolatedVertices) {
    coordIndexOut = coordIndexRemovingIsolatedVertices;
    vIndexMap = coordMap;
    return;
  }


  int iCOut = 0;
  nVOut = coordMap.size();
  // Assign the new vertices into the coordIndexOut
  for (int iC = 0; iC < nC; iC++) {

    if (coordIndexRemovingIsolatedVertices[iC] == -1) {
      coordIndexOut[iCOut] = -1; // face separator
      iCOut++;
      continue;
    }

    int iCroot = partition.find(iC);

    coordIndexOut[iCOut] = rootCornerToNewVertex[iCroot];

    iCOut++;
  }

  _hasSingularEdges = false;

  _create_nPartsVertex(coordMap.size(), nC, partition);
  _createDualGraph();
  _createIsBoundaryVertex(coordMap.size(), nE);

}

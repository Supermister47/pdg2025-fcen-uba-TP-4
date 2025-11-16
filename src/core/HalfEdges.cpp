//------------------------------------------------------------------------
//  Copyright (C) Gabriel Taubin
//  Time-stamp: <2025-08-05 23:13:01 taubin>
//------------------------------------------------------------------------
//
// HalfEdges.cpp (Assignment 3)
//
// Written by: <Your Name>
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


#include <math.h>
#include "HalfEdges.hpp"
#include "Graph.hpp"
#include <io/StrException.hpp>

// #include "Graph.hpp"

// TODO Mon Mar 6 2023
// - merge your code from Assignment 2

HalfEdges::HalfEdges(const int nVertices, const vector<int>& coordIndex) : Edges(nVertices), // a graph with no edges is created here
_coordIndex(coordIndex),
_twin(),
_face(),
_firstCornerEdge(),
_cornerEdge() {

  int nV = nVertices;
  int nC = static_cast<int>(_coordIndex.size()); // number of corners

  // 0) just to be safe, verify that for each corner iC that
  //    -1<=iV && iV<nV, where iV=coordIndex[iC]
  //
  //    if you find a violation, you can increment and the variable
  //    nV, and then use the method Edges::_reset() to adjust the
  //    number of vertices of the graph, if necessary; or you can
  //    abort throwing an StrException

  int iV;

  for (int iC = 0; iC < nC; iC++) {
    iV = _coordIndex[iC];

    if (!(-1 <= iV && iV < nV)) {
      throw StrException("There is at least one vertex which is not a valid index.");
    }
  }

  printf("All vertices ok\n");

  // 1) create an empty vector<int> to count the number of incident
  //    faces per edge; size is not known at this point because the
  //    edges have not been created yet
  vector<int> nFacesEdge;

  // 2) insert all the edges in the graph; at the same time initialize
  //    the _twin array so that all the half edges are boundary, count
  //    the number of incident faces per edge, fill the _face
  //    array, and count the number of faces incident to each edge
  _face.resize(nC, -1);

  int iV0, iV1, iF, iE, iC, iC0, iC1;
  for (iF = iC0 = iC1 = 0; iC1 < nC; iC1++) {
    if (_coordIndex[iC1] >= 0)
      continue;
    // face iF comprises corners iC0<=iC<iC1
    // - each corner in this range corresponds to one half edge
    // - find the next corner within the face
    // - get the two vertex indices and insert an edge in the graph if
    //   not already there
    for (iC = iC0; iC < iC1 - 1; iC++) {
      iV0 = _coordIndex[iC];
      iV1 = _coordIndex[iC + 1];

      iE = _insertEdge(iV0, iV1);

      // If the edge is new, add it to nFacesEdge with 1 incident face
      if (iE >= nFacesEdge.size()) {
        nFacesEdge.push_back(1);
      }
      // Otherwise, increase by one the amount of incident faces
      else {
        nFacesEdge[iE] = nFacesEdge[iE] + 1;
      }

      // assign face index to corner
      _face[iC] = iF;
    }

    // Add the edge between the last corner and the first
    // It repeats the code in the loop above, but well...
    _face[iC] = iF;
    iE = _insertEdge(_coordIndex[iC1 - 1], _coordIndex[iC0]);

    // If the edge is new, add it to nFacesEdge with 1 incident face
    if (iE >= nFacesEdge.size()) {
      nFacesEdge.push_back(1);
    }
    // Otherwise, increase by one the amount of incident faces
    else {
      nFacesEdge[iE] = nFacesEdge[iE] + 1;
    }

    // increment variables to continue processing next face
    iC0 = iC1 + 1;
    iF++;
  }

  _nF = iF;
  int nE = getNumberOfEdges();

  // 3) create an array to hold the first twin corner for each edge
  // - the size of this array should be equal to the number of edges
  // - initialize it with -1's
  vector<int> twinCorner;
  twinCorner.resize(nE, -1);
  _twin.resize(nC, -1);

  // 4) fill the _twin array
  // - visit all the half edges using a loop similar to the one used in step 2)
  // - for each half edge iC, get the src and dst vertex indices, and
  //   from them the index iE of the corresponding edge
  // - if twinCorner[iE]<1 save iC in twinCorner[iE]
  //   otherwise save the value stored in twinCorner[iE] in _twin[iC]
  //   and iC in _twin[_twin[iC]]

  for (iF = iC0 = iC1 = 0; iC1 < nC; iC1++) {
    if (_coordIndex[iC1] >= 0)
      continue;

    for (iC = iC0; iC < iC1 - 1; iC++) {
      iV0 = _coordIndex[iC];
      iV1 = _coordIndex[iC + 1];

      iE = getEdge(iV0, iV1);

      if (iE == -1) // This should not happen, but anyway
        throw StrException("Primal graph is not properly initialized");

      // If it is the first time we see this edge, save the corner
      if (twinCorner[iE] == -1) {
        twinCorner[iE] = iC;
      }
      // Otherwise, make them twins
      else {
        _twin[iC] = twinCorner[iE];
        _twin[twinCorner[iE]] = iC;
      }
    }

    // Edge between last corner and first corner
    // It repeats the code in the loop above AGAIN, but well...
    iV0 = _coordIndex[iC1 - 1];
    iV1 = _coordIndex[iC0];

    iE = getEdge(iV0, iV1);

    if (iE == -1) // This should not happen, but anyway
      throw StrException("Primal graph is not properly initialized");

    // If it is the first time we see this edge, save the corner
    if (twinCorner[iE] == -1) {
      twinCorner[iE] = iC1 - 1;
    }
    // Otherwise, make them twins
    else {
      _twin[iC] = twinCorner[iE];
      _twin[twinCorner[iE]] = iC1 - 1;
    }

    // increment variables to continue processing next face
    iC0 = iC1 + 1;
    iF++;
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

  // a decision has to be made about inconsistently oriented half
  // edges incident to the same edge, as well as how to deal with
  // singular edges; for the moment let's assume that the mesh does
  // not have singular edges, but inconsistently oriented half edges
  // incident to the same edge are made twins (i.e. we do not have to
  // check for orientation here); later on we may want to modify this
  // class to have the option to do one thing or the other, and
  // methods to indicate which case we have.

  // get everything up to here implemented, debugged, and commited
  // before continuing

  // 5) initialize the array of arrays representing the half-edge to
  //    edge incident relationships
  //    _firstCornerEdge, and _cornerEdge
  //    - the size of _firstCornerEdge should be equal to nE+1
  //    - the size of _cornerEdge should be equal to the number of valid corners
  //      nC-nF
  //    - set boundaries
  //      _firstCornerEdge[0]=0
  //      _firstCornerEdge[iE+1] = _firstCornerEdge[iE]+nFacesEdge[iE] (1<=iE<nE)
  _firstCornerEdge.resize(nE + 1, -1);
  _cornerEdge.resize(nC - _nF, -1);

  _firstCornerEdge[0] = 0;
  for (iE = 0; iE < nE; iE++) {
    _firstCornerEdge[iE + 1] = _firstCornerEdge[iE] + nFacesEdge[iE];
  }

  // 6) fill the array of arrays - the indices of corners incident to
  //    edge iE (1 if boundary, 2 if regular, >2 if singular) should
  //    be stored consecutively in _cornerEdge starting at the
  //    location _firstCornerEdge[iE]

  for (iF = iC0 = iC1 = 0; iC1 < nC; iC1++) {
    if (_coordIndex[iC1] >= 0)
      continue;

    for (iC = iC0; iC < iC1 - 1; iC++) {
      iV0 = _coordIndex[iC];
      iV1 = _coordIndex[iC + 1];

      iE = getEdge(iV0, iV1);

      int iE_i = 0;

      // Look for the next available position
      while (_cornerEdge[_firstCornerEdge[iE] + iE_i] != -1) {
        iE_i++;
      }

      // Save the corner index
      _cornerEdge[_firstCornerEdge[iE] + iE_i] = iC;

      // I make use of this cycle so i can see if there are regular, singular or boundary edges
      _classifyEdge(iE);
    }

    // Edge between last corner and first corner
    // It repeats the code in the loop above AGAIN x2, but well...

    iV0 = _coordIndex[iC1 - 1];
    iV1 = _coordIndex[iC0];

    iE = getEdge(iV0, iV1);

    int iE_i = 0;
    while (_cornerEdge[_firstCornerEdge[iE] + iE_i] != -1) {
      iE_i++;
    }
    _cornerEdge[_firstCornerEdge[iE] + iE_i] = iC1 - 1;

    iC0 = iC1 + 1;
    iF++;
  }
}

void HalfEdges::_classifyEdge(int iE) {
  if (isBoundaryEdge(iE)) {
    _hasBoundaryEdges = true;
  }
  else if (isRegularEdge(iE)) {
    _hasRegularEdges = true;
  }
  else if (isSingularEdge(iE)) {
    _hasSingularEdges = true;
  }
}

int HalfEdges::getNumberOfCorners() const {
  return static_cast<int>(_coordIndex.size());
}

// in all subsequent methods check that the arguments are valid, and
// return -1 if any argument is out of range

bool HalfEdges::_isValidCorner(const int iC) const {
  return 0 <= iC && iC < getNumberOfCorners() && _coordIndex[iC] != -1;
}

bool HalfEdges::_isValidEdge(const int iE) const {
  return 0 <= iE && iE < getNumberOfEdges();
}

// half-edge method srcVertex()
int HalfEdges::getFace(const int iC) const {
  if (!_isValidCorner(iC)) {
    return -1;
  }

  return _face[iC];
}

// half-edge method srcVertex()
int HalfEdges::getSrc(const int iC) const {
  if (!_isValidCorner(iC)) {
    return -1;
  }

  return iC;
}

// half-edge method dstVertex()
int HalfEdges::getDst(const int iC) const {

  return getSrc(getNext(iC));
}

// half-edge method next()
int HalfEdges::getNext(const int iC) const {
  // if iC is the last corner of its face, use the face size
  // stored in _twin[iC+1] to locate the first corner of the face

  if (!_isValidCorner(iC)) {
    return -1;
  }

  // Case if it is the last corner
  if (iC == getNumberOfCorners() - 1 || _coordIndex[iC + 1] == -1) {
    // Search backwards for the first corner of the face
    for (int iC_search = iC - 1; iC_search >= 0; iC_search--) {
      if (_coordIndex[iC_search] == -1) {
        return iC_search + 1;
      }
    }
    // If we reach here, the first corner is at index 0
    return 0;
  }
  else {
    return iC + 1;
  }
}

// half-edge method prev()
int HalfEdges::getPrev(const int iC) const {

  // if iC is the first corner of its face, since the face size is
  // stored at the end of the face in the _twin array, you will have
  // to search forward for the last corner of the face; you can use
  // the fact that all the faces have at least 3 corners to start the
  // search for the face separator at iC+3

  if (!_isValidCorner(iC)) {
    return -1;
  }



  int iF = _face[iC];
  // Check if the corner is the first corner of a face
  if (iC == 0 || _coordIndex[iC - 1] == -1) {
    // Search forwards for the last corner of the face
    for (int iC_search = iC + 3; iC_search < getNumberOfCorners(); iC_search++) {
      if (_coordIndex[iC_search] == -1) {
        return iC_search - 1;
      }
    }
    // If we reach here, the last corner is at the last index
    return getNumberOfCorners() - 1;
  }
  else {
    return iC - 1;
  }
}

int HalfEdges::getTwin(const int iC) const {
  if (!_isValidCorner(iC)) {
    return -1;
  }

  return _twin[iC];
}

// represent the half edge as an array of lists, with one list
// associated with each edge

int HalfEdges::getNumberOfEdgeHalfEdges(const int iE) const {
  if (!_isValidEdge(iE)) {
    return 0;
  }

  // Compute the number of half-edges for this edge
  return _firstCornerEdge[iE + 1] - _firstCornerEdge[iE];
}

int HalfEdges::getEdgeHalfEdge(const int iE, const int j) const {
  if (!(_isValidEdge(iE) && 0 <= j && j < getNumberOfEdgeHalfEdges(iE))) {
    return -1;
  }

  // Return the j-th half-edge for this edge
  return _cornerEdge[_firstCornerEdge[iE] + j];
}
// TODO Mon Mar 6 2023
// - new functions to implement

bool HalfEdges::isOriented(const int iC) const {
  // iC     : iV00->iV01
  // iCtwin : iV10->iV11

  // if (twin half edges are consistently oriented)
  /* \                  / */
  /*  \ iV01 <-- iV00  /  */
  /*   X ------------ X   */
  /*  / iV10 --> iV11  \  */
  /* /                  \ */
  // return true;
  
  // else if (twin half edges are not consistently oriented)
  /* \                  / */
  /*  \ iV01 <-- iV00  /  */
  /*   X ------------ X   */
  /*  / iV11 <-- iV10  \  */
  /* /                  \ */
  // return false;

  if (!_isValidCorner(iC)) {
    return false;
  }
  int iCtwin = getTwin(iC);
  if (iCtwin == -1) {
    return true;
  }

  // Get its vertices
  int iV00 = _coordIndex[getSrc(iC)];
  int iV01 = _coordIndex[getDst(iC)];
  int iV10 = _coordIndex[getSrc(iCtwin)];
  int iV11 = _coordIndex[getDst(iCtwin)];

  return iV00 == iV11 && iV01 == iV10;
}

// half-edge method getFaceSize()
int HalfEdges::getFaceSize(const int iC) const {
  if (!_isValidCorner(iC)) {
    return -1;
  }
  int faceSize = 0;
  for (int iCAnt=iC; _isValidCorner(iCAnt) && _coordIndex[iCAnt] != -1; iCAnt--) {
    faceSize++;
  }

  for (int iCSig=iC+1; _isValidCorner(iCSig) && _coordIndex[iCSig] != -1; iCSig++) {
    faceSize++;
  }

  return faceSize;
}
  
int HalfEdges::getNumberOfFacesEdge(const int iE) const {
  if (!_isValidEdge(iE)) {
    return -1;
  }

  int firstIndex = _firstCornerEdge[iE];
  int lastIndex = _firstCornerEdge[iE + 1];
  int nFaces = lastIndex - firstIndex;
  return nFaces;
}

bool HalfEdges::isBoundaryEdge(const int iE) const {
  if (!_isValidEdge(iE)) {
    return false;
  }

  return getNumberOfFacesEdge(iE) == 1;
}


bool HalfEdges::isRegularEdge(const int iE) const {
  if (!_isValidEdge(iE)) {
    return false;
  }

  return getNumberOfFacesEdge(iE) == 2;
}

bool HalfEdges::isSingularEdge(const int iE) const {
  if (!_isValidEdge(iE)) {
    return false;
  }

  return getNumberOfFacesEdge(iE) > 2;
}

// HINTS

// - it is best to determine the return values of these methods in the
//   constructor and save them in private variables

bool HalfEdges::hasBoundaryEdges() const {
  return _hasBoundaryEdges;
}

bool HalfEdges::hasRegularEdges() const {
  return _hasRegularEdges;
}

bool HalfEdges::hasSingularEdges() const {
  return _hasSingularEdges;
}

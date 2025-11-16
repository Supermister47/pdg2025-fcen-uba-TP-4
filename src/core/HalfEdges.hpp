//------------------------------------------------------------------------
//  Copyright (C) Gabriel Taubin
//  Time-stamp: <2025-08-05 23:13:02 taubin>
//------------------------------------------------------------------------
//
// HalfEdges.hpp (Assignment 3)
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

#ifndef _HALF_EDGES_HPP_
#define _HALF_EDGES_HPP_

#include <vector>
#include "Edges.hpp"

using namespace std;

class HalfEdges : public Edges {

public:

  // methods inherited from Edges
  //
  // int     getNumberOfVertices()                     const;
  // int     getNumberOfEdges()                        const;
  // int     getEdge(const int iV0, const int iV1)     const;
  // int     getVertex0(const int iE)                  const;
  // int     getVertex1(const int iE)                  const;

  // constructor performs most of the work

    HalfEdges(const int nV, const vector<int>& coordIndex);

    void _classifyEdge(int iE);

  // returns the number of elements of the coordIndex array

  int     getNumberOfCorners() const;

  // half-edges are in one-to-one correspondence with the corners of a
  // mesh, i.e., with the indices of the coordIndex array which do not
  // correspond to face separators; if the corner index iC is out of
  // the range 0<=iC<coordIndex.size(), or coordIndex[iC]<0, these
  // methods return -1; these two methods return vertex indices;

  int     getSrc(const int iC) const;
  int     getDst(const int iC) const;

  // the mesh faces define loops of half edges; these two methods can
  // be used to move back and forth along these loops;

  int     getNext(const int iC) const;
  int     getPrev(const int iC) const;

  // a regular edge of a mesh has exactly two incident half-edges; if
  // the half-edge associated with corner iC corresponds to a regular
  // edge of the mesh, this methods returns the other half edge;
  // otherwie it returns -1

  int     getTwin(const int iC) const;

  // if the edge index iE is in range, this method returns the number
  // of half edges incident to the given edge; otherwise it returns 0
                                   
  int     getNumberOfEdgeHalfEdges(const int iE) const;

  // if the edge index iE is in range, and
  // 0<=j<getNumberOfEdgeHalfEdges(iE), this method returns the j-th
  // corner corresponding to a half edge incident to the given edge

  int     getEdgeHalfEdge(const int iE, const int j) const;

  // TODO Mon Mar 6 2023
  // - new functions 

  // returns the index of the face containing the half edge
  // corresponding to the corner index iC; if the corner index is out
  // of range, or it corresponds to a face separator, this method
  // returns -1;

  int     getFace(const int iC) const;

  // returns the size of the face containing the half edge
  // corresponding to the corner index iC; if the corner index is out
  // of range, or it corresponds to a face separator, this method
  // returns -1;
  
  int     getFaceSize(const int iC) const;

  // returns true if the corresponding edge is regular and the twin is
  // consistently oriented; otherwise it returns false

  bool    isOriented(const int iC) const;

  // returns the number of faces incident to an edge if the edge index
  // is in the valid range, i.e., if 0<=iE<getNumberOfEdges();
  // otherwise it returns -1

  int     getNumberOfFacesEdge(const int iE) const;

  // classification of edges

  bool    hasBoundaryEdges() const;
  bool    hasRegularEdges() const;
  bool    hasSingularEdges() const;

  bool    isBoundaryEdge(const int iE) const;
  bool    isRegularEdge(const int iE) const;
  bool    isSingularEdge(const int iE) const;
  
protected:

    bool _isValidCorner(const int iC) const;
    bool _isValidEdge(const int iE) const;

    int _nF;

    // reference to the coordIndex passed as argument
    const vector<int>& _coordIndex;

    // - consider these private variables are just a hint
    // - feel free to use different private variables

    // array of twin corners
    vector<int> _twin;

    // mapping from corners to faces
    vector<int> _face;

    // the half-edge to edge incidence relations is represented as an
    // arrray of arrays
    vector<int> _firstCornerEdge;
    vector<int> _cornerEdge;


  bool _hasBoundaryEdges;
  bool _hasRegularEdges;
  bool _hasSingularEdges;

};

#endif /* _HALF_EDGES_HPP_ */

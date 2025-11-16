//------------------------------------------------------------------------
//  Copyright (C) Gabriel Taubin
//  Time-stamp: <2025-08-05 23:13:05 taubin>
//------------------------------------------------------------------------
//
// PolygonMesh.hpp
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

#ifndef _POLYGON_MESH_HPP_
#define _POLYGON_MESH_HPP_

#include <vector>
#include "HalfEdges.hpp"
#include "Graph.hpp"
#include "Faces.hpp"
#include "Partition.hpp"

using namespace std;

class PolygonMesh : public HalfEdges {

public:

  // inherits from Edges
  //
  // void    reset(const int nV);
  // int     getNumberOfVertices()                     const;
  // int     getNumberOfEdges()                        const;
  // int     getEdge(const int iV0, const int iV1)     const;
  // int     getVertex0(const int iE)                  const;
  // int     getVertex1(const int iE)                  const;

  // inherits from HalfEdges
  //
  // int     getNumberOfCorners();
  // int     getFace(const int iC) const;
  // int     getSrc(const int iC) const;
  // int     getDst(const int iC) const;
  // int     getNext(const int iC) const;
  // int     getPrev(const int iC) const;
  // int     getTwin(const int iC) const;
  // int     getNumberOfFaceEdges(const int iE);

    PolygonMesh(const int nV, const vector<int>& coordIndex);

    void _createIsBoundaryVertex(int nV, int nE);

    void _create_nPartsVertex(const int nVertices, int nC, Partition& partition);

    void _createDualGraph();

  // number of -1's in the coordIndex argument

     int     getNumberOfFaces()                        const;

  // number of faces incident to each edge; note that this is equal to
  // the number of half edges incident to each edge

     int     getNumberOfEdgeFaces(const int iE)        const;

  // if the arguments fall within their respective ranges, this method
  // returns the j-th face in the list of faces incident to the edge
  // iE; and it returns -1 if either argument is out of range

     int     getEdgeFace(const int iE, const int j)    const;

  // if the arguments fall within their respective ranges, this method
  // returns returns true if iF is found in the list of faces incident
  // to the edge iE; otherwise it returns false

     bool    isEdgeFace(const int iE, const int iF)    const;

  // a vertex is boundary if and only if it is the end of a boundary
  // edge

     bool    isBoundaryVertex(const int iV)            const;

  // a vertex is internal if and only if it is not a boundary edge

     bool    isInternalVertex(const int iV)            const;

  // a vertex is singular if the number of connected components in the
  // subgraph of the dual graph defined by the subset of faces
  // incident to the vertex is larger than 1; otherwise it is regular

     bool    isSingularVertex(const int iV)            const;

  // a way to determine which vertices are singular and which are
  // regular is to construct a partition of the corners of the mesh,
  // i.e. the elements of the coordIndex array, including the -1
  // separators. For each regular edge there are exactly two
  // half-edges incident to the edge.

  // If the two incident faces are consistently oriented, the
  // neighborhood of the edge will look like this

  //
  /* \                  / */
  /*  \ iC01 <-- iC00  /  */
  /*   X ------------ X   */
  /*  / iC10 --> iC11  \  */
  /* /                  \ */

  // If the two incident faces are not consistently oriented, the
  // neighborhood of the edge will look like this

  //
  /* \                  / */
  /*  \ iC00 --> iC01  /  */
  /*   X ------------ X   */
  /*  / iC10 --> iC11  \  */
  /* /                  \ */

  // In either case the two pairs of corners which are opposite to
  // each other accross the edge jshould be joined. Once all of these
  // join operations have been applied, the corners are partitioned
  // into parts. Some of these parts are singletons corresponding to
  // the face separators. Any other part is composed of corners which
  // point to the same vertex. Counting how many of these parts point
  // to the same vertex determines which vertices are regular and
  // which are singular. If only one of these parts point to a given
  // vertex, then the vertex is regular. If two or more of these parts
  // point to a given vertex, then the vertex is singular. 

  // the polygon mesh is regular if and only if it does not have any
  // singular edges and it does not have any singular vertices

     bool    isRegular()                               const;

  // the polygon mesh has boundary if and only if it has at least 1
  // boundary edge

     bool    hasBoundary()                             const;

  // CONNECTED COMPONENTS

  // connected components of the primal graph
  // - returns number of connected components nCC
  // - fills the faceLabel array with connected component number iCC
  // - for each face; 0<=iCC<nCC
  int computeConnectedComponentsPrimal(vector<int>& faceLabel) const;

  // connected components of the dual graph
  // - returns number of connected components nCC
  // - fills the faceLabel array with connected component number iCC
  // - for each face; 0<=iCC<nCC
  int computeConnectedComponentsDual(vector<int>& faceLabel) const;

  // ORIENTATION
  
  // determines if the mesh is oriented
  // - a mesh is oriented if all regular edges are consistently oriented
  // - by definition, a mesh with singular edges is not oriented
  // - note that isolated singular vertices play no role in this
  //   definition (since cuting through them does not affect
  //   orientation)
  bool isOriented() const;

  // determines if the mesh is orientable
  // - a mesh is orientable if a choice of orientation can be made for
  //   each face so that, after the face orientation changes are made,
  //   the resulting mesh is oriented
  // - by definition, a mesh with singular edges is not orientable
  // - note that isolated singular vertices play no role in this
  //   definition (since cuting through them does not affect
  //   orientation)
  bool isOrientable() const;
  void _addAllHalfEdgesTo(int& iC_j, std::vector<int>& cornerStack) const;
  bool _isOrientableInvertingFaces(vector<bool>& invertFace) const;
 

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
  int orient(vector<int>& ccIndex, vector<bool>& invert_face);

  // MANIFOLD

  // determine how many isolated vertices the mesh has
  int numberOfIsolatedVertices();

  // get array of isolated vertex indices
  void getIsolatedVertices(vector<int>& isolated_vertex);

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
  // - returns true if successful, and false otherwise
  // - if not successful, the output arrays should be empty as well
  bool removeIsolatedVertices
  (vector<int>& coordMap, vector<int>& coordIndexOut);

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
  void cutThroughSingularVertices
  (vector<int>& vIndexMap, vector<int>& coordIndexOut);

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
  void convertToManifold
  (vector<int>& vIndexMap, vector<int>& coordIndexOut);

private:
   bool _isValidFace(const int iF) const;
  vector<int>      _nPartsVertex; // if _nPartsVertex[iV]>1 => vertex is singular 
  vector<bool> _isBoundaryVertex;
  Graph _dualGraph;
  
};

#endif /* _POLYGON_MESH_HPP_ */

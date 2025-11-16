//------------------------------------------------------------------------
//  Copyright (C) Gabriel Taubin
//  Time-stamp: <2025-08-06 09:57:23 taubin>
//------------------------------------------------------------------------
//
// LoaderOff.cpp
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

#include <cstdio>
#include <cstring>
#include "TokenizerFile.hpp"
#include "TokenizerString.hpp"
#include "LoaderOff.hpp"
#include "StrException.hpp"

#include "core/Geometry.hpp"
#include "wrl/Shape.hpp"
#include "wrl/Appearance.hpp"
#include "wrl/Material.hpp"
#include "wrl/IndexedFaceSet.hpp"

// reference
// https://en.wikipedia.org/wiki/OFF_(file_format)

const char* LoaderOff::_ext = "off";

// static
bool LoaderOff::_addFaceNormals = true;

// static
void LoaderOff::setAddFaceNormals(const bool value) {
  _addFaceNormals = value;
}

// static
IndexedFaceSet* LoaderOff::_initializeSceneGraph
(SceneGraph& wrl, const char* filename, const string shapeName) {
  // 0) clear the container
  wrl.clear();
  wrl.setUrl(filename);
  // 1) the SceneGraph should have a single Shape node a child
  Shape* shape = new Shape();
  wrl.addChild(shape);
  // 2) the Shape node should have an Appearance node in its appearance field
  Appearance* appearance = new Appearance();
  shape->setAppearance(appearance);
  shape->setName(shapeName);
  // 3) the Appearance node should have a Material node in its material field
  Material* material = new Material();

  // TODO Fri Mar 10 15:24:04 2023
  // - set as a static class variable or pass it as argument
  Color c(1.0,0.7,0.3);

  material->setDiffuseColor(c);
  appearance->setMaterial(material);
  // 4) the Shape node should have an IndexedFaceSet node in its geometry node
  IndexedFaceSet* ifs = new IndexedFaceSet();
  shape->setGeometry(ifs);
  // return the IndexedFaceSet pointer
  return ifs;
}

bool LoaderOff::load(const char* filename, SceneGraph& wrl) {
  bool success = false;

  FILE* fp = (FILE*)0;
  try {
    // open the file
    if(filename==(char*)0) throw new StrException("filename==null");

    fp = fopen(filename,"r");
    if(fp==(FILE*)0)
      throw new StrException("unable to OFF open file");
        
    // use the io/TokenizerFile class to parse the input ascii file
    TokenizerFile tkn(fp);
    // read the first token
    if(tkn.get()==false)
      throw new StrException("unable to read first token");

    // TODO Fri Mar 10 14:35:24 2023
    // - generalize to [ST][C][N][4][n]OFF
    if(tkn!="OFF")
      throw new StrException("not an OFF file");
    int nVertices=0;
    if(tkn.getInt(nVertices)==false || nVertices<=0)
      throw new StrException("error parsing nVertices");
    int nFaces=0;
    if(tkn.getInt(nFaces)==false || nFaces<=0)
      throw new StrException("error parsing nFaces");
    int nEdges=0;
    if(tkn.getInt(nEdges)==false /* || nEdges<=0 */)
      throw new StrException("error parsing nEdges");
      
    // create the scene graph structure :
    IndexedFaceSet* ifs        = _initializeSceneGraph(wrl,"SURFACE");
    vector<float>&  coord      = ifs->getCoord();
    vector<int>&    coordIndex = ifs->getCoordIndex();

    int iV;
    float x,y,z;
    for(iV=0;iV<nVertices;iV++) {
      if(tkn.getline()==false)
        throw new StrException("vertex getline");
      // parse the line
      TokenizerString tkns(tkn);
      if(!(tkns.getFloat(x) && tkns.getFloat(y) && tkns.getFloat(z)))
        throw new StrException("parsing vertices");
      coord.push_back(x); coord.push_back(y); coord.push_back(z);
    }

    int iF,niF,j;
    for(iF=0;iF<nFaces;iF++) {
      // one line per face
      if(tkn.getline()==false)
        throw new StrException("face getline");
      // parse the line
      TokenizerString tkns(tkn);
      if(tkns.getInt(niF)==false || niF<3)
        throw new StrException("parsing face size");
      for(j=0;j<niF;j++) {
        if(tkns.getInt(iV)==false || niF<3)
          throw new StrException("parsing face vertex index");
        coordIndex.push_back(iV);
      }
      coordIndex.push_back(-1);

      // TODO Fri Mar 10 14:44:45 2023
      // parse colorspec from rest of the line, which may have 0 to 4 tokens
    }

    if(_addFaceNormals) {
      ifs->clearNormal();
      ifs->setNormalPerVertex(false);
      vector<float>& normal = ifs->getNormal();
      Geometry::computeNormalsPerFace(coord,coordIndex,normal);
    }

    success = true;

    // close the file (this statement may not be reached)
    fclose(fp);
 
  } catch(StrException* e) { 

    if(fp!=(FILE*)0) fclose(fp);
    fprintf(stderr,"LoaderOff | ERROR | %s\n",e->what());
    delete e;
    wrl.clear();
    wrl.setUrl("");

  }

  return success;
}

//------------------------------------------------------------------------
//  Copyright (C) Gabriel Taubin
//  Time-stamp: <2025-08-06 09:57:27 taubin>
//------------------------------------------------------------------------
//
// LoaderStl.cpp
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
#include "LoaderStl.hpp"
#include "StrException.hpp"

#include "wrl/Shape.hpp"
#include "wrl/Appearance.hpp"
#include "wrl/Material.hpp"
#include "wrl/IndexedFaceSet.hpp"

// reference
// https://en.wikipedia.org/wiki/STL_(file_format)

const char* LoaderStl::_ext = "stl";

IndexedFaceSet* LoaderStl::_initializeSceneGraph
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
  Color c(1.0,0.0,0.0); // RED
  material->setDiffuseColor(c);
  appearance->setMaterial(material);
  // 4) the Shape node should have an IndexedFaceSet node in its geometry node
  IndexedFaceSet* ifs = new IndexedFaceSet();
  shape->setGeometry(ifs);
  // return the IndexedFaceSet pointer
  return ifs;
}

bool LoaderStl::_loadFacetAscii
(TokenizerFile& tkn, Vec3f& n, Vec3f& v1, Vec3f& v2, Vec3f& v3) {

  // - parse one facet :
  //
  // facet normal ni nj nk
  //   outer loop
  //     vertex v1x v1y v1z
  //     vertex v2x v2y v2z
  //     vertex v3x v3y v3z
  //   endloop
  // endfacet
  //
  // - return true if successful, false if not

  // facet normal ni nj nk
  if(tkn.expecting("facet")==false) return false;
  if(tkn.expecting("normal")==false) return false;
  if(tkn.getVec3f(n)==false) return false;
  //   outer loop
  if(tkn.expecting("outer")==false) return false;
  if(tkn.expecting("loop")==false) return false;
  //     vertex v1x v1y v1z
  if(tkn.expecting("vertex")==false) return false;
  if(tkn.getVec3f(v1)==false) return false;
  //     vertex v2x v2y v2z
  if(tkn.expecting("vertex")==false) return false;
  if(tkn.getVec3f(v2)==false) return false;
  //     vertex v3x v3y v3z
  if(tkn.expecting("vertex")==false) return false;
  if(tkn.getVec3f(v3)==false) return false;
  //   endloop
  if(tkn.expecting("endloop")==false) return false;
  // endfacet
  if(tkn.expecting("endfacet")==false) return false;
  return true;
}

bool LoaderStl::_loadFacetBinary
(FILE* fp, Vec3f& n, Vec3f& v1, Vec3f& v2, Vec3f& v3, uint16_t* abc) {

  if(fread(&(n[0]),1,12,fp)<12)
    throw new StrException("unable to read normal vector");
  if(fread(&(v1[0]),1,12,fp)<12)
    throw new StrException("unable to read vertex 0");
  if(fread(&(v2[0]),1,12,fp)<12)
    throw new StrException("unable to read vertex 1");
  if(fread(&(v3[0]),1,12,fp)<12)
    throw new StrException("unable to read vertex 2");
  if(fread(abc,1,2,fp)<2)
    throw new StrException("unable to read atribute byte count");

  return true;
}

bool LoaderStl::load(const char* filename, SceneGraph& wrl) {
  bool success = false;

  FILE* fp = (FILE*)0;
  try {
    // open the file
    if(filename==(char*)0) throw new StrException("filename==null");

    // allocate binary header and initialize to zero
    char header[80];
    memset(header,0x00,80);
    // determine if file is ascii or binary
    fp = fopen(filename,"rb");
    if(fp==(FILE*)0)
      throw new StrException("unable to open file for binary read");
    if(fread(header,1,5,fp)<5)
      throw new StrException("unable to read first characters of file");
    bool binary = (strncmp(header,"solid",5)!=0);
    if(binary) {
      // read the rest of the header
      if(fread(header+5,1,75,fp)<75)
        throw new StrException("unable to read 75 next characters of file");
      // read number of triangles
      uint32_t nTriangles = 0;
      if(fread(&nTriangles,1,4,fp)<4)
        throw new StrException("unable to read number of triangles");

      IndexedFaceSet* ifs = _initializeSceneGraph(wrl,"SURFACE");
      // get references to the coordIndex, coord, and normal arrays
      vector<int>& coordIndex = ifs->getCoordIndex();
      vector<float>& coord    = ifs->getCoord();
      vector<float>& normal   = ifs->getNormal();
      // 6) set the normalPerVertex variable to false (i.e., normals per face)  
      ifs->setNormalPerVertex(false);

      int   iV0,iV1,iV2,iT;
      Vec3f n,v1,v2,v3;
      uint16_t abc; // attribute byte count
        for(iT=0;iT<static_cast<int>(nTriangles);iT++) {
        _loadFacetBinary(fp,n,v1,v2,v3,&abc);
        normal.push_back(n[0]);
        normal.push_back(n[1]);
        normal.push_back(n[2]);
        iV0 = static_cast<int>(coord.size()/3);
        coord.push_back(v1[0]);
        coord.push_back(v1[1]);
        coord.push_back(v1[2]);
        iV1 = static_cast<int>(coord.size()/3);
        coord.push_back(v2[0]);
        coord.push_back(v2[1]);
        coord.push_back(v2[2]);
        iV2 = static_cast<int>(coord.size()/3);
        coord.push_back(v3[0]);
        coord.push_back(v3[1]);
        coord.push_back(v3[2]);
        coordIndex.push_back(iV0);
        coordIndex.push_back(iV1);
        coordIndex.push_back(iV2);
        coordIndex.push_back(-1);
      }
      
      success = true;

      fclose(fp);
    } else /* if(ascii) */ {
      // close the binary file and reopen it
      fclose(fp);
      fp = fopen(filename,"r");
      if(fp==(FILE*)0)
        throw new StrException("unable to open ASCII STL file");
        
      // use the io/TokenizerFile class to parse the input ascii file
      TokenizerFile tkn(fp);
      // first token should be "solid"
      if(tkn.expecting("solid")==false)
        throw new StrException("not an ASCII STL file");
      // second token should be the solid name
      if(tkn.get()==false)
        throw new StrException("unable to get solid name");
      string stlName = tkn; // second token should be the solid name

      // create the scene graph structure :
      IndexedFaceSet* ifs = _initializeSceneGraph(wrl,filename);
      // get references to the coordIndex, coord, and normal arrays
      vector<int>& coordIndex = ifs->getCoordIndex();
      vector<float>& coord    = ifs->getCoord();
      vector<float>& normal   = ifs->getNormal();
      // set the normalPerVertex variable to false (i.e., normals per face)  
      ifs->setNormalPerVertex(false);

      int   iV0,iV1,iV2;
      Vec3f n,v1,v2,v3;
      while(_loadFacetAscii(tkn,n,v1,v2,v3)) {
        normal.push_back(n[0]);
        normal.push_back(n[1]);
        normal.push_back(n[2]);
        iV0 = static_cast<int>(coord.size()/3);
        coord.push_back(v1[0]);
        coord.push_back(v1[1]);
        coord.push_back(v1[2]);
        iV1 = static_cast<int>(coord.size()/3);
        coord.push_back(v2[0]);
        coord.push_back(v2[1]);
        coord.push_back(v2[2]);
        iV2 = static_cast<int>(coord.size()/3);
        coord.push_back(v3[0]);
        coord.push_back(v3[1]);
        coord.push_back(v3[2]);
        coordIndex.push_back(iV0);
        coordIndex.push_back(iV1);
        coordIndex.push_back(iV2);
        coordIndex.push_back(-1);
      }

      success = true;

      // close the file (this statement may not be reached)
      fclose(fp);
    }
 
  } catch(StrException* e) { 

    if(fp!=(FILE*)0) fclose(fp);
    fprintf(stderr,"LoaderStl | ERROR | %s\n",e->what());
    delete e;
    wrl.clear();
    wrl.setUrl("");

  }

  return success;
}

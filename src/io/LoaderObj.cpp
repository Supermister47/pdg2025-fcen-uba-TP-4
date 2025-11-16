//------------------------------------------------------------------------
//  Copyright (C) Gabriel Taubin
//  Time-stamp: <2025-08-06 09:57:22 taubin>
//------------------------------------------------------------------------
//
// LoaderObj.cpp
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
#include "LoaderObj.hpp"
#include "StrException.hpp"

#include "core/Geometry.hpp"
#include "wrl/Shape.hpp"
#include "wrl/Appearance.hpp"
#include "wrl/Material.hpp"
#include "wrl/IndexedFaceSet.hpp"

// reference
// https://en.wikipedia.org/wiki/Wavefront_.obj_file

// TODO Sat Mar 11 14:35:35 2023
// - support full spec

const char* LoaderObj::_ext = "obj";

// static
bool LoaderObj::_addFaceNormals = true;

// static
void LoaderObj::setAddFaceNormals(const bool value) {
  _addFaceNormals = value;
}

// static
IndexedFaceSet* LoaderObj::_initializeSceneGraph
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

bool LoaderObj::load(const char* filename, SceneGraph& wrl) {
  bool success = false;

  FILE* fp = (FILE*)0;
  try {
    // open the file
    if(filename==(char*)0) throw new StrException("filename==null");

    fp = fopen(filename,"r");
    if(fp==(FILE*)0)
      throw new StrException("unable to OBJ open file");
        
    // use the io/TokenizerFile class to parse the input ascii file
    TokenizerFile tknf(fp);

    // TODO Sat Mar 11 14:32:35 2023
    // - verify that file starts with a comment block which states
    //   that this is an OBJ file
    // - the following commented code is not general enough
    //
    // bool isOBJ = false;
    // while(isOBJ==false) {
    //   if(tknf.getline()==false)
    //     throw new StrException("unable to read line of file");
    //   if(tkn[0]!='#')
    //     throw new StrException("first line of file does not start with #");
    //   TokenizerString tkns(tknf.substr(1));
    //   int nTokensLine = 0;
    //   while(tkns.get()) {
    //     nTokensLine++;
    //     if(tkns=="OBJ") {
    //       isOBJ = true;
    //       break;
    //     }
    //   } // while(tkns.get())
    //   if(nTokensLine>0) break;
    // }
    // if(isOBJ==false)    
    //   throw new StrException("doesn not look like an OBJ file");

    tknf.setSkipComments(true);
      
    // create the scene graph structure :
    IndexedFaceSet* ifs        = _initializeSceneGraph(wrl,"SURFACE");
    vector<float>&  coord      = ifs->getCoord();
    vector<int>&    coordIndex = ifs->getCoordIndex();

    vector<int> face;
    while(tknf.get()) {
      if(tknf=="v") {
        float x,y,z;
        if(!(tknf.getFloat(x) && tknf.getFloat(y) && tknf.getFloat(z)))
          throw new StrException("unable to parse x,y,z");
        // if a w coordinate is present, it is skipped
        coord.push_back(x); coord.push_back(y); coord.push_back(z);
      } else if(tknf=="f") {
        // since the number of face corners could be of any length >=3
        // get the rest of the line
        if(tknf.getline()==false)
          throw new StrException("unable to get rest of face line");
        TokenizerString tkns(tknf);
        int k;
        face.clear();
        while(tkns.getInt(k))
          face.push_back(k-1); // vertex indices start at 1 in OBJ files
        if(face.size()<3)
          throw new StrException("found face with less than 3 vertices");
        for(k=0;k<static_cast<int>(face.size());k++)
          coordIndex.push_back(face[k]);
        coordIndex.push_back(-1);
      } else {
        // TODO Sat Mar 11 14:36:15 2023
        // - parse other elements, such as vn and vt
        throw new StrException("found unsupported token");
      }
    }

    // while(tknf.getline()) {
    //   if(tknf[0]=='#') continue;  // Tokenizer::getline() does not skip comments
    //   TokenizerString tkns(tknf);
    //   if(tkns.get()==false)
    //     throw new StrException("unable to get first token from line");
    //   if(tkns=="v") {
    //     float x,y,z;
    //     if(!(tkns.getFloat(x) && tkns.getFloat(y) && tkns.getFloat(z)))
    //       throw new StrException("unable to parse x,y,z");
    //     // if a w coordinate is present, it is skipped
    //     coord.push_back(x); coord.push_back(y); coord.push_back(z);
    //   } else if(tkns=="f") {
    //     int k;
    //     face.clear();
    //     while(tkns.getInt(k))
    //       face.push_back(k-1); // vertex indices start at 1 in OBJ files
    //     if(face.size()<3)
    //       throw new StrException("found face with less than 3 vertices");
    //     for(k=0;k<static_cast<int>(face.size());k++)
    //       coordIndex.push_back(face[k]);
    //     coordIndex.push_back(-1);
    //   } else {
    //
    //     // TODO Sat Mar 11 14:36:15 2023
    //     // - parse other element types
    //
    //     throw new StrException("found unsupported token");
    //   }
    // }

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
    fprintf(stderr,"LoaderObj | ERROR | %s\n",e->what());
    delete e;
    wrl.clear();
    wrl.setUrl("");

  }

  return success;
}

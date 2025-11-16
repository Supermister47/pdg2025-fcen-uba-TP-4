//------------------------------------------------------------------------
//  Copyright (C) Gabriel Taubin
//  Time-stamp: <2025-10-14 17:23:34 taubin>
//------------------------------------------------------------------------
//
// dgpTest4.cpp
//
// Software developed for the course
// Digital Geometry Processing
// Copyright (c) 2025, Gabriel Taubin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Brown University nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL GABRIEL TAUBIN BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <string>
#include <iostream>

using namespace std;

#include <wrl/SceneGraphTraversal.hpp>

#include <io/AppLoader.hpp>
#include <io/AppSaver.hpp>
#include <io/LoaderObj.hpp>
#include <io/LoaderOff.hpp>
#include <io/LoaderPly.hpp>
#include <io/LoaderStl.hpp>
#include <io/LoaderWrl.hpp>
#include <io/SaverOff.hpp>
#include <io/SaverPly.hpp>
#include <io/SaverStl.hpp>
#include <io/SaverWrl.hpp>

// #include <core/PolygonMesh.hpp>
#include <core/PolygonMeshTest.hpp>
#include <core/Optimization.hpp>
#include "dgpPrt.hpp"

// Build your own command line program to test your code

int v;
string method;
float mu = 0.5f;
float lambda = 0.5f;
int steps = 2;
float jacobiWeightSmoothing = 0.5f;
float integrateNormalsWeightData = 1.0f;
float integrateNormalsWeightLink = 1.0f;
float integrateNormalsWeightSmoothing = 1.0f;

enum Operation {
  NONE,
  LAPLACIAN_SMOOTHING_VERTEX_COORDINATES,
  LAPLACIAN_SMOOTHING_FACE_NORMALS,
  JACOBI,
  INTEGRATE_NORMALS,
  COLLAPSE_EDGES,
  SPLIT_EDGES,
};

class Data {
public:
  bool      _debug;
  bool      _binaryOutput;
  bool      _removeProperties;
  Operation _operation;
  string    _inFile;
  string    _outFile;
public:
  Data():
    _debug(false),
    _binaryOutput(false),
    _removeProperties(false),
    _operation(NONE),
    _inFile(""),
    _outFile("")
  { }
};

void options(Data& D) {
  cout << "   -d|-debug               [" << tv(D._debug)            << "]" << endl;
  cout << "   -b|-binaryOutput        [" << tv(D._binaryOutput)     << "]" << endl;
  cout << "   -r|-removeProperties    [" << tv(D._removeProperties) << "]" << endl;

  cout << "   -op| POSSIBLE OPERATIONS:" << endl;

    cout << "   laplacianSmoothingVertexCoordinates" << endl;
    cout << "        -mu <float> : smoothing factor (default 0.5)" << endl;
    cout << "        -lambda <float> : smoothing factor (default 0.5)" << endl;
    cout << "        -steps <int> : number of smoothing iterations (default 2)" << endl;

    cout << "   laplacianSmoothingFaceNormals" << endl;
    cout << "        -mu <float> : smoothing factor (default 0.5)" << endl;
    cout << "        -lambda <float> : smoothing factor (default 0.5)" << endl;
    cout << "        -steps <int> : number of smoothing iterations (default 2)" << endl;

    cout << "   jacobi" << endl;
    cout << "       requires additional parameters:" << endl;
    cout << "              -mu <float> : smoothing factor (default 0.5)" << endl;
    cout << "              -lambda <float> : smoothing factor (default 0.5)" << endl;
    cout << "              -steps <int> : number of smoothing iterations (default 2)" << endl;
    cout << "              -jacobiWeightSmoothing <float> : weight for the smoothing term (default 0.5)" << endl;
    cout << "                   *jacobiWeightRegularization will be 1 - jacobiWeightSmoothing" << endl;
    cout << "   integrateNormals" << endl;
    cout << "       requires additional parameters:" << endl;
    cout << "              -mu <float> : smoothing factor (default 0.5)" << endl;
    cout << "              -lambda <float> : smoothing factor (default 0.5)" << endl;
    cout << "              -steps <int> : number of smoothing iterations (default 2)" << endl;
    cout << "              -integrateNormalsWeightData <float> : weight for the data term (default 1.0)" << endl;
    cout << "              -integrateNormalsWeightLink <float> : weight for the link term (default 1.0)" << endl;
    cout << "              -integrateNormalsWeightSmoothing <float> : weight for the smoothing term (default 1.0)" << endl;
    cout << "   collapseEdges" << endl;
    cout << "       requires additional parameters:" << endl;
    cout << "              -m <edge-lengths/garland-heckbert>: whether it uses edge lengths or Garland-Heckbert error metric." << endl;
    cout << "              -v <2/4/8>: to specify the size of the independent vertices that have to be collapsed." << endl;         
    cout << "   splitEdges" << endl;
    cout << "       requires additional parameters:" << endl;
    cout << "              -t <long-threshold/all-edges>: whether it should subdivided long edges or all edges." << endl;
  

}

void usage(Data& D) {
  cout << "USAGE: dgpTest4 [options] inFile outFile" << endl;
  cout << "   -h|-help" << endl;
  options(D);
  cout << endl;
  exit(0);
}

void error(const char *msg) {
  cout << "ERROR: dgpTest4 | " << ((msg)?msg:"") << endl;
  exit(0);
}

//////////////////////////////////////////////////////////////////////
int main(int argc, char **argv) {

  Data D;

  if(argc==1) usage(D);

  for(int i=1;i<argc;i++) {
    if(string(argv[i])=="-h" || string(argv[i])=="-help") {
      usage(D);
    } else if(string(argv[i])=="-d" || string(argv[i])=="-debug") {
      D._debug = !D._debug;
    } else if(string(argv[i])=="-b" || string(argv[i])=="-binaryOutput") {
      D._binaryOutput = !D._binaryOutput;
    } else if(string(argv[i])=="-r" || string(argv[i])=="-removeProperties") {
      D._removeProperties = !D._removeProperties;

    } else if(string(argv[i])=="-op") {
      i++;
      if(i>=argc) error("missing operation after -op");
      string operation = string(argv[i]);
      i++; // Move to the first parameter

      if(operation=="laplacianSmoothingVertexCoordinates" || operation=="laplacianSmoothingFaceNormals") {
        if (operation=="laplacianSmoothingVertexCoordinates") {
          D._operation = LAPLACIAN_SMOOTHING_VERTEX_COORDINATES;
        } else {
          D._operation = LAPLACIAN_SMOOTHING_FACE_NORMALS;
        }

        while(i<argc && string(argv[i])[0]=='-') {
          if(string(argv[i])=="-mu") {
            i++; if(i>=argc) error("missing value after -mu");
            mu = stof(string(argv[i]));
          } else if(string(argv[i])=="-lambda") {
            i++; if(i>=argc) error("missing value after -lambda");
            lambda = stof(string(argv[i]));
          } else if(string(argv[i])=="-steps") {
            i++; if(i>=argc) error("missing value after -steps");
            steps = stoi(string(argv[i]));
          } else {
            error("unknown option for this operation");
          }
          i++;
        }
        i--; // Adjust because main loop will increment

      } else if (operation=="jacobi") {
        D._operation = JACOBI;
        while(i<argc && string(argv[i])[0]=='-') {
          if(string(argv[i])=="-mu") {
            i++; if(i>=argc) error("missing value after -mu");
            mu = stof(string(argv[i]));
          } else if(string(argv[i])=="-lambda") {
            i++; if(i>=argc) error("missing value after -lambda");
            lambda = stof(string(argv[i]));
          } else if(string(argv[i])=="-steps") {
            i++; if(i>=argc) error("missing value after -steps");
            steps = stoi(string(argv[i]));
          } else if (string(argv[i])=="-jacobiWeightSmoothing") {
            i++; if(i>=argc) error("missing value after -jacobiWeightSmoothing");
            jacobiWeightSmoothing = stof(string(argv[i]));
          } else {
            error("unknown option for jacobi");
          }
          i++;
        }
        i--; // Adjust because main loop will increment

      } else if (operation=="integrateNormals") {
        D._operation = INTEGRATE_NORMALS;
        while(i<argc && string(argv[i])[0]=='-') {
          if(string(argv[i])=="-mu") {
            i++; if(i>=argc) error("missing value after -mu");
            mu = stof(string(argv[i]));
          } else if(string(argv[i])=="-lambda") {
            i++; if(i>=argc) error("missing value after -lambda");
            lambda = stof(string(argv[i]));
          } else if(string(argv[i])=="-steps") {
            i++; if(i>=argc) error("missing value after -steps");
            steps = stoi(string(argv[i]));
          } else if (string(argv[i])=="-integrateNormalsWeightData") {
            i++; if(i>=argc) error("missing value after -integrateNormalsWeightData");
            integrateNormalsWeightData = stof(string(argv[i]));
          } else if (string(argv[i])=="-integrateNormalsWeightLink") {
            i++; if(i>=argc) error("missing value after -integrateNormalsWeightLink");
            integrateNormalsWeightLink = stof(string(argv[i]));
          } else if (string(argv[i])=="-integrateNormalsWeightSmoothing") {
            i++; if(i>=argc) error("missing value after -integrateNormalsWeightSmoothing");
            integrateNormalsWeightSmoothing = stof(string(argv[i]));
          } else {
            error("unknown option for integrateNormals");
          }
          i++;
        }
        i--; // Adjust because main loop will increment

      } else if(operation=="collapseEdges") {
        D._operation = COLLAPSE_EDGES;
        bool m_set = false, v_set = false;
        while(i<argc && string(argv[i])[0]=='-') {
          if (string(argv[i])=="-m") {
            i++; if(i>=argc) error("missing value after -m");
            if (string(argv[i])=="edge-lengths") {
              method = "edge-lengths";
            } else if (string(argv[i])=="garland-heckbert") {
              method = "garland-heckbert";
            } else {
              error("unknown method after -m");
            }
            m_set = true;
          } else if (string(argv[i])=="-v") {
            i++; if(i>=argc) error("missing value after -v");
            v = stoi(string(argv[i]));
            if (v!=2 && v!=4 && v!=8) {
              error("unknown value after -v");
            }
            v_set = true;
          } else {
            error("unknown option for collapseEdges");
          }
          i++;
        }
        if (!m_set || !v_set) {
          error("collapseEdges requires both -m and -v parameters");
        }
        i--; // Adjust because main loop will increment

      } else if(operation=="splitEdges") {
        D._operation = SPLIT_EDGES;
        bool t_set = false;
        while(i<argc && string(argv[i])[0]=='-') {
          if (string(argv[i])=="-t") {
            i++; if(i>=argc) error("missing value after -t");
            if (string(argv[i])=="long-threshold") {
              method = "long-threshold";
            } else if (string(argv[i])=="all-edges") {
              method = "all-edges";
            } else {
              error("unknown method after -t");
            }
            t_set = true;
          } else {
            error("unknown option for splitEdges");
          }
          i++;
        }
        if (!t_set) {
          error("splitEdges operation requires the -t parameter");
        }
        i--; // Adjust because main loop will increment

      } else {
        error("unknown operation after -op");
      }

      // - add code to parse the desired operation to be performed
      // - from the command line

    } else if(string(argv[i])[0]=='-') {
      error("unknown option");
    } else if(D._inFile=="") {
      D._inFile = string(argv[i]);
    } else if(D._outFile=="") {
      D._outFile = string(argv[i]);
    }
  }

  if(D._inFile =="") error("no inFile");

  // if D._outFile is not specified then no output file will be written
  // if(D._outFile=="") error("no outFile");

  if(D._debug) {
    cout << "dgpTest4 {" << endl;
    cout << endl;
    options(D);
    cout << endl;
    cout << "  inFile  = " << D._inFile << endl;
    cout << "  outFile = " << D._outFile << endl;
    cout << endl;
    fflush(stderr);
  }

  bool success;

  //////////////////////////////////////////////////////////////////////
  // create loader and saver factories
  AppLoader loaderFactory;
  AppSaver  saverFactory;

  // register input file loaders
  LoaderObj* objLoader = new LoaderObj();
  loaderFactory.registerLoader(objLoader);
  LoaderOff* offLoader = new LoaderOff();
  loaderFactory.registerLoader(offLoader);
  LoaderPly* plyLoader = new LoaderPly();
  loaderFactory.registerLoader(plyLoader);
  LoaderStl* stlLoader = new LoaderStl();
  loaderFactory.registerLoader(stlLoader);
  LoaderWrl* wrlLoader = new LoaderWrl();
  loaderFactory.registerLoader(wrlLoader);

  //  If SaverPly::setDefaultDataType is used, it must be called
  //  before the Saver constructor; otherwise SaverPly::setDataType
  //  should be called after to set the proper value for the private
  //  variable SaverPly::_dataType before this instance of SaverPly is
  //  used

  // register output file savers  
  SaverOff* offSaver = new SaverOff();
  saverFactory.registerSaver(offSaver);
  SaverPly* plySaver = new SaverPly();
  saverFactory.registerSaver(plySaver);
  SaverStl* stlSaver = new SaverStl();
  saverFactory.registerSaver(stlSaver);
  SaverWrl* wrlSaver = new SaverWrl();
  saverFactory.registerSaver(wrlSaver);

  SaverStl::FileType stlFt =
    (D._binaryOutput)?SaverStl::FileType::BINARY:SaverStl::FileType::ASCII;
  stlSaver->setFileType(stlFt);

  Ply::DataType plyDt =
    (D._binaryOutput)?Ply::DataType::BINARY_LITTLE_ENDIAN:Ply::DataType::ASCII;
  plySaver->setDataType(plyDt);

  if(D._debug) {
    SaverPly::setOstream(&cout);
    SaverPly::setIndent("    ");
  }

  //////////////////////////////////////////////////////////////////////
  // read SceneGraph

  SceneGraph wrl; // create empty scene graph

  if(D._debug) {
    cout << "  loading inFile {" << endl;
  }

  success = loaderFactory.load(D._inFile.c_str(),wrl);

  if(D._debug) {
    cout << "    success        = " << tv(success)          << endl;
    cout << "  } loading inFile" << endl;
    cout << endl;
  }

  if(success==false) return -1;

  // I will get the first IndexedFaceSet found in the scene graph
  SceneGraphTraversal sgt(wrl);
  IndexedFaceSet* ifsInput = (IndexedFaceSet*)0;
  Node* node;
  for (int iIfs=0; (node=sgt.next())!=(Node*)0; iIfs++) {
    Shape* shape = dynamic_cast<Shape*>(node);
    if (shape==(Shape*)0) continue;
    IndexedFaceSet* ifs = dynamic_cast<IndexedFaceSet*>(shape->getGeometry());
    if (ifs==(IndexedFaceSet*)0) continue;
    // do something with ifs

    ifsInput = ifs;
    break;
  }

  if(D._removeProperties) {
    if(D._debug) cout << "  removing properties {" << endl;
    Node* node;
    SceneGraphTraversal sgt(wrl);
    for(int iIfs=0;(node=sgt.next())!=(Node*)0;iIfs++) {
      Shape* shape = dynamic_cast<Shape*>(node);
      if(shape==(Shape*)0) continue;
      IndexedFaceSet* ifs = dynamic_cast<IndexedFaceSet*>(shape->getGeometry());
      if(ifs==(IndexedFaceSet*)0) continue;
      ifs->setNormalPerVertex(true);
      ifs->getNormal().clear();
      ifs->getNormalIndex().clear();
      ifs->setColorPerVertex(true);
      ifs->getColor().clear();
      ifs->getColorIndex().clear();
      ifs->getTexCoord().clear();
      ifs->getTexCoordIndex().clear();
    }
    if(D._debug) cout << "  } removing properties" << endl;
    if(D._debug) cout << endl;
  }

  // print PolygonMesh info before processing
  if(D._debug) {
    cout << "  before processing" << endl;
    PolygonMeshTest(wrl,"  ");
    cout << endl;
  }
  
  // process
  
  if(D._debug) cout << "  processing {" << endl;

  // I will consider only the first IndexedFaceSet found
  Optimization optimization = Optimization();
  optimization.setInput(ifsInput);
  IndexedFaceSet* ifsOptimized = new IndexedFaceSet();
  optimization.setOptimized(ifsOptimized);

  switch (D._operation)
  {
  case LAPLACIAN_SMOOTHING_VERTEX_COORDINATES:
    {
      cout << "    LAPLACIAN_SMOOTHING_VERTEX_COORDINATES" << endl; 
      cout << "      mu      = " << mu      << endl;
      cout << "      lambda  = " << lambda  << endl;
      cout << "      steps   = " << steps   << endl;
      optimization.setMu(mu);
      optimization.setLambda(lambda);
      optimization.setSteps(steps);
      optimization.laplacianSmoothingVertexCoordinatesRun();
    }
    break;
  case LAPLACIAN_SMOOTHING_FACE_NORMALS:
    {
      cout << "    LAPLACIAN_SMOOTHING_FACE_NORMALS" << endl;
      cout << "      mu      = " << mu      << endl;
      cout << "      lambda  = " << lambda  << endl;
      cout << "      steps   = " << steps   << endl;
      optimization.setMu(mu);
      optimization.setLambda(lambda);
      optimization.setSteps(steps);
      optimization.laplacianSmoothingFaceNormalsRun(true);
    }
    break;
  case JACOBI:
    {
      cout << "    JACOBI" << endl;
      cout << "      mu                      = " << mu                      << endl;
      cout << "      lambda                  = " << lambda                  << endl;
      cout << "      steps                   = " << steps                   << endl;
      cout << "      jacobiWeightSmoothing   = " << jacobiWeightSmoothing   << endl;
      optimization.setMu(mu);
      optimization.setLambda(lambda);
      optimization.setSteps(steps);
      optimization.setJacobiWeightSmoothing(jacobiWeightSmoothing);
      optimization.jacobiRun();
    }
    break;
  case INTEGRATE_NORMALS:
    {
      cout << "    INTEGRATE_NORMALS" << endl;
      cout << "      mu                              = " << mu                              << endl;
      cout << "      lambda                          = " << lambda                          << endl;
      cout << "      steps                           = " << steps                           << endl;
      cout << "      integrateNormalsWeightData      = " << integrateNormalsWeightData      << endl;
      cout << "      integrateNormalsWeightLink      = " << integrateNormalsWeightLink      << endl;
      cout << "      integrateNormalsWeightSmoothing = " << integrateNormalsWeightSmoothing << endl;
      optimization.setMu(mu);
      optimization.setLambda(lambda);
      optimization.setSteps(steps);
      optimization.setIntegrateNormalsWeightData(integrateNormalsWeightData);
      optimization.setIntegrateNormalsWeightLink(integrateNormalsWeightLink);
      optimization.setIntegrateNormalsWeightSmoothing(integrateNormalsWeightSmoothing);
      optimization.integrateNormalsRun(true);
    }
    break;
  case COLLAPSE_EDGES:
    {
      cout << "    COLLAPSE_EDGES" << endl;
      cout << "      method = " << method << endl;
      cout << "      v      = " << v      << endl;
      Optimization::EdgeCollapseErrorMetric errorMetric;
      if (method=="edge-lengths") {
        errorMetric = Optimization::EdgeCollapseErrorMetric::EDGE_LENGTH;
      } else if (method=="garland-heckbert") {
        errorMetric = Optimization::EdgeCollapseErrorMetric::GARLAND_HECKBERT;
      }
      Optimization::EdgeCollapseIndependentSet indepSet;
      if (v==2) {
        indepSet = Optimization::EdgeCollapseIndependentSet::VERTICES_2;
      } else if (v==4) {
        indepSet = Optimization::EdgeCollapseIndependentSet::VERTICES_4;
      } else if (v==8) {
        indepSet = Optimization::EdgeCollapseIndependentSet::VERTICES_8;
      }
      optimization.collapseEdgesApply(errorMetric, indepSet);
    }
    break;
  case SPLIT_EDGES:
    {
      cout << "    SPLIT_EDGES" << endl;
      cout << "      method = " << method << endl;
      Optimization::SplitEdgesMode splitMode;
      if (method=="long-threshold") {
        splitMode = Optimization::SplitEdgesMode::LONG;
      } else if (method=="all-edges") {
        splitMode = Optimization::SplitEdgesMode::ALL;
      }
      optimization.adaptiveSubdivisionApply(splitMode, true);
    }
    break;
  default:
    break;
  }

  if(D._debug) cout << "  } processing" << endl;

  // update scene graph with optimized IndexedFaceSet
  SceneGraphTraversal sgt2(wrl);
  for(int iIfs=0;(node=sgt2.next())!=(Node*)0;iIfs++) {
    Shape* shape = dynamic_cast<Shape*>(node);
    if(shape==(Shape*)0) continue;
    IndexedFaceSet* ifs = dynamic_cast<IndexedFaceSet*>(shape->getGeometry());
    if(ifs==(IndexedFaceSet*)0) continue;
    ifs = optimization.getOptimized();
    shape->setGeometry(ifs);
  }
  
  // print PolygonMesh info after processing
  if(D._debug) {
    cout << "  after processing" << endl;
    PolygonMeshTest(wrl,"  ");
    cout << endl;
  }
  
  // write output file
  if(D._outFile!="") {

    if(D._debug) {
      cout << "  saving outFile {" << endl;
    }
      
    success = saverFactory.save(D._outFile.c_str(),wrl);
        
    if(D._debug) {
      cout << "    success        = " << tv(success)          << endl;
      cout << "  }" << endl;
      cout << endl;
    }

  } else {
    if(D._debug) {
      cout << "  no outFile written" << endl;
    }
  }
    
  //////////////////////////////////////////////////////////////////////

  if(D._debug) {
    cout << "} dgpTest4" << endl;
    fflush(stderr);
  }

  return 0;
}

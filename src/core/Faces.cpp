//------------------------------------------------------------------------
//  Copyright (C) Gabriel Taubin
//  Time-stamp: <2025-08-04 22:09:56 gtaubin>
//------------------------------------------------------------------------
//
// Faces.cpp
//
// Written by: <Your Name>
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

#include <math.h>
#include <iostream>

#include "Faces.hpp"
  
Faces::Faces(const int nV, const vector<int>& coordIndex) {
    _coordIndex = vector<int>();
    _faceIndex = vector<pair<int, int>>();
    bool isFaceFirstCorner = true;
    pair<int, int> faceData;
    int nF = 0;
    int faceSize = 0;

    _nV = 0;

    for (int i=0; i < coordIndex.size(); i++) {
        int coord = coordIndex[i];

        if (coord == -1) {
            nF++;
            _coordIndex.push_back(-nF);

            isFaceFirstCorner = true;
            faceData.second = faceSize;
            _faceIndex.push_back(faceData);
        }
        else {
            _nV++;
            if (isFaceFirstCorner) {
                faceSize = 0;
                faceData.first = i;
                isFaceFirstCorner = false;
            }
            _coordIndex.push_back(coord);
            faceSize++;
        }
    }

    _nF = nF;

    if (nV != _nV) {
        printf("Faces constructor: The nV passed as parameter is %u, but differs from the actual amount of vertices, which is %u.\n", nV, _nV);
    }
}

int Faces::getNumberOfVertices() const {
  return _nV;
}

int Faces::getNumberOfFaces() const {
  return _nF;
}

int Faces::getNumberOfCorners() const {
  return _coordIndex.size();
}

int Faces::getFaceSize(const int iF) const {
    if (!isValidFaceIndex(iF))
        return 0;

    return _faceIndex[iF].second;
}

int Faces::getFaceFirstCorner(const int iF) const {
    if (!isValidFaceIndex(iF))
        return -1;

    return _faceIndex[iF].first;
}

int Faces::getFaceVertex(const int iF, const int j) const {
    if (!isValidFaceIndex(iF) || j >= getFaceSize(iF))
        return -1;

    int indexFirstFaceCorner = _faceIndex[iF].first;
    return _coordIndex[indexFirstFaceCorner+j];

}

int Faces::getCornerFace(const int iC) const {
    if (!isValidCorner(iC))
        return -1;

    int i = iC+1;
    while (i < _coordIndex.size() && _coordIndex[i] >= 0) {
        i++;
    }

    int faceIndex = abs(_coordIndex[i]);
    return faceIndex-1;

}

int Faces::getNextCorner(const int iC) const {
    if (!isValidCorner(iC))
        return -1;

    int faceIndex = getCornerFace(iC);
    int nextCorner = _coordIndex[iC+1];

    if (nextCorner >= 0) {
        return nextCorner;
    }
    else {
        return getFaceFirstCorner(faceIndex);
    }
}




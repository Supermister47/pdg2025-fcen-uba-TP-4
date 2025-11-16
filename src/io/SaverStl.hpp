//------------------------------------------------------------------------
//  Copyright (C) Gabriel Taubin
//  Time-stamp: <2025-08-06 09:57:34 taubin>
//------------------------------------------------------------------------
//
// SaverStl.hpp
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

#ifndef _SAVER_STL_HPP_
#define _SAVER_STL_HPP_

#include <cstdio>
#include "Saver.hpp"
#include "wrl/IndexedFaceSet.hpp"
#include <stdint.h>


class SaverStl : public Saver {

private:

  const static char* _ext;

public:

  enum FileType {
    ASCII = 0,
    BINARY
  };

  SaverStl()  {};
  ~SaverStl() {};

  bool  save(const char* filename, SceneGraph& wrl) const;
  const char* ext() const { return _ext; }

  static void setFileType(const FileType ft);
  
private:

  static FileType _fileType; // default : ASCII

private:
  
  bool _saveAscii
  (FILE* fp, const char* solidname, IndexedFaceSet& ifs) const;
  bool _saveBinary
  (FILE* fp, const char* solidname, IndexedFaceSet& ifs) const;

};

#endif /* _SAVER_STL_HPP_ */

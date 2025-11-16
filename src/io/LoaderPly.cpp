//------------------------------------------------------------------------
//  Copyright (C) Gabriel Taubin / 3D Shape Tech LLC
//  Time-stamp: <2025-08-06 09:57:26 taubin>
//------------------------------------------------------------------------
//
// LoaderPly.cpp
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

using namespace std;

#include "LoaderPly.hpp"
#include "TokenizerFile.hpp"
#include "TokenizerString.hpp"
#include "StrException.hpp"
#include <wrl/Shape.hpp>
#include <wrl/Appearance.hpp>
#include <wrl/Material.hpp>
#include <wrl/ImageTexture.hpp>
#include <wrl/IndexedFaceSetPly.hpp>

const char* LoaderPly::_ext = "ply";

//////////////////////////////////////////////////////////////////////
// static
Ply::DataType LoaderPly::systemEndian() {
  return
    (Endian::isLittleEndianSystem())?
    Ply::DataType::BINARY_LITTLE_ENDIAN:
    Ply::DataType::BINARY_BIG_ENDIAN;
}

//////////////////////////////////////////////////////////////////////
// static
bool LoaderPly::sameAsSystemEndian(Ply::DataType fileEndian) {
   return (fileEndian==systemEndian());
}

//////////////////////////////////////////////////////////////////////
// static
void LoaderPly::addBinaryValue
(Endian::SingleValueBuffer& buff,
 const Ply::Element::Property::Type propertyType,
 const bool swapBytes,
 void* value) {

  switch(propertyType) {
  case Ply::Element::Property::CHAR:
  case Ply::Element::Property::INT8:
    {
      vector<char>* valueChar= static_cast<vector<char>*>(value);
      char v = buff.c[0];
      valueChar->push_back(v);
    }
    break;
  case Ply::Element::Property::UCHAR:
  case Ply::Element::Property::UINT8:
    {
      vector<uchar>* valueUChar= static_cast<vector<uchar>*>(value);
      uchar v = buff.uc[0];
      valueUChar->push_back(v);
    }
    break;
  case Ply::Element::Property::SHORT:
  case Ply::Element::Property::INT16:
    {
      vector<short>* valueShort= static_cast<vector<short>*>(value);
      if(swapBytes) Endian::swapShort(buff);
      short v = buff.s[0];
      valueShort->push_back(v);
    }
    break;
  case Ply::Element::Property::USHORT:
  case Ply::Element::Property::UINT16:
    {
      vector<ushort>* valueUShort= static_cast<vector<ushort>*>(value);
      if(swapBytes) Endian::swapUShort(buff);
      ushort v = buff.us[0];
      valueUShort->push_back(v);
    }
    break;
  case Ply::Element::Property::INT:
  case Ply::Element::Property::INT32:
    {
      vector<int>* valueInt= static_cast<vector<int>*>(value);
      if(swapBytes) Endian::swapInt(buff);
      int v = buff.i[0];
      valueInt->push_back(v);
    }
    break;
  case Ply::Element::Property::UINT:
  case Ply::Element::Property::UINT32:
    {
      vector<uint>* valueUInt= static_cast<vector<uint>*>(value);
      if(swapBytes) Endian::swapUInt(buff);
      uint v = buff.ui[0];
      valueUInt->push_back(v);
    }
    break;
  case Ply::Element::Property::FLOAT:
  case Ply::Element::Property::FLOAT32:
  case Ply::Element::Property::FLOAT32_2:
  case Ply::Element::Property::FLOAT32_3:
    {
      vector<float>* valueFloat = static_cast<vector<float>*>(value);
      if(swapBytes) Endian::swapFloat(buff);
      float v = buff.f[0];
      valueFloat->push_back(v);
    }
    break;
  case Ply::Element::Property::DOUBLE:
  case Ply::Element::Property::FLOAT64:
    {
      vector<double>* valueDouble = static_cast<vector<double>*>(value);
      if(swapBytes) Endian::swapDouble(buff);
      double v = buff.d[0];
      valueDouble->push_back(v);
    }
    break; 
  case Ply::Element::Property::NONE:
    {
      throw new StrException("unexpected NONE binary value tyde");
    }
  }
}

//////////////////////////////////////////////////////////////////////
// static
void LoaderPly::addAsciiValue
(const string& token,
 const Ply::Element::Property::Type propertyType,
 void* value) {

  switch(propertyType) {
  case Ply::Element::Property::CHAR:
  case Ply::Element::Property::INT8:
    {
      vector<char>* valueChar= static_cast<vector<char>*>(value);
      char v = static_cast<char>(atoi(token.c_str()));
      valueChar->push_back(v);
    }
    break;
  case Ply::Element::Property::UCHAR:
  case Ply::Element::Property::UINT8:
    {
      vector<uchar>* valueUChar= static_cast<vector<uchar>*>(value);
      uchar v = static_cast<uchar>(atoi(token.c_str()));
      valueUChar->push_back(v);
    }
    break;
  case Ply::Element::Property::SHORT:
  case Ply::Element::Property::INT16:
    {
      vector<short>* valueShort= static_cast<vector<short>*>(value);
      short v = static_cast<short>(atoi(token.c_str()));
      valueShort->push_back(v);
    }
    break;
  case Ply::Element::Property::USHORT:
  case Ply::Element::Property::UINT16:
    {
      vector<ushort>* valueUShort= static_cast<vector<ushort>*>(value);
      ushort v = static_cast<ushort>(atoi(token.c_str()));
      valueUShort->push_back(v);
    }
    break;
  case Ply::Element::Property::INT:
  case Ply::Element::Property::INT32:
    {
      vector<int>* valueInt= static_cast<vector<int>*>(value);
      int v = static_cast<int>(atoi(token.c_str()));
      valueInt->push_back(v);
    }
    break;
  case Ply::Element::Property::UINT:
  case Ply::Element::Property::UINT32:
    {
      vector<uint>* valueUInt= static_cast<vector<uint>*>(value);
      uint v = static_cast<uint>(atol(token.c_str()));
      valueUInt->push_back(v);
    }
    break;
  case Ply::Element::Property::FLOAT:
  case Ply::Element::Property::FLOAT32:
  case Ply::Element::Property::FLOAT32_2:
  case Ply::Element::Property::FLOAT32_3:
    {
      vector<float>* valueFloat = static_cast<vector<float>*>(value);
      float v = static_cast<float>(atof(token.c_str()));
      valueFloat->push_back(v);
    }
    break;
  case Ply::Element::Property::DOUBLE:
  case Ply::Element::Property::FLOAT64:
    {
      vector<double>* valueDouble = static_cast<vector<double>*>(value);
      double v = static_cast<double>(atof(token.c_str()));
      valueDouble->push_back(v);
    }
    break;
  case Ply::Element::Property::NONE:
    {
      throw new StrException("unexpected NONE ascii value type");
    }
  }
}

//////////////////////////////////////////////////////////////////////
// returns number of bytes read
size_t LoaderPly::readHeader(FILE* fp, Ply& ply, const string indent) {

  (void) indent;

  size_t nBytes = 0;
  if(fp) {
    TokenizerFile ftkn(fp);

    // read first line
    if(ftkn.getline()==false)
      throw new StrException("cannot read first line");

    // first line should be "ply"
    if(ftkn.compare(0,3,"ply")!=0)
      throw new StrException("this is not a ply file");

    // second line should be "format"
    if(ftkn.getline()==false)
      throw new StrException("cannot read scond line");

    if(ftkn == "format ascii 1.0")
      ply._dataType = Ply::DataType::ASCII;
    else if(ftkn == "format binary_little_endian 1.0")
      ply._dataType = Ply::DataType::BINARY_LITTLE_ENDIAN;
    else if(ftkn == "format binary_big_endian 1.0")
      ply._dataType = Ply::DataType::BINARY_BIG_ENDIAN;
    else
      throw new StrException("cannot parse second line");

    Ply::Element* element = nullptr;

    while(ftkn.get()) {

      if(ftkn.equals("end_header")) {

        break;

      } else if(ftkn=="obj_info") {

        if(ftkn.getline()==false)
          throw new StrException("cannot read rest of line");

        ply.addObjInfo(ftkn);

      } else if(ftkn=="comment") {

        if(ftkn.getline()==false)
          throw new StrException("cannot read rest of line");

        // length("TextureFile")=11
        if(ftkn.compare(0,11,"TextureFile")==0) {

          // skip white space
          ulong i=12; for(;isspace(ftkn[i]);i++);
          // rest of line is file name
          string textureFile = ftkn.substr(i,ftkn.length()-i);
          ply.setTextureFile(textureFile);

        } else {

          if(ply._skipComments==false) {
            ply.addComment(ftkn);
          }

        }

      } else if(ftkn=="element") {

        if(ftkn.get()==false)
          throw new StrException("expecting element name");
         string elementName = ftkn;

        if(ftkn.get()==false)      
          throw new StrException("expecting element nRecords");
        int nRecords = atoi(ftkn.c_str());
        if(nRecords<0)
          throw new StrException("expecting non-negative element nRecords");
    
        element = ply.addElement(elementName,nRecords);

      } else if(ftkn=="property") {

        bool                         list;
        Ply::Element::Property::Type listType;
        Ply::Element::Property::Type propertyType;
        string                       propertyName;
    
        if(ftkn.get()==false)      
          throw new StrException("early end of property");

        if(ftkn.equals("list")) {

          list = true;

          if(ftkn.get()==false)      
            throw new StrException("expecting listType");

          listType = Ply::Element::Property::parseType(ftkn);

          if(ftkn.get()==false)      
            throw new StrException("expecting propertyType");

          propertyType = Ply::Element::Property::parseType(ftkn);

          if(ftkn.get()==false)      
            throw new StrException("expecting property name");

          propertyName = ftkn;

        } else {

          list     = false;
          listType = Ply::Element::Property::Type::NONE;

          propertyType = Ply::Element::Property::parseType(ftkn);

          if(ftkn.get()==false)      
            throw new StrException("expecting property name");

          propertyName = ftkn;
        }

        element->addProperty(propertyName,list,listType,propertyType);

      } else {
        if(ftkn.getline()==false)
          throw new StrException("cannot read rest of line");
      }
    }
    nBytes = static_cast<size_t>(ftell(fp));
  }

  return nBytes;
}

//////////////////////////////////////////////////////////////////////
// static
size_t LoaderPly::readBinaryData(FILE* fp, Ply& ply, const string indent) {

  (void)indent;

  size_t nBytesData = 0;
  if(fp) {
    long fp0 = ftell(fp);

    int                     nElements,iElement,nProperties,iProperty;
    int                     nRecords,iRecord,k0,k1,i;
    int                     nList,nBytesListCount, nBytesListValue;
    int                     nBytesValue,nBytesRead/*,nBytesRecord*/;
    string                  name;
    string                  propertyName;
    void*                   value        = nullptr;
    Ply::DataType           dataType     = ply.getDataType();
    Ply::Element*           element      = nullptr;
    Ply::Element::Property* property     = nullptr;
    Ply::Element::Property::Type propertyType = Ply::Element::Property::Type::NONE;
    Ply::Element::Property::Type listType     = Ply::Element::Property::Type::NONE;

    Endian::SingleValueBuffer buff;

    bool swapBytes = (sameAsSystemEndian(dataType)==false);
    // bool little  = (dataType==Ply::DataType::BINARY_LITTLE_ENDIAN);

    nElements = ply.getNumberOfElements();

    bool wrlMode = ply.getWrlMode();

    for(iElement=0;iElement<nElements;iElement++) {
      element = ply.getElement(iElement);
      name    = element->getName();

      nProperties = element->getNumberOfProperties();
      for(iProperty=0;iProperty<nProperties;iProperty++) {
        property     = element->getProperty(iProperty);
        propertyName = property->getName();
        propertyType = property->getPropertyType();
        if(property->isList()==true) {
          listType = property->getListType();
        }
      }

      nRecords = element->getNumberOfRecords();

      k0 = 0;
      for(iRecord=0;iRecord<nRecords;iRecord++) {
        // nBytesRecord = 0;
        for(iProperty=0;iProperty<nProperties;iProperty++) {

          property     = element->getProperty(iProperty);
          propertyName = property->getName();
          propertyType = property->getPropertyType();
          value        = property->getValue();

          if(property->isList()) {
            nBytesListCount = property->getListTypeSize();
            nBytesListValue = property->getPropertyTypeSize();

            // number of elements in the list

            nBytesRead =
              static_cast<int>
              (fread(&(buff.c),1,static_cast<size_t>(nBytesListCount),fp));

            if(nBytesRead<nBytesListCount) {
              char s[128]; snprintf(s,128,"end of file in record %d",iRecord);
              throw new StrException(string(s));
            }

            // nBytesRecord += nBytesRead;

            nList = 0;
            listType = property->getListType();
            switch(listType) {
            case Ply::Element::Property::Type::CHAR:
            case Ply::Element::Property::Type::INT8:
              nList = static_cast<int>(buff.c[0]&0xff);
              break;
            case Ply::Element::Property::Type::UCHAR:
            case Ply::Element::Property::Type::UINT8:
              nList = static_cast<int>(buff.uc[0]&0xff);
              break;
            case Ply::Element::Property::Type::SHORT:
            case Ply::Element::Property::Type::INT16:
              nList = static_cast<int>(buff.s[0]&0xff);
              break;
            case Ply::Element::Property::Type::USHORT:
            case Ply::Element::Property::Type::UINT16:
              nList = static_cast<int>(buff.us[0]&0xff);
              break;
            case Ply::Element::Property::Type::INT:
            case Ply::Element::Property::Type::INT32:
              nList = static_cast<int>(buff.i[0]&0xff);
              break;
            case Ply::Element::Property::Type::UINT:
            case Ply::Element::Property::Type::UINT32:
              nList = static_cast<int>(buff.ui[0]&0xff);
              break;
            default:
              throw new StrException("unexpected list type");
            }

            // ???
            if(wrlMode && propertyName=="coordIndex")
              property->pushBackList(nList+1);
            else
              property->pushBackList(nList);

            // read nList values, each of length nBytesListValue

            for(i=0;i<nList;i++) {
              nBytesRead =
                static_cast<int>
                (fread(&(buff.c),1,static_cast<size_t>(nBytesListValue),fp));
              if(nBytesRead<nBytesListValue) {
                char s[128]; snprintf(s,128,"end of file in record %d",iRecord);
                throw new StrException(string(s));
              }

              // nBytesRecord += nBytesRead;

              addBinaryValue(buff,propertyType,swapBytes,value);
            }

            if(wrlMode && propertyName=="coordIndex")
              static_cast<vector<int>*>(value)->push_back(-1);

          } else /* if(property.isList()==false) */ {

            nBytesValue = property->getPropertyTypeSize();
            if(wrlMode) {
              int n =
                (propertyType==Ply::Element::Property::Type::FLOAT32_3)?3:
                (propertyType==Ply::Element::Property::Type::FLOAT32_2)?2:1;
            
              if(propertyName=="color")
                nBytesValue = 1;
              else
                nBytesValue /= n;

              while((--n)>=0) {
                nBytesRead =
                  static_cast<int>
                  (fread(&(buff.c),1,static_cast<size_t>(nBytesValue),fp));
                if(nBytesRead<nBytesValue) {
                  char s[128]; snprintf(s,128,"end of file in record %d",iRecord);
                  throw new StrException(string(s));
                }

                // nBytesRecord += nBytesRead;

                if(wrlMode && propertyName=="color") {
                  vector<float>* colorValue = static_cast<vector<float>*>(value);
                  colorValue->push_back(static_cast<float>(buff.uc[0])/255.0f);
                } else {
                  addBinaryValue(buff,propertyType,swapBytes,value);
                }
              }

            } else {

              nBytesRead =
                static_cast<int>
                (fread(&(buff.c),1,static_cast<size_t>(nBytesValue),fp));
              if(nBytesRead<nBytesValue) {
                char s[128]; snprintf(s,128,"end of file in record %d",iRecord);
                throw new StrException(string(s));
              }

              // nBytesRecord += nBytesRead;

              addBinaryValue(buff,propertyType,swapBytes,value);
            }
          }

        } // for(iProperty=0;iProperty<nProperties;iProperty++)

        // report progress
        k1 = (10*(iRecord+1))/nRecords;
        if(k1>k0) {
          char s[128];
          snprintf(s,128,"        %3d %%",10*k1);
          k0 = k1;
        }

      } // } for(iRecord=0;iRecord<nRecords;iRecord++)
    } // } for(iElement=0;iElement<nElements;iElement++)

    long fp1 = ftell(fp);
    nBytesData = static_cast<size_t>(fp1-fp0);
  }

  return nBytesData;
}

//////////////////////////////////////////////////////////////////////
// static
size_t LoaderPly::readAsciiData(FILE* fp, Ply& ply, const string indent) {

  (void)indent;

  size_t nBytes = 0;
  if(fp) {
    long fp0 = ftell(fp);
    TokenizerFile ftkn(fp);

    int nElements = ply.getNumberOfElements();
    Ply::Element* element;
    Ply::Element::Property* property;
    // Ply::Element::Property::Type listType = Ply::Element::Property::Type::NONE;
    Ply::Element::Property::Type propertyType = Ply::Element::Property::Type::NONE;
    void* value;
    string line,name,propertyName,token;
    int i,iElement,iProperty,iRecord,k0,k1,nList,nProperties,nRecords;

    bool wrlMode = ply.getWrlMode();

    for(iElement=0;iElement<nElements;iElement++) {
       element = ply.getElement(iElement);

       name    = element->getName();
       nProperties = element->getNumberOfProperties();
       nRecords = element->getNumberOfRecords();

       k0 = 0;
       for(iRecord=0;iRecord<nRecords;iRecord++) {

          // one record per line
          if(ftkn.getline()==false) {
            char s[128]; snprintf(s,128,"found empty record %d",iRecord);
            throw new StrException(string(s));
          }

          TokenizerString stkn(ftkn);

          for(iProperty=0;iProperty<nProperties;iProperty++) {

            property     = element->getProperty(iProperty);
            propertyName = property->getName();
            propertyType = property->getPropertyType();
  
            if(property->isList()==true) {
 
              nList = 0;

              if(stkn.get()==false) {
                char s[128];
                snprintf(s,128,"end of line in property record %d",iRecord);
                throw new StrException(string(s));
              }

              nList = atoi(stkn.c_str());

              // Sun Feb 26 17:31:14 2023 ???
              if(wrlMode && propertyName=="coordIndex")
                property->pushBackList(nList+1);
              else
                property->pushBackList(nList);
 
              value = property->getValue();
  
               for(i=0;i<nList;i++) {
                 if(stkn.get()==false) {
                   char s[128];
                   snprintf(s,128,"end of line in property record %d",iRecord);
                   throw new StrException(string(s));
                 }
                 addAsciiValue(stkn,propertyType,value);
               }

               if(wrlMode && propertyName=="coordIndex")
                 static_cast<vector<int>*>(value)->push_back(-1);

            } else /* if(property.isList()==false) */ {

              value = property->getValue();
 
              int n =
                (propertyType==Ply::Element::Property::Type::FLOAT32_3)?3:
                (propertyType==Ply::Element::Property::Type::FLOAT32_2)?2:1;

              while(--n>=0) {
                if(stkn.get()==false) {
                  char s[128];
                  snprintf(s,128,"end of line in property record %d",iRecord);
                  throw new StrException(string(s));
                }
                addAsciiValue(stkn,propertyType,value);
                if(wrlMode && propertyName=="color") {
                    static_cast<vector<float>*>(value)->back() /= 255.0;
                }
              }
            }
          }

          // report progress
          k1 = (10*(iRecord+1))/nRecords;
          if(k1>k0) {
            char s[128];
            snprintf(s,128,"        %3d %%",10*k1);
            k0 = k1;
          }

      } // for(iRecord=0;iRecord<nRecords;iRecord++)
    } // for(iElement=0;iElement<nElements;iElement++)

    long fp1 = ftell(fp);
    nBytes = static_cast<size_t>(fp1-fp0);
  }

  return nBytes;
}

//////////////////////////////////////////////////////////////////////
// static
bool LoaderPly::load(const char* filename, Ply & ply, const string indent) {

  bool success = false;

  FILE* fp  = nullptr;
  ply.clear();
  try {

    // open the file for ascii reading
    if(filename==nullptr)
      throw new StrException("no filename");
    fp = fopen(filename,"r");
    if(fp==nullptr)
      throw new StrException("unable to open file for ascii reading");

    size_t nBytesHeader = readHeader(fp,ply,indent+"  ");

    /* size_t nBytesData   = 0; */
    if(ply.getDataType()==Ply::DataType::ASCII) {
      // continue reading ascii data from the same FileInputStream
      /* nBytesData = */ readAsciiData(fp,ply,indent+"  ");
      fclose(fp);
    } else /* if(ply.getDataType()==Ply::DataType::BINARY_LITTLE_ENDIAN ::
                 ply.getDataType()==Ply::DataType::BINARY_BIG_ENDIAN) */ {

      fclose(fp);
      fp = fopen(filename,"rb");
      if(fp==nullptr)
        throw new StrException("unable to open file to read binary data");

      // skip header
      if(fseek(fp,static_cast<long>(nBytesHeader),SEEK_SET)!=0)
        throw new StrException("failed to skip header to read binary data");

      /* nBytesData = */ readBinaryData(fp,ply,indent+"  ");

      fclose(fp);
    }

    ply.logInfo(std::cout,indent+"  ");

    success = true;

  } catch(StrException* e) { 
    ply.clear();
    if(fp) fclose(fp);
    delete e;
  }

  return success;
}

//////////////////////////////////////////////////////////////////////
bool LoaderPly::load
(const char* filename, SceneGraph& wrl) {
  (void) wrl;

  const string indent = "";

  bool success = false;

  Ply*  ply = nullptr;
  try {

    ply = new Ply();

    if(load(filename,*ply,"  ")==false)
      throw new StrException("load(const char*,Ply&)==false");

    // TODO Fri Mar 10 15:43:43 2023
    // - what if the scene graph already has a POINTS node ?

    // insert into scene graph
    Shape* s = new Shape();
    s->setName("SURFACE");
    Appearance* a = new Appearance();
    if(ply->getTextureFile()!="") {
      ImageTexture* it = new ImageTexture();
      // TODO : set ImageTexture properties from _ply
      a->setTexture(it);
    } else {
      Material* m = new Material();
      // TODO : set material properties from _ply
      a->setMaterial(m);
    }
    s->setAppearance(a);
    IndexedFaceSetPly* ifsPly = new IndexedFaceSetPly(ply,"  ");
    s->setGeometry(ifsPly);

    wrl.addChild(s);

    success = true;

  } catch(StrException* e) { 
    if(ply) delete ply;
    delete e;
  }

  return success;
}

/*
* @file
* This file is part of the developer version of GADGET3 and contains
* the license conditions for its usage.
*
* @author GADGET-development team, led by Volker Springel and Klaus Dolag.
*
* @section LICENSE
* Copyright (c) 2016, Volker Springel, Klaus Dolag, and all contributing authors
* (see change logs). All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*
* 1. Received source code may be modified and used as convenient.
*
* 2. Redistributions of source code or in binary form is only possible with
*    explicit agreement of the copyright holders.
*
* 3. Redistributions of source code must retain the above copyright notice,
*    this list of conditions and the following disclaimer.
*
* 4. Redistributions in binary form must reproduce the above copyright notice,
*    this list of conditions and the following disclaimer in the documentation
*    and/or other materials provided with the distribution.
*
* 5. Neither the name of the copyright holder nor the names of its
*    contributors may be used to endorse or promote products derived from this
*    software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
* IMPLIED WARRANTIES OF  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
* ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
* LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
* CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
* SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
* CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
* ARISING IN ANY WAY OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
* POSSIBILITY OF SUCH DAMAGE.
*
*/
#ifndef __SOA_H__
#define __SOA_H__

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <string>
#include <map>
#include <memory>
#include <vector>
#include <exception>
#include <stddef.h>


using namespace std;

/*
  data of a single array, e.g. T=float, {.name='rho', .is_density_in=true}
*/

typedef MyAtLeastDouble MyAtLeastDouble3[3];
typedef MyLongDouble MyLongDouble3[3];
typedef float float3[3];


// types of assignment from datastructures
enum SoAAssignType{
  ASSIGN_ADD,
  ASSIGN_MAX
  //to be continued if you need more..
}

//SoA descriptor will store data element (type, dimensions, name, pointer) ecc..
struct SoADescriptor{

  string name;


  size_t size=0;

  bool is_density_out = false;
  size_t densdata_out_displacement;
  SoAAssignType densdata_out_assign;


  bool is_hydro_out = false;
  size_t hydrodata_out_displacement;
  SoAAssignType hydrodata_out_assign;

  shared_ptr<void>* data = NULL;
  bool allocated = false;
};

// manage erros when accessing data
struct MyException : public std::exception
{
  std::string s;
 MyException(std::string ss) : s(ss) {}
  ~MyException() throw () {} // Updated
  const char* what() const throw() { return s.c_str(); }
};

/*
  SoA class. Pass an array of SoAData to the constrcutor.
  Allocate with soa.alloc(size); soa.free();
*/
class SoA{
  map<string, SoADescriptor*> m;
  bool allocated = false;
 public:

/*
  we initialise our SoA with an array of SoADescriptor
*/
  SoA(vector <SoADescriptor> & arrays){
    this->allocated = false;

    for(auto const& v: arrays) {
      this->m[v.name] = &v;
    }
  };
  // access array from key
  shared_ptr<void>*  get(string &prop){
    if (m.count(prop))
      return m[prop]->data;
    else
      throw MyException("SoA: key not found: "+prop);
  };

  //assign from DensDataOut our equivalent
  template <typename DataOut, typename DataType>
    void assign(string &prop, T* DataOut, size_t displacement, int  range_min, int range_max, data_index *DataIndexTable, int mode){
    DataType *a = (DataType *)get(prop);
    for(int i = range_min;i<range_max;i++){
      int place = DataIndexTable[i].Index;
      //we extract the value from the DataOut
      DataType v = (DataType)((byte*)(DataOut+i)+displacement);

      if(mode==0){
	a[place] = v;
      }else{
	a[place]+= v;
      }
    }
  };


  //template <typename DataOut>
  //  void out2particle(DataOut *out, int i, int mode){

  //  };
  /*
    we allocate each mapped arrays with n elements * size bytes;
    f decides which field we allocate
  */
  void alloc(int n, function <void(SoADescriptor*)> f){
    if(this->allocated)
      throw "SoA already allocated.";
    
    for(auto const& v: this->m) {
      SoADescriptor *e = v.second;
      if(f(e)){
	if(e->size==0)
	  throw MyException("SoA: zero size datatype "+e->name);
	
	size_t bytes =  e->size  * n;
	
	e->data = mymalloc(e->name, bytes);
	e->allocated = true;
      }
    }
    this->allocated = true;
  };

  /*
    we allocate all mapped arrays
  */
  void alloc(int n){
    alloc(n, [] (SoADescriptor *d) -> void { return true });
  };
  /*
    deallocate all mapped arrays
  */
  void free(){
    if(!this->allocated)
      throw "SoA not allocated.";

    //we free in reverse order
    for (auto iter = m.rbegin(); iter != m.rend(); ++iter) {
      if(iter->allocated){
	::myfree(iter->second->data);
      }
    }

    this->allocated = false;

  };
};


static void soa_init(){
  vector<SoADescriptor> fields = {
#ifdef AR_FIX_TIME_DEP_ART_COND
    {
      .name="RhoOut", 
      .size=sizeof(MyAtLeastDouble), 
      .is_density_out = true, 
      .densdata_out_displacement = offsetof(DensDataOut, Rho),
      .densdata_out_assign = ASSIGN_ADD
    },
    {
      .name="GradAOut", 
      .size=sizeof(MyAtLeastDouble3),
      .is_density_out = true,
      .densdata_out_displacement = offsetof(DensDataOut, GradAOut),
      .densdata_out_assign = ASSIGN_ADD
    },
#endif
#ifdef AR_FIX_WAKEUP
    {
      .name="MaxSignalVelOut",
      .size=sizeof(MyFloat),
      .is_hydro_out = true,
      .densdata_out_displacement = offsetof(HydroDataOut, MaxSignalVelOut),
      .densdata_out_assign = ASSIGN_MAX
    }
#endif
  };

  soa = new SoA(fields);
};

  static void soa_end(){};

#endif

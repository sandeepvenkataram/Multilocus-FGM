/*
 * locus.cpp
 *
 *  Created on: Mar 16, 2013
 *      Author: sandeep
 */


#include "locus.h"
#include <iostream>
#include "modelFunctions.h"
#include "randomv.h"
using namespace std;

size_t hash_locus::GetHashCodeForBytes(const char * bytes, int numBytes) const
{
   size_t h = 0, g;
   for (int i=0; i<numBytes; i++)
   {
      h = ( h << 4 ) + bytes[i];
      if (g = h & 0xF0000000L) {h ^= g >> 24;}
      h &= ~g;
   }
   return h;
}
size_t hash_locus::GetHashForDouble(double v) const{
   return GetHashCodeForBytes((const char *)&v, sizeof(v));
}
size_t hash_locus::operator()(locus &l ) const{
	size_t ret = 0;
	Matrix<double,Dynamic,1> v = l.getPhenotype();
	for (int i=0; i<v.rows(); i++) ret += ((i+1)*(GetHashForDouble(v(i+1,1))));
	return ret;
	
}


locus locus::mutate(){
	//cout<<"Locus Mutate A"<<endl;
	Matrix<double,Dynamic,1> mutationVector = modelFunctions::getMutationVector(numDimensions);
	//cout<<"Locus Mutate C"<<endl;	
	Matrix<double,Dynamic,1> newPosition = modelFunctions::add(phenotype, mutationVector);
	//cout<<"Locus Mutate D"<<endl;
	locus newLocus(numDimensions);
	//cout<<"Locus Mutate E"<<endl;
	newLocus.setLocusOrder(getLocusOrder());
	newLocus.setPhenotype(newPosition);
	newLocus.setParentId(id);
	//cout<<"Locus Mutate F"<<endl;
	return newLocus;
}


void locus::setPhenotype(Matrix<double,Dynamic,1> phenotype) {
	if(phenotype.rows()!=numDimensions){
		cerr<<"Incorrect Number of Dimensions!"<<endl;		
	}else{
		this->phenotype = phenotype; 
	}
	
}

bool locus::operator== (const locus &x) const{
	return id==x.id;
}

bool locus::operator!= (const locus &x) const{
	return id!=x.id;
}
/*
 * locus.h
 *
 *  Created on: Mar 16, 2013
 *      Author: sandeep
 */

#ifndef LOCUS_H_
#define LOCUS_H_

#include <vector>
#include "modelFunctions.h"
#include "randomv.h"
#include <Eigen/Core>
using namespace Eigen;
using namespace std;

class locus{
	public:
		locus(int numDim){id = modelFunctions::getNextLocusId(); numDimensions=numDim;}
		
		locus mutate(void);
		int getLocusOrder(){	return locusOrder;	}
		void setLocusOrder(int locusOrder) { this->locusOrder = locusOrder;	}
		Matrix<double,Dynamic,1> getPhenotype() { return phenotype;	}
		void setPhenotype(Matrix<double,Dynamic,1> phenotype);
		bool operator== (const locus &x) const;
		bool operator!= (const locus &x) const;
		int getId(void) {return id;}
		void setId(int i) {id = i;}
		int getParentId(void) {return parentId;}
		void setParentId(int i){parentId = i;}

	private:
		Matrix<double,Dynamic,1> phenotype;
		int locusOrder;
		int numDimensions;
		int id;
		int parentId;

};

struct hash_locus{
	size_t operator()(locus &l ) const;
	size_t GetHashForDouble(double v) const;
	size_t GetHashCodeForBytes(const char * bytes, int numBytes) const;
};

#endif /* LOCUS_H_ */

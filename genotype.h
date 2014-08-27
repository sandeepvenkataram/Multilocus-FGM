/*
 * genotype.h
 *
 *  Created on: Mar 16, 2013
 *      Author: sandeep
 */

#ifndef GENOTYPE_H_
#define GENOTYPE_H_
#include "locus.h"
#include "allele.h"
#include "modelFunctions.h"
#include "environment.h"
#include <vector>
#include<Eigen/Core>
using namespace Eigen;
using namespace std;


class genotype{
	public:
		genotype(vector<allele *> al, double myFreq);
		vector<allele *> getAlleles (void) const;
		void setAlleles(vector<allele *> alleles);
		double getFreq(void) const;
		void setFreq(double freq);
		double getW(void) const;
		Matrix<double,Dynamic,1> getPhenotype() const;
		void computeFitness();
		allele * getRecombinedAllele();
		std::pair<genotype,allele *> mutate(int currentAlleleCount, int currentGen);
		size_t getHashValue(void) const {return hashValue;}
		bool operator==(const genotype &x) const;
		bool operator!=(const genotype &x) const;
		string toString() const;

	private:
		vector<allele *> alleles;
		Matrix<double,Dynamic,1> phenotype;
		double freq;
		double fW;
		string myString;
		size_t hashValue;
};

struct hash_genotype{
	size_t operator()(const genotype &g ) const;
	size_t operator()(const genotype *g ) const;
};

#endif /* GENOTYPE_H_ */

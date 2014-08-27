/*
 * allele.h
 *
 *  Created on: Mar 16, 2013
 *      Author: sandeep
 */

#ifndef ALLELE_H_
#define ALLELE_H_
#include <vector>
#include "locus.h"
#include "environment.h"
#include "modelFunctions.h"
#include<Eigen/Core>
using namespace Eigen;
using namespace std;


class allele{
	public:
		allele(void);
		allele(vector<locus> myLoci);

		allele* mutate();

		
		int getBirthday();
		void setBirthday(int birthday);
		bool isFixed();
		void setFixed(bool fixed);
		double getFrequency();
		void setFrequency(double frequency);
		double getW();
		void setW(double w);
		int getId() const;
		void setId(int id);
		vector<locus> getLoci();
		Matrix<double,Dynamic,1> getPosition();
		size_t getHashValue(void) const {return hashValue;}
		bool operator==(const allele &x) const;
		bool operator!=(const allele &x) const;
		string getStringId(void) const {return uniqueId;}
		void computeStringId(void);

	private:
		Matrix<double,Dynamic,1> position;
		vector<locus> loci;
		int id;
		double frequency;
		int birthday;
		double fW;
		bool fixed;
		string uniqueId;
		size_t hashValue;
		void setLoci(vector<locus> loci);
		void computePositionFromLoci();
		void setPosition(Matrix<double,Dynamic,1> position);
		
};

struct hash_allele{
	size_t operator()(const allele &a ) const;
	size_t operator()(const allele *a ) const;
};



#endif /* ALLELE_H_ */

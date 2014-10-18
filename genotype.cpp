/*
 * genotype.cpp
 *
 *  Created on: Mar 16, 2013
 *      Author: sandeep
 */

#include "genotype.h"
#include "randomv.h"
#include "allele.h"
#include <unordered_map>
#include <iostream>
#include <Eigen/Core>
using namespace Eigen;
using namespace std;

genotype::genotype(vector<allele *> al, double myFreq, environment &envRef) {
	alleles = al;
	freq = myFreq;
	computeFitness(envRef); 
	
	vector<int> ids;
	for(vector<allele *>::const_iterator iter = alleles.begin(); iter!=alleles.end(); ++iter){
		ids.push_back((*iter)->getId());
	}
	std::sort(ids.begin(),ids.end());
	myString = "";
	for(int i=0; i<ids.size(); i++){
		std::ostringstream oss2;
		oss2<<ids[i];
		myString += "_" + oss2.str();
	}	
	
	hashValue=0;
	for(int i = 0; i < myString.length(); ++i) {
		hashValue = 65599 * hashValue + myString[i];
	}
}
vector<allele *> genotype::getAlleles(void) const {  return alleles; }
void genotype::setAlleles(vector<allele *> alleles)  { this->alleles = alleles;  }
double genotype::getFreq(void) const {  return freq;  }
void genotype::setFreq(double freq) {  this->freq = freq;  }
double genotype::getW(void)  const {  return fW;   }
Matrix<double,Dynamic,1> genotype::getPhenotype() const { return phenotype;}

size_t hash_genotype::operator()(const genotype &g ) const {
	size_t hash = g.getHashValue();
	return hash ^ (hash >> 16);
}

size_t hash_genotype::operator()(const genotype *g ) const {
	size_t hash = g->getHashValue();
	return hash ^ (hash >> 16);
}

void genotype::computeFitness(environment &envRef){
	Matrix<double,Dynamic,1> res=alleles[0]->getPosition();
	for(int i=1; i<alleles.size(); i++){
		res = modelFunctions::add(res,alleles[i]->getPosition());
	}
	res(0,0)=res(0,0)/alleles.size();
	fW=envRef.fW(res);
	phenotype = res;
}

allele * genotype::getRecombinedAllele(){
	if(alleles.size()!=2){
		int allele = randomv::sampleUniformInt(alleles.size());
		return alleles[allele];
	}else{
		vector<locus> alleleOne;
		vector<locus> alleleTwo;
		alleleOne.push_back(alleles[0]->getLoci()[0]);
		alleleTwo.push_back(alleles[1]->getLoci()[0]);
		bool flipped = false;
		for(int i=1; i<modelFunctions::getNumberLoci();i++){
			if(randomv::sampleUniform() < modelFunctions::getRecombinationDistances()[i-1]){
				flipped = !flipped;
			}
			if(flipped){
				alleleOne.push_back(alleles[1]->getLoci()[i]);
				alleleTwo.push_back(alleles[0]->getLoci()[i]);
			}else{
				alleleOne.push_back(alleles[0]->getLoci()[i]);
				alleleTwo.push_back(alleles[1]->getLoci()[i]);
			}
		}
		if(randomv::sampleUniformInt(2)==1){
			allele * ret = new allele(alleleOne);
			/*if(std::isnan(ret->getW())){
				cout<<"allele is nan in getRecombinedAllele if"<<endl;
			}*/
			return ret;
		}else{
			allele * ret = new allele(alleleTwo);
			/*if(std::isnan(ret->getW())){
				cout<<"allele is nan in getRecombinedAllele else"<<endl;
			}*/
			return ret;
		}
	}

}

std::pair<genotype,allele *> genotype::mutate(int currentAlleleCount, int currentGen, environment &envRef){
	//cout<<"Genotype mutate A"<<endl;
	int al = randomv::sampleUniformInt(alleles.size());
	//cout<<"Genotype mutate B"<<endl;
	allele * mutant = (alleles[al])->mutate();
	//cout<<"Genotype mutate C"<<endl;
	mutant->setId(currentAlleleCount);
	mutant->setBirthday(currentGen);
	//cout<<"Genotype mutate D"<<endl;
	int count=0;
	vector<allele *> newAlleles;
	//cout<<"Genotype mutate E"<<endl;
	for(vector<allele *>::iterator iter = alleles.begin(); iter!=alleles.end(); ++iter){
		if(count==al){
			newAlleles.push_back(mutant);
		}else{
			newAlleles.push_back((*iter));
		}
		count++;
	}
	//cout<<"Genotype mutate F"<<endl;
	genotype g(newAlleles, 0, envRef);
	//cout<<"Genotype mutate G"<<endl;
	std::pair<genotype,allele *> pair (g,mutant);
	//cout<<"Genotype mutate H"<<endl;
	return pair;
}

bool genotype::operator==(const genotype &otherGenotype) const{
	return toString() == otherGenotype.toString();
}

bool genotype::operator!=(const genotype &otherGenotype) const{
	return !(toString()==otherGenotype.toString());
}


string genotype::toString() const{
	return myString;
}
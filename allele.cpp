/*
 * allele.cpp
 *
 *  Created on: Mar 16, 2013
 *      Author: sandeep
 */
#include "allele.h"
#include "locus.h"
#include <iostream>
#include <cmath>
#include <sstream>
using namespace std;

allele::allele(void){

}

int allele::getBirthday() {	return birthday;}
void allele::setBirthday(int birthday) {this->birthday = birthday;}
bool allele::isFixed() {	return fixed;	}
void allele::setFixed(bool fixed) {	this->fixed = fixed;	}
double allele::getFrequency() {	return frequency;	}
void allele::setFrequency(double frequency) {	
	this->frequency = frequency;	
	if(frequency==1.0){
		setFixed(1);
	}else{
		setFixed(0);		
	}	
}
double allele::getW() {	return fW;	}
void allele::setW(double w) {	fW = w;	}
int allele::getId() const{	return id;	}
void allele::setId(int id) {	this->id = id;	}
vector<locus> allele::getLoci() {	return loci;	}
Matrix<double,Dynamic,1> allele::getPosition() {	return position;	}
void allele::setLoci(vector<locus> loci) {	this->loci = loci;	computePositionFromLoci(); computeStringId();}
void allele::setPosition(Matrix<double,Dynamic,1> position) {	this->position = position;	}


size_t hash_allele::operator()(const allele &a ) const{
	//cout<<"Hash Allele Reference A"<<endl;
	size_t hash = a.getHashValue();
	return hash ^ (hash >> 16);
}

size_t hash_allele::operator()(const allele *a ) const{
	//cout<<"Hash Allele Pointer A"<<endl;
	size_t hash = a->getHashValue();
	return hash ^ (hash >> 16);
}

allele::allele(vector<locus> myLoci){
	setLoci(myLoci);
	//cout<<getStringId()<<endl;
	//computeStringId();
}

void allele::computeStringId(){
	//cout<<"Compute String Id A"<<endl;
	uniqueId="";
	std::ostringstream oss;
	oss<<loci[0].getId();
	uniqueId += oss.str();
	//cout<<"Compute String Id B"<<endl;
	for(int i=1; i<loci.size(); i++){
		uniqueId+="_";
		std::ostringstream oss2;
		oss2<<loci[i].getId();
		uniqueId+=oss2.str();
	}
	
	hashValue=0;
	for(int i = 0; i < uniqueId.length(); ++i) {
		hashValue = 65599 * hashValue + uniqueId[i];
	}
	
	//cout<<"Compute String Id C "<<getId()<<"\t"<<getStringId()<<"\t"<<hashValue<<endl;
}

allele* allele::mutate(){ //guaranteed mutation, just need to pick a locus to mutate
	//cout<<"Allele mutate A"<<endl;
	int numberLoci =modelFunctions::getNumberLoci();
	int mutatedLocus = randomv::sampleUniformInt(numberLoci);
	//mutatedLocus=0;
	//cout<<"Allele mutate B"<<endl;
	locus mutant = loci[mutatedLocus].mutate();
	vector<locus> newLoci;
	//cout<<"Allele mutate C"<<endl;
	for(int i=0; i<numberLoci; i++){
		if(i==mutatedLocus){
			newLoci.push_back(mutant);
		}else{
			newLoci.push_back(loci[i]);
		}
	}
	//cout<<"Allele mutate D"<<endl;
	allele * temp = new allele (newLoci);
	/*if(std::isnan(temp->getW())){
		cout<<"allele mutate has NAN!"<<endl;
	}*/
	//cout<<"Allele mutate E"<<endl;
	return temp;
}

void allele::computePositionFromLoci(){
	//cout<<"ComputePositionFromLoci A"<<endl;
	
	int numDimensions = modelFunctions::getNumDimensions();
	int numberLoci =modelFunctions::getNumberLoci();
	Matrix<double,Dynamic,1> cartesianPhenotype(numDimensions,1);
	cartesianPhenotype.fill(0);
	//cout<<"ComputePositionFromLoci B"<<endl;
	//cout<<"ComputePositionFromLoci C "<<numberLoci<<endl;
	for(int i=0; i<numberLoci; i++){
		Matrix<double,Dynamic,1> cartPhen = modelFunctions::polarToCartesian(loci[i].getPhenotype());
		//cout<<"ComputePositionFromLoci C-1"<<endl;
		//cout<<"LOCUS: "<<i<<endl;
		//cout<<"PolarLocus:"<<endl;
		//cout<<loci[i].getPhenotype()<<endl;
		//cout<<"Cart Locus:"<<endl;
		//cout<<cartPhen<<endl;
		cartesianPhenotype = cartesianPhenotype + (cartPhen.transpose() * modelFunctions::getCovarianceMatrices()[i]).transpose();
		//cout<<"ComputePositionFromLoci C-2"<<endl;
	}
	//cout<<"ComputePositionFromLoci D"<<endl;
	
	/*if(cartesianPhenotype != modelFunctions::polarToCartesian(modelFunctions::cartesianToPolar(cartesianPhenotype))){
		cout<<"Polar does not equal cartesian!"<<endl;
		cout<<"ORIGINAL:"<<endl;
		cout<<position<<endl;
		cout<<"Cartesian:"<<endl;
		cout<<cartesianPhenotype<<endl;
		cout<<"BACK TO POLAR:"<<endl;
		cout<<modelFunctions::polarToCartesian(modelFunctions::cartesianToPolar(cartesianPhenotype))<<endl;
		exit(1);
	}*/
	
	position = modelFunctions::cartesianToPolar(cartesianPhenotype);
	//cout<<"ComputePositionFromLoci E"<<endl;
	setW(environment::fW(position));
	
	
	
	/*if(std::isnan(getW())){
		cout<<"allele computing pos from loci has NAN!"<<endl;
		cout<<"CARTESIAN:"<<endl;
		cout<<cartesianPhenotype<<endl;
		cout<<"POLAR:"<<endl;
		cout<<getPosition()<<endl;
		int locNum=0;
		for(vector<locus>::iterator locIter = loci.begin(); locIter!=loci.end(); locIter++){
			cout<<"LOCUS: "<<locNum<<endl;
			cout<<(*locIter).getPhenotype()<<endl;
			locNum++;
		}
	}*/
	//cout<<"ComputePositionFromLoci F"<<endl;
}

bool allele::operator==(const allele &x) const{
	return loci == x.loci;
}
bool allele::operator!=(const allele &x) const{
	return loci != x.loci;
}

#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <gsl/gsl_math.h>
#include "randomv.h"
#include <algorithm>
#include <string>
#include <vector>
#include <iterator>
#include "modelFunctions.h"
#include <Eigen/Core>
using namespace Eigen;
using namespace std;

char modelFunctions::fitnessFunction ='g';
double modelFunctions::maxR = 2.0; 
double modelFunctions::b = 1;
double modelFunctions::a = 1;
double modelFunctions::c = 1;
double modelFunctions::d = 2;
char modelFunctions::mutationFunction = 'e';
bool modelFunctions::walk = false;
double modelFunctions::sigma = 1.0;
double modelFunctions::alpha = 1.0;
double modelFunctions::beta = 1.0;
bool modelFunctions::allowRecombination = false;
int modelFunctions::numberLoci = 1;
int modelFunctions::numDimensions = 2;
vector<double> modelFunctions::parameters;
vector<double> modelFunctions::recombinationDistances;
int modelFunctions::locusId=0;
Matrix<double,Dynamic,Dynamic,RowMajor>  modelFunctions::selectionMatrix;
vector<Matrix<double,Dynamic,Dynamic,RowMajor>>  modelFunctions::covarianceMatrices;

char   modelFunctions::getFitnessFunction(void)  {return fitnessFunction;}
void   modelFunctions::setFitnessFunction(char x){fitnessFunction = x;}
double modelFunctions::getMaxR(void)    { return maxR;}
void   modelFunctions::setMaxR(double x){ maxR = x; }
// for Gaussian
double modelFunctions::getA(void)    { return a; }
void   modelFunctions::setA(double x){ a = x; }
double modelFunctions::getB(void)    { return b; }
void   modelFunctions::setB(double x){ b = x; }
double modelFunctions::getC(void)    { return c; }
void   modelFunctions::setC(double x){ c = x; }
double modelFunctions::getD(void)    { return d; }
void   modelFunctions::setD(double x){ d = x; }
//parameters for new mutation functions
vector<double> modelFunctions::getParameters(void)   { return parameters; }
void   modelFunctions::setParameters(vector<double> x) { parameters = x; }
void   modelFunctions::setMutationFunction(char a)     { mutationFunction = a; }
char   modelFunctions::getMutationFunction(void)       { return mutationFunction; }
//environmental
void   modelFunctions::setWalk(  bool x )   { walk = x; }
void   modelFunctions::setSigma( double x ) { sigma = x; }
void   modelFunctions::setAlpha( double x ) { alpha = x; }
void   modelFunctions::setBeta(  double x ) { beta = x; }
bool   modelFunctions::getWalk(  void )     { return walk; }
double modelFunctions::getSigma( void )     { return sigma; }
double modelFunctions::getAlpha( void )     { return alpha; }
double modelFunctions::getBeta(  void )     { return beta; }
int modelFunctions::getNumberLoci() {	return numberLoci;	}
int modelFunctions::getNumDimensions() {	return numDimensions;	}
void modelFunctions::setNumberLoci(int i){ numberLoci=i;}
void modelFunctions::setNumDimensions(int i) {numDimensions = i;}
vector<double> modelFunctions::getRecombinationDistances() {	return recombinationDistances;	}
void modelFunctions::setRecombinationDistances(vector<double> recombDists){ recombinationDistances = recombDists;}
bool modelFunctions::isAllowRecombination() {	return allowRecombination;	}
void modelFunctions::setAllowRecombination(bool allowR) {allowRecombination = allowR;	}

vector<Matrix<double,Dynamic,Dynamic,RowMajor>> modelFunctions::getCovarianceMatrices(void){	return covarianceMatrices;	}
void modelFunctions::setCovarianceMatrices(vector<Matrix<double,Dynamic,Dynamic,RowMajor>> m){	covarianceMatrices = m;		}
Matrix<double,Dynamic,Dynamic,RowMajor> modelFunctions::getSelectionMatrix(void){	return selectionMatrix;		}
void modelFunctions::setSelectionMatrix(Matrix<double,Dynamic,Dynamic,RowMajor> m){	selectionMatrix=m;	}


void modelFunctions::initialize(){
  //for Gaussian
  b = 1; //mean
  c = 1; //shape
  a = 1; //maximum for s
  d = 2;
}

vector<string> modelFunctions::split(const string& s, const string& delim, const bool keep_empty) {
    vector<string> result;
    if (delim.empty()) {
        result.push_back(s);
        return result;
    }
    string::const_iterator substart = s.begin(), subend;
    while (true) {
        subend = search(substart, s.end(), delim.begin(), delim.end());
        string temp(substart, subend);
        if (keep_empty || !temp.empty()) {
            result.push_back(temp);
        }
        if (subend == s.end()) {
            break;
        }
        substart = subend + delim.size();
    }
    return result;
}


double modelFunctions::fR(){
  double newR;
  switch (mutationFunction){
  case 'g': //GAMMA
    {
      double shape = parameters.front();
      newR = randomv::sampleGamma(shape,1);
    }
    break;
  case 'u': //UNIFORM
    newR = randomv::sampleUniform()*parameters.front();
    break;
  case 'e': //EXPONENTIAL
    newR = randomv::sampleExponential(parameters.front());
    break;
  case 'n': //NORMAL
    {
      newR = randomv::sampleNormal(parameters.at(0),parameters.at(1));
      if (newR<0){
	newR = -newR;
      }
    }
    break;
  default:
    cout<<"unknown function: "<<mutationFunction<<endl;
    exit(1);
  }
  return newR;
}

double modelFunctions::fTheta(){
  double newTheta = randomv::sampleUniform();
  return newTheta*2*M_PI - M_PI; // (-PI,PI]
}

Matrix<double,Dynamic,1> modelFunctions::getMutationVector(int numDimensions){
	double fr = fR();
	
	if(numDimensions == 1){
		Matrix<double,Dynamic,1> mutationVector(numDimensions,1);
		mutationVector(0,0)=fr;
		if(randomv::sampleUniform()<0.5){
				mutationVector(0,0) = -1*mutationVector(0,0);
		}
		return mutationVector;
	}
	
	double direction[numDimensions];
	randomv::sampleDirVector(numDimensions, direction);
	std::vector<double> result (direction, direction+sizeof(direction)/sizeof(double));
	Matrix<double,Dynamic,1> mutationVector(numDimensions,1);
	for(int i=0; i<numDimensions;i++){
		mutationVector(i,0)=result[i];
		//mutationVector.push_back(3.14159);
	}
	mutationVector = cartesianToPolar(mutationVector);
	mutationVector(0,0) =fr;	
	return mutationVector;	
}

//from http://en.wikipedia.org/wiki/N-sphere#Hyperspherical_coordinates
Matrix<double,Dynamic,1> modelFunctions::cartesianToPolar(Matrix<double,Dynamic,1> cart){
	double r = 0;
	int numDim = cart.rows();
	if(numDim==1){
		return cart;
	}
	Matrix<double,Dynamic,1> tempMat(numDim-1,1);
	tempMat.fill(cart(numDim-1,0)*cart(numDim-1,0));
	Matrix<double,Dynamic,1> result(numDim,1);
	result.fill(0);
	for(int i=0; i<cart.rows();i++){
		r+=cart(i,0)*cart(i,0);
	}
	r = pow(r,0.5);
	if(r==0){
		return result;
	}
	result(0,0)=r;
	
	for(int i=cart.rows()-2; i>=0;i--){
		double val = cart(i,0)*cart(i,0);
		for(int j=0;j<=i;j++){
			tempMat(j,0)+=val;
		}
	}
	//cout<<"cartToPolar tempMat"<<endl;
	//cout<<tempMat<<endl;
	for(int j=1;j<cart.rows()-1;j++){
		if(tempMat(j,0)==0){
			result(j,0)=0;
		}else{
			result(j,0)=acos(cart(j-1,0)/pow(tempMat(j-1,0),0.5));
			//cout<<cart(j-1,0)<<"\t"<<pow(tempMat(j-1,0),0.5)<<"\t"<<result(j,0)<<endl;
		}
	}
	
	double denom = pow(pow(cart(numDim-1,0),2)+pow(cart(numDim-2,0),2),0.5);
	if(denom!=0){
		result(numDim-1,0)=acos(cart(numDim-2,0)/denom);
	}
	
	if(cart(numDim-1,0)<0){
		result(numDim-1,0) = 2*M_PI-result(numDim-1,0);
	}
	
	return result;

}

Matrix<double,Dynamic,1> modelFunctions::polarToCartesian(Matrix<double,Dynamic,1> polar){
	double r = polar(0,0);
	int numDim = polar.rows();
	if(numDim==1){
		return polar;
	}
	
	//cout<<"PolarToCart A"<<endl;
	Matrix<double,Dynamic,1> Cartesian(numDim,1);
	//cout<<"PolarToCart B"<<endl;
	for(int i=1;i<polar.rows();i++){
		double res = r;
		//cout<<"PolarToCart C"<<endl;
		for(int j=1;j<=i;j++){
			if(j==i){
				res=res*cos(polar(j,0));
			}else{
				res=res*sin(polar(j,0));
			}
		}
		//cout<<"PolarToCart D"<<endl;
		Cartesian(i-1,0)=res;
		//cout<<"PolarToCart E"<<endl;
	}
	double res = r;
	for(int j=1;j<polar.rows();j++){
		res=res*sin(polar(j,0));
	}
	Cartesian(polar.rows()-1,0)=res;
	//cout<<"PolarToCart F"<<endl;
	return Cartesian;
}

Matrix<double,Dynamic,1> modelFunctions::add(Matrix<double,Dynamic,1> a,Matrix<double,Dynamic,1> b){
  //variables with convenient names
  Matrix<double,Dynamic,1> cartA = polarToCartesian(a);
  Matrix<double,Dynamic,1> cartB = polarToCartesian(b);
  return cartesianToPolar(cartA + cartB);
}

double modelFunctions::tau(double p, double N){
  double t = 0;
  t = 4*N*(-p*log(p)/(1-p));
  return t;
}

double modelFunctions::getStep(){
  double step;
  if (walk == false){
    step = randomv::sampleNormal(0,sigma);
  }else{
    step = randomv::samplePareto(alpha, beta);
  }
  return step;
}

double modelFunctions::vdistance(Matrix<double,Dynamic,1> a, Matrix<double,Dynamic,1> b, int print){
  Matrix<double,Dynamic,1> cartA = polarToCartesian(a);
  Matrix<double,Dynamic,1> cartB = polarToCartesian(b);
  double result;
  if(print!=0){
	cout<<"Starting vdistance computation!"<<endl;
  }
  for(int i=0; i<cartA.rows();i++){
	  if(print!=0){
	  cout<<"vdistance: "<<a(i+1,1)<<" "<<cartA(i+1,1)<<" "<<b(i+1,1)<<" "<<cartB(i+1,1)<<endl;
	  }
	result+=pow(cartA(i+1,1)-cartB(i+1,1),2);
  }
  result = pow(result,0.5);
  if(print!=0){
	cout<<"End Distance is:\t"<<result<<endl;
  }
  return result;

}

int modelFunctions::getStates(double wa, double wb, double wab){
  int state = -1;
  /*
    1 Underdominance
    2 Recesive
    3 Dominant
    4 Overdominance
    5 Codominant
    6 Complete dominance
    7 Complete recesiveness
  */
  if (wab>wa && wab>wb){
    state = 4;
  }else if (wa > wb){
    double h = (wa-wab)/(wa-wb);
    if(h>1){
      state = 1;
    }else if(h>0.5 && h <1){
      state = 3;
    }else if (h == 0.5){
      state = 5;
    }else if(h>0 && h < 0.5){
      state = 2;
    }else{
      state = 6;
    }
  }else if (wb > wa){
    double h = (wb-wab)/(wb-wa);
    if(h>1){
      state = 1;
    }else if(h>0.5 && h <1){
      state = 2;
    }else if (h == 0.5){
      state = 5;
    }else if(h>0 && h < 0.5){
      state = 3;
    }else{
      state = 7; //Complete recesiveness
    }

  }else{
    state = 5; // waa == wbb codominance
  }
  return state;
}

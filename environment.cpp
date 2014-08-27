#include <iostream>
#include <vector>
#include <gsl/gsl_math.h>
#include "randomv.h"
#include "modelFunctions.h"
#include "environment.h"
#include <Eigen/Core>
using namespace Eigen;
using namespace std;

double environment::f;               //! Frequency
double environment::phi;             //! Phase
double environment::fixedTheta;      //! Fixed angle
Matrix<double,Dynamic,1> environment::optimum; //! Coordinates for the optimum
bool environment::periodic;          //! True if the environmental change is periodic
int environment::timeOfChange;       //! Most recent time the optimum changed
bool environment::changing;
double environment::stepProb;

void environment::setOptimum( Matrix<double,Dynamic,1> x ){ optimum = x; }
void environment::setPeriodic( bool x ){ periodic = x; }
Matrix<double,Dynamic,1> environment::getOptimum () { return optimum; }
bool environment::getPeriodic() { return periodic; }
int environment::getTimeOfChange() {return timeOfChange; }
void   environment::setF(double x){f = x;}
double environment::getF(void){return f;}
void   environment::setPhi(double x){phi=x;}
double environment::getPhi(void){return phi;}

void   environment::setFixedTheta(double x){ fixedTheta = x; }
double environment::getFixedTheta() { return fixedTheta; }
bool environment::isChanging(){return changing;}
double environment::getStepProb(){return stepProb;}
void environment::setStepProb(double p){stepProb=p;}
void environment::setChanging(bool c){changing=c;}


void environment::change(int t){ //this is broken with the current implementation as fitness is only calculated when the new allele is created.
  return;
  double r;
  double theta;
  Matrix<double,Dynamic,1> move;
  if(getPeriodic() == true){
    /*
      Simple harmonic motion
      f:   frequency
      phi: phase
      d:   displacement
    */
    double d = 1; //in order b:[0,2]
    theta = fixedTheta; //move in a straight line
    double td  = (double) t;
    r = sin(2*M_PI*f*td + phi) + d;
    move(0,0)=r;
    for(int i=1; i<modelFunctions::getNumDimensions();i++){
    	move(i,0)=theta;
    }
    optimum = move;
  }else{
    r = modelFunctions::getStep();

    move(0,0)=r;
    for(int i=1; i<modelFunctions::getNumDimensions();i++){
		theta = modelFunctions::fTheta();
		move(i,0)=theta;
    }
    optimum = modelFunctions::add(move, optimum);
  }
  timeOfChange = t;
}

double environment::fW(Matrix<double,Dynamic,1> rv){
  /*
    'Gaussian function'
    f(x)=a*exp(-(x-b)^d/(2*c^2))
    b: mean
    a: maximum
    c: shape
    b = xd
    a = 1
    c = 1
    d = 2
  */
  double w = 0;
  //find distance from optimum
 // cerr<<"env.fW";
 // for(int i=0; i<rv.size();i++){
//	cerr<<" "<<rv[i];
 // }
  Matrix<double,Dynamic,1> cart = modelFunctions::polarToCartesian(rv) - modelFunctions::polarToCartesian(optimum);
  
  double x = cart.transpose() * modelFunctions::getSelectionMatrix() * cart;
  //apply fitness function
  if (modelFunctions::getFitnessFunction() == 'g'){
    //gaussian
    w = modelFunctions::getA()*exp((-1*x)/(2*pow(modelFunctions::getC(),2)));
  }else if (modelFunctions::getFitnessFunction() == 'l'){
    //linear
    if(x>=modelFunctions::getMaxR()){
      w = 0;
    }else{
      w = modelFunctions::getA()*(1-(x/modelFunctions::getMaxR()));
    }
  }else{
    cout<<"unknown fitness function: "<<modelFunctions::getFitnessFunction() <<endl;
    exit(1);

  }
// cerr<<" "<<w<<endl;
  return w;
}

double environment::diploidFitness(Matrix<double,Dynamic,1> a, Matrix<double,Dynamic,1> b){
  //todo: finding the diploid phenotype should be separated from finding the diploid fitness
  //add vectords
  Matrix<double,Dynamic,1> newRv = modelFunctions::add(a, b);
//  for(int i=0; i<newRv.size();i++){
//	cout<<"diploidFitness: "<<newRv[i]<<endl;
//  }
  //divide by 2
  newRv(0,0) = newRv(0,0)/2.;
  //find fitness (distance from optimum is computed within fW)
  double wab = fW(newRv);

  /*
   *  double wa = fW(a,m);
   *  double wb = fW(b,m);
   *  double wab = (wa + wb)/2; //h=1/2 always
   */

  return wab;
}

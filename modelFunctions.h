/*
 * modelFunctions.h
 *
 *  Created on: Mar 16, 2013
 *      Author: sandeep
 */

#ifndef MODELFUNCTIONS_H_
#define MODELFUNCTIONS_H_

#include <vector>
#include "randomv.h"
#include <string>
#include <Eigen/Core>
//! Class containing the model functions for organisms and the environment
//! version = 1;
using namespace Eigen;
using namespace std;
class modelFunctions {
	public:
		static void initialize(); //! constructor
		static double fR();
		static double fTheta();
		static Matrix<double,Dynamic,1> add(Matrix<double,Dynamic,1> a, Matrix<double,Dynamic,1> b);
		static double vdistance(Matrix<double,Dynamic,1> a, Matrix<double,Dynamic,1> b, int print);
		static double tau(double p, double N);
		static char   getFitnessFunction(void);
		static void   setFitnessFunction(char x);
		static double getMaxR(void);
		static void   setMaxR(double x);
		static double getA(void);
		static void   setA(double x);
		static double getB(void);
		static void   setB(double x);
		static double getC(void);
		static void   setC(double x);
		static double getD(void);
		static void   setD(double x);
		static vector<double> getParameters(void);
		static void   setParameters(vector<double> x);
		static void   setMutationFunction(char a);
		static char   getMutationFunction(void);
		static void   setWalk(  bool x );
		static void   setSigma( double x );
		static void   setAlpha( double x );
		static void   setBeta(  double x );
		static bool   getWalk(  void );
		static double getSigma( void );
		static double getAlpha( void );
		static double getBeta(  void );
		static double getStep();
		static int    getStates(double waa, double wbb, double wab);
		static Matrix<double,Dynamic,1> polarToCartesian(Matrix<double,Dynamic,1> polar);
		static Matrix<double,Dynamic,1> cartesianToPolar(Matrix<double,Dynamic,1> cart);
		static int getNumberLoci();
		static int getNumberLociOverlap();
		static int getNumDimensions();
		static void setNumberLoci(int i);
		static void setNumDimensions(int i);
		static vector<double> getRecombinationDistances();
		static void setRecombinationDistances(vector<double> recombDists);
		static bool isAllowRecombination();
		static void setAllowRecombination(bool allowR);
		static vector<string> split(const string& s, const string& delim, const bool keep_empty=true);
		static int getNextLocusId(void){
			locusId+=1;
			return locusId;			
		}
		static int getLocusId(void){return locusId;}
		static vector<Matrix<double,Dynamic,Dynamic,RowMajor>> getCovarianceMatrices(void);
		static void setCovarianceMatrices(vector<Matrix<double,Dynamic,Dynamic,RowMajor>> m);
		static Matrix<double,Dynamic,Dynamic,RowMajor> getSelectionMatrix(void);
		static void setSelectionMatrix(Matrix<double,Dynamic,Dynamic,RowMajor> m);
		static Matrix<double,Dynamic,1> getMutationVector(int numDimensions);

		

	//global printing
	//  bool printTable;   //table of balanced pairs
	private:
		static char fitnessFunction; //! linear or Gaussian fitness function
		static double maxR; //! maximum distance from optimum
		static double b; //! mean
		static double a; //! maximum
		static double c; //! shape
		static double d; //! shape of curve
		static vector<double> parameters;
		static char mutationFunction;
		static bool walk; // false: Wiener process (Brownian motion)
		static double sigma; // normal step for Wiener
		static double alpha; // pareto step for Levy flight
		static double beta; //
		static bool allowRecombination;
		static int numberLoci;
		static vector<double> recombinationDistances;
		static int numDimensions;
		static int locusId;
		static Matrix<double,Dynamic,Dynamic,RowMajor> selectionMatrix;
		static vector<Matrix<double,Dynamic,Dynamic,RowMajor>> covarianceMatrices;
};


#endif /* MODELFUNCTIONS_H_ */

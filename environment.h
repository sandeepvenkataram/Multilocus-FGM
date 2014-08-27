/*
 * environment.h
 *
 *  Created on: Mar 16, 2013
 *      Author: sandeep
 */

#ifndef ENVIRONMENT_H_
#define ENVIRONMENT_H_
#include <vector>
#include "randomv.h"
#include "modelFunctions.h"
#include<Eigen/Core>
using namespace Eigen;
using namespace std;
class environment{
	public:
		static void setOptimum( Matrix<double,Dynamic,1> x );
		static void setPeriodic( bool x );
		static Matrix<double,Dynamic,1> getOptimum ();
		static bool getPeriodic();
		static int getTimeOfChange();
		static void change(int t);
		static double fW(Matrix<double,Dynamic,1> rv);
		//find diploidFitness from
		static double diploidFitness(Matrix<double,Dynamic,1> a, Matrix<double,Dynamic,1> b); //find diploidFitness from coordinates
		static bool isChanging();
		static double getStepProb();
		static void setStepProb(double p);
		static void setChanging(bool c);

		//for Harmonic
		static void   setF(double x);
		static double getF(void);
		static void   setPhi(double x);
		static double getPhi(void);
		static void   setFixedTheta(double x);
		static double getFixedTheta();
	private:
		//for the harmonic motion of the mean
		static double f;               //! Frequency
		static double phi;             //! Phase
		static double fixedTheta;      //! Fixed angle
		static Matrix<double,Dynamic,1> optimum; //! Coordinates for the optimum
		static bool periodic;          //! True if the environmental change is periodic
		static int timeOfChange;       //! Most recent time the optimum changed
		static bool changing;
		static double stepProb;
};



#endif /* ENVIRONMENT_H_ */

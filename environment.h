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
		environment(void);
		void setOptimum( Matrix<double,Dynamic,1> x );
		void setPeriodic( bool x );
		Matrix<double,Dynamic,1> getOptimum ();
		bool getPeriodic();
		int getTimeOfChange();
		void change(int t);
		double fW(Matrix<double,Dynamic,1> rv);
		//find diploidFitness from
		double diploidFitness(Matrix<double,Dynamic,1> a, Matrix<double,Dynamic,1> b); //find diploidFitness from coordinates
		bool isChanging();
		double getStepProb();
		void setStepProb(double p);
		void setChanging(bool c);

		//for Harmonic
		void   setF(double x);
		double getF(void);
		void   setPhi(double x);
		double getPhi(void);
		void   setFixedTheta(double x);
		double getFixedTheta();
	private:
		//for the harmonic motion of the mean
		double f;               //! Frequency
		double phi;             //! Phase
		double fixedTheta;      //! Fixed angle
		Matrix<double,Dynamic,1> optimum; //! Coordinates for the optimum
		bool periodic;          //! True if the environmental change is periodic
		int timeOfChange;       //! Most recent time the optimum changed
		bool changing;
		double stepProb;
};



#endif /* ENVIRONMENT_H_ */

/*
 * randomVariables.h
 *
 *  Created on: Mar 15, 2013
 *      Author: sandeep
 */

#ifndef RANDOMV_H_
#define RANDOMV_H_
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

//! random variable class
/*! an object oriented interface for sampling from distributions using gsl
 */
class randomv{
 public:
  static void initialize(void);                          //! default constructor
  static double sampleUniform(void);                     //! sample from the uniform distribution [0,1)
  static int sampleUniformInt(int n);                    //! sample from the uniform distribution [0,n-1)
  static unsigned int sampleBinomial(double p, int N);   //! sample from the binomial distribution
  static double sampleNormal(double mu, double sigma);   //! sample from a normal distribution with mean mu and standard deviation sigma
  static double sampleNormal(double sigma);              //! sample from a normal distribution with mean 0 and standard deviation sigma
  static double sampleExponential(double mu);            //! sample from an exponential distribution with mean mu
  static double sampleGamma(double alpha, double beta);  //! sample from a gamma distribution with parameters alpha and beta
  static void sampleMultinomial(size_t K, unsigned int N, const double p[], unsigned int n[]); //sample from a multinomial distribution
  static double samplePareto(double alpha, double beta); //! sample from a pareto distribution with parameters alpha and beta \fp(x) dx = (a/b) / (x/b)^{a+1} dx

  static double samplePoisson(double lambda);
  static void sampleDirVector(size_t n, double * x);
  static gsl_rng *Gr;
};

#endif /* RANDOMV_H_ */

/*
 * runSimulation.cpp
 *
 *  Created on: Mar 16, 2013
 *      Author: sandeep
 */

#include "population.h"
#include "environment.h"
#include "modelFunctions.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <sstream>
#include <limits>
#include <gsl/gsl_math.h>
#include <gsl/gsl_statistics.h>
#include <Eigen/Core>
using namespace Eigen;
using namespace std;

population * pop;

int              generation;
int              balancedGenerations = 0;
//initial allele parameters


//population parameters
int    n;       //population size
double u;       //mutation rate
int    ploidy;  //ploidy
int    maxTime; //simulation time
char   f;       //new mutation function

//for fitness function
double         a;           //mean
double         b;           //shape
double         c;           //maximum
double         d;           //shape of curve
vector<double> parameters;  //new mutation function parameters

//environmental parameters
double  F;           //frequency of motion
double  P;           //phase
bool    changing;    //changing environment
double  stepProb;    //step probability
bool    periodic;    //periodic change
double  fixedTheta;  // theta for optimum
                     // (if not -PI<fixedTheta<=PI then randomized)
bool     walk;       //random walk (0 Wiener process, 1 Levy flight)
double   sigma;      // sigma for gausian step in Wiener process
double   alpha;      // alpha (scale) for Pareto step in Levy flight
double   beta;       // beta (shape) for Pareto step in Levy flight

//simulation run parameters
int burnIn;          //generations to run before starting to print status
int step;            //print status every step generations
bool allowRecombination;
int numDimensions;
int numberLoci;

void importArguments(char* popFile, char* parFile, char* envFile){
	string line;
	ifstream populationInit(parFile);
	if(populationInit.is_open()){
		while(! populationInit.eof()){
			getline(populationInit, line);
			if (line.size() < 1 || line.at(0) == '#'){
				continue; //skip comments and empty lines
			}
			stringstream sstr(line);
			char parameter;
			sstr>>parameter;
			switch (parameter){
				case 'N':
				 sstr >> n;
				 break;
				case 'u':
				 sstr >> u;
				 break;
				case 'p':
				 sstr >> ploidy;
				 break;
				case 't':
				 sstr >> maxTime;
				 break;
				case 'a':
				 sstr >> a;
				 break;
				case 'b':
				 sstr >> b;
				break;
				case 'c':
				 sstr >> c;
				 break;
				case 'd':
				 sstr >> d;
				 break;
				case 'f':
				 sstr >> f;
				 break;
				case 'P':
				 sstr >>P;
				 break;
				case 'm':
				 {
					double myparameter;
					while(sstr >> myparameter){
					 parameters.push_back(myparameter);
					}
				 }
				 break;
				case 'F':
				 sstr >> F;
				break;
				case 'C':
				 sstr >> changing;
				 break;
				case 'S':
				 sstr >> stepProb;
				 break;
				case 'E':
				 sstr >> periodic;
				 break;
				case 'T':
				 sstr >> fixedTheta;
				 break;
				case 'W':
				 sstr >> walk;
				 break;
				case 's':
				 sstr >> sigma;
				 break;
				case 'A':
				 sstr >> alpha;
				 break;
				case 'B':
				 sstr >> beta;
				break;
				case 'I':
				 sstr >> burnIn;
				 break;
				case 'J':
				 sstr >> step;
				 break;
				case 'L':
				 sstr >> numberLoci;
				 break;
				case 'D':
				 sstr >> numDimensions;
				 break;
				case 'R':
				 sstr >> allowRecombination;
				 break;
				default:
				 {
					cerr<<"unknown parameter "<<parameter<<endl;
					exit(1);
				 }
			}
		}
		populationInit.close();
	}else{
		cerr << "Unable to open file "<<parFile<<endl;
		exit(1);
	}
	
	ifstream environmentInit(envFile);
	vector<double> recombDistances;
	if(environmentInit.is_open()){
		while(! environmentInit.eof()){
			getline(environmentInit, line);
			if (line.size() < 1 || line.at(0) == '#'){
				continue; //skip comments and empty lines
			}
			stringstream sstr(line);
			while(!sstr.eof()){
				double d;
				sstr >> d;
				recombDistances.push_back(d);
			}
		}
	}
	if(recombDistances.size()!=numberLoci-1){
		cerr<<"Number of recombination distances not equal to number of loci - 1!"<<endl;
		exit(1);
	}
	
	
	//read input files
	//environment
	Matrix<double,Dynamic,1> opt(numDimensions,1); // optimum always at origin without loss of generality.

	environment::setOptimum(opt);
	environment::setPeriodic(periodic);
	environment::setF(F);
	environment::setFixedTheta(fixedTheta);
	environment::setPhi(P);
	environment::setChanging(changing);
	environment::setStepProb(stepProb);
	
	
	modelFunctions::setFitnessFunction('g');
	modelFunctions::setA(a);
	modelFunctions::setB(b);
	modelFunctions::setC(c);
	modelFunctions::setD(d);
	modelFunctions::setMutationFunction(f);
	modelFunctions::setParameters(parameters);
	modelFunctions::setWalk(walk);
	modelFunctions::setSigma(sigma);
	modelFunctions::setAlpha(alpha);
	modelFunctions::setBeta(beta);
	modelFunctions::setAllowRecombination(allowRecombination);
	modelFunctions::setRecombinationDistances(recombDistances);
	modelFunctions::setNumberLoci(numberLoci);
	modelFunctions::setNumDimensions(numDimensions);
}

void importSelectionMatrix(char* selFile){
	//cout<<"selection matrix A"<<endl;
	int numDim = modelFunctions::getNumDimensions();
	Matrix<double,Dynamic,Dynamic, RowMajor> selMat(numDim,numDim);
	ifstream selInit(selFile);
	//cout<<"selection matrix B"<<endl;
	if (selInit.is_open()){
		while (! selInit.eof()){
			//cout<<"selection matrix C"<<endl;
			string line;
			getline (selInit,line);
			//cout<<"selection matrix D"<<endl;
			if (line.size() < 1 || line.at(0) == '#'){
				continue; //skip comments and empty lines
			}
			vector<string> tokens;
			tokens = modelFunctions::split(line, " ",true);
			//cout<<"selection matrix E"<<endl;
			for(int k=0; k<tokens.size(); k++){
				selMat(floor(k/numDim),k%numDim) = atof(tokens[k].c_str());
			}
			//cout<<"selection matrix F"<<endl;
			break;
		}
	}
	//cout<<"selection matrix G"<<endl;
	modelFunctions::setSelectionMatrix(selMat);
	//cout<<"selection matrix H"<<endl;
}

void importCovarianceMatrices(char* covFile){
	//cout<<"covariance matrix A"<<endl;
	int numDim = modelFunctions::getNumDimensions();
	//cout<<"covariance matrix B"<<endl;
	Matrix<double,Dynamic,Dynamic, RowMajor> selMat(numDim,numDim);
	vector<Matrix<double,Dynamic,Dynamic, RowMajor>> covMatrices;
	//cout<<"covariance matrix C"<<endl;
	ifstream covInit(covFile);
	if (covInit.is_open()){
		while (! covInit.eof()){
			//cout<<"covariance matrix D"<<endl;
			string line;
			getline (covInit,line);
			//cout<<"covariance matrix E"<<endl;
			if (line.size() < 1 || line.at(0) == '#'){
				continue; //skip comments and empty lines
			}
			//cout<<"covariance matrix F"<<endl;
			vector<string> tokens;
			Matrix<double,Dynamic,Dynamic,RowMajor> covMat(numDim,numDim);
			//cout<<"covariance matrix G"<<endl;
			tokens = modelFunctions::split(line, " ",true);
			//cout<<"covariance matrix H"<<endl;
			for(int k=0; k<tokens.size(); k++){
				covMat(floor(k/numDim),k%numDim) = atof(tokens[k].c_str());
			}
			//cout<<"covariance matrix I"<<endl;
			covMatrices.push_back(covMat);
			//cout<<"covariance matrix J"<<endl;
		}
	}
	//cout<<"covariance matrix K"<<endl;
	modelFunctions::setCovarianceMatrices(covMatrices);
	//cout<<"covariance matrix L"<<endl;
}


void initializePopulation(char* popFile, char* parFile, char* envFile){
	vector<allele *> initial;
	ifstream alleleInit(popFile);
	double totalFreq = 0;
	int alleleID=0;
	//cout<<"Initialize Pop A"<<endl;
	pop = new population(ploidy);
	//cout<<"Initialize Pop B"<<endl;
	
	if (alleleInit.is_open()){
		double totalAlleleFreq=0;
		//cout<<"Initialize Pop C"<<endl;
		while (! alleleInit.eof() ){
			//cout<<"Initialize Pop D"<<endl;
			string line;
			getline (alleleInit,line);
			if (line.size() < 1 || line.at(0) == '#'){
				continue; //skip comments and empty lines
			}
			//cout<<"Initialize Pop E"<<endl;
			//split line to words
			stringstream sstr(line);
			double freq;
			sstr>>freq;
			totalFreq+=freq;
			vector<locus> loci;
			//cout<<"Initialize Pop F"<<endl;
			//use same order as printStatus output
			// so the end of a .table file can be used
			int locusID =0;
			while(!sstr.eof()){ //get thetas for the 
				//cout<<"Initialize Pop G"<<endl;
				string str;
				sstr>>str;
				vector<string> tokens;
				tokens = modelFunctions::split(str, "_",true);
				
				double r = atof(tokens[0].c_str());
				Matrix<double,Dynamic,1> loc(numDimensions,1);
				loc(0,0)=r;
				//cout<<"Initialize Pop H"<<endl;
				for(int k=1; k<tokens.size(); k++){
					double theta = atof(tokens[k].c_str());
					loc(k,0)=theta;
				}
				//cout<<"Initialize Pop I"<<endl;
				if(loc.rows()!=numDimensions){
					cerr << "Wrong number of dimensions for locus "<<str<<endl;
					exit(1);
				}
				//cout<<"Initialize Pop J"<<endl;
				locus l(loc.rows());
				l.setPhenotype(modelFunctions::cartesianToPolar(loc));
				l.setLocusOrder(locusID);
				locusID++;
				loci.push_back(l);
				//cout<<"Initialize Pop K"<<endl;
				
			}
			
			if(loci.size()!=numberLoci){
				cerr<< "Wrong number of loci!"<<endl;
				exit(1);
			}
			//cout<<"Initialize Pop L"<<endl;
			allele * newAllele = new allele(loci);
			//cout<<"Initialize Pop M"<<endl;
			newAllele->setBirthday(0);
			newAllele->setFrequency(freq);
			newAllele->setId(pop->getAlleleCounter());
			pop->setAlleleCounter(pop->getAlleleCounter()+1);
			//cout<<"Initialize Pop N"<<endl;
			initial.push_back(newAllele);
			//cout<<"Initialize Pop O"<<endl;
		}
		alleleInit.close();
	}else{
		cerr << "Unable to open file "<<popFile<<endl;
		exit(1);
	}
	
	if(totalFreq!=1){
		cerr<< "Total Frequency of initial alleles do not add up to one"<<endl;
		exit(1);
	}
	//cout<<"Initialize Pop P"<<endl;
	pop->setAlleles(initial);
	//cout<<"Initialize Pop P1"<<endl;
	pop->setNumGenerations(0);
	pop->setPloidy(ploidy);
	pop->setN(n);
	pop->setU(u);
	//cout<<"Initialize Pop Q"<<endl;
	
	
}


void evolvePopulation(string outPrefix){
	string   out = outPrefix;
	ofstream fp_out(out.append(".table").c_str(), ios::out);
	out = outPrefix;
	ofstream fts(out.append(".ts").c_str(), ios::out);
	out = outPrefix;
	ofstream fse(out.append(".alleles").c_str(),ios::out); //edges
	out = outPrefix;
	ofstream fsl(out.append(".loci").c_str(),ios::out); //edges
	out = outPrefix;
	ofstream fsg(out.append(".genotypes").c_str(),ios::out);
	fp_out<<"NumGenerations Frequency Allele_position Allele_id Allele_fitness mean_fitness age s isFixed"<<endl;
	fts<<"NumGenerations meanFitness numAlleles numPolymorphic Optimum MaxGenW"<<endl;
	fse<<"Id loci numGenerations newPosition newFitness"<<endl;
	fsl<<"Id parentId locusOrder phenotype"<<endl;
	fsg<<"NumGenerations meanFitness genotype frequency fitness"<<endl;
	pop->printInitialAlleles(fse);
	pop->printStatus(fp_out,fts, fsg);
	for(int i=1; i<=maxTime; i++){
		pop->evolve(fse);
		pop->printStatus(fp_out,fts,fsg);
	}
	pop->printAllLoci(fsl);
}

void print1DMatrix(Matrix<double,Dynamic,1> d){
	for(int i=0; i<d.rows(); i++){
		//cout<<d(i,0)<<"_";
	}
}

void printPopInfo(ofstream &fp_out, ofstream &fts, ofstream &fse, ofstream &fsg){
	pop->printStatus(fp_out,fts, fsg);
	//cout<<pop->getAlleleCounter()<<" "<<pop->getMeanFitness()<<" "<<pop->getN()<<" "<<pop->getNumGenerations()<<" "<<pop->getPloidy()<<" "<<pop->getR0()<<" "<<pop->getU()<<endl;
	vector<genotype> genotypes = pop->getGenotypes();
	//cout<<"Printing Genotypes!"<<endl;
	for(int i=0; i<genotypes.size(); i++){
		//cout<<"\tPrinting Genotype "<<i<<endl;
		//cout<<"\t"<<genotypes[i].getFreq()<<" "<<genotypes[i].getW()<<endl;
		//cout<<"\t";
		print1DMatrix(genotypes[i].getPhenotype());
		//cout<<endl;
		vector<allele *> alleles = genotypes[i].getAlleles();
		//cout<<"\tPrinting Alleles in Genotype!"<<endl;
		for(int j=0; j<alleles.size(); j++){
			//cout<<"\t\tPrinting Allele "<<j<<endl;
			allele * al = alleles[j];
			//cout<<"\t\t"<<al->getId()<<" "<<al->getFrequency()<<" "<<al->getW()<<" "<<al->isFixed()<<" "<<al->getBirthday()<<endl;
			//cout<<"\t\t";
			print1DMatrix(al->getPosition());
			//cout<<endl;
			vector<locus> loci = al->getLoci();
			//cout<<"\t\tPrinting Loci in Allele!"<<endl;
			for(int k=0; k<loci.size(); k++){
				//cout<<"\t\t\tPrinting Locus "<<k<<endl;
				locus loc = loci[k];
				//cout<<"\t\t\t"<<loc.getId()<<" "<<loc.getParentId()<<" "<<loc.getLocusOrder()<<" ";
				print1DMatrix(loc.getPhenotype());
				//cout<<endl;
			}
		}
		
	}
	
	
	//cout<<"Printing alleles in population!"<<endl;
	
	unordered_map<string,allele *> alleles = pop->getAlleles();
	int c =0;
	for(unordered_map<string,allele *>::iterator iter = alleles.begin(); iter != alleles.end(); ++iter){
		//cout<<"\t\tPrinting Allele "<<c<<endl;
		allele * al = (*iter).second;
		//cout<<"\t\t"<<al->getId()<<" "<<al->getFrequency()<<" "<<al->getW()<<" "<<al->isFixed()<<" "<<al->getBirthday()<<endl;
		//cout<<"\t\t";
		print1DMatrix(al->getPosition());
		//cout<<endl;
		vector<locus> loci = al->getLoci();
		//cout<<"\t\tPrinting Loci in Allele!"<<endl;
		for(int k=0; k<loci.size(); k++){
			//cout<<"\t\t\tPrinting Locus "<<k<<endl;
			locus &loc = loci[k];
			//cout<<"\t\t\t"<<loc.getId()<<" "<<loc.getParentId()<<" "<<loc.getLocusOrder()<<" ";
			print1DMatrix(loc.getPhenotype());
			//cout<<endl;
		}
		c++;
	}
	
}

int main(int argc, char **argv){ //args = population file, environment file, parameters file, output prefix
	if(argc<7){
		cerr<<"Not Enough Arguments\nCorrect Arguments are: population file, parameters file, environment file, selection Matrix, covariance Matrices, output prefix\n";
		return(1);
	}
	char* popFile = argv[1];
	char* parFile = argv[2];
	char* envFile = argv[3];
	char* selFile = argv[4];
	char* covFile = argv[5];
	string outPrefix = argv[6];
	
	gsl_rng_env_setup();
	gsl_rng_default_seed = (unsigned long)time(0)*getpid();
	const gsl_rng_type *Gt;
	Gt = gsl_rng_default;
	randomv::Gr = gsl_rng_alloc (Gt);
	//cout<<"Starting execution!"<<endl;
	importArguments(popFile, parFile, envFile);
	//cout<<"importing selection matrix"<<endl;
	importSelectionMatrix(selFile);
	//cout<<"importing covariance matrix"<<endl;
	importCovarianceMatrices(covFile);
	//cout<<"Done importing arguments!"<<endl;
	initializePopulation(popFile, parFile, envFile);
	//cout<<"Done initializing population!"<<endl;
	evolvePopulation(outPrefix);
	
	
	/*
	string   out = outPrefix;
	ofstream fp_out(out.append(".table").c_str(), ios::out);
	out = outPrefix;
	ofstream fts(out.append(".ts").c_str(), ios::out);
	out = outPrefix;
	ofstream fse(out.append(".edges").c_str(),ios::out); //edges
	out = outPrefix;
	ofstream fsg(out.append(".genotypes").c_str(),ios::out); //edges
	
	fp_out<<"NumGenerations Frequency Allele_position Allele_id Allele_fitness mean_fitness age s isFixed"<<endl;
	fts<<"NumGenerations meanFitness numAlleles numPolymorphic Optimum MaxGenW"<<endl;
	fse<<"Id numGenerations newPosition newFitness"<<endl;
	
	printPopInfo(fp_out,fts,fse,fsg);
	pop->propagate();
	//cout<<"Done Running Propagate!"<<endl;
	//printPopInfo(fp_out,fts,fse,fsg);
	pop->updateAlleleFrequencies();
	//cout<<"Done Updating Allele Freqs!"<<endl;
	//printPopInfo(fp_out,fts,fse,fsg);
	pop->computeMeanFitness();
	//cout<<"Done Computing Mean Fitness!"<<endl;
	printPopInfo(fp_out,fts,fse,fsg);
	pop->mutate(fse);
	//cout<<"Done Mutating!"<<endl;
	//printPopInfo(fp_out,fts,fse,fsg);
	pop->updateAlleleFrequencies();
	//cout<<"Done Updating Allele Freqs!"<<endl;
	//printPopInfo(fp_out,fts,fse,fsg);
	pop->computeMeanFitness();
	//cout<<"Done Computing Mean Fitness!"<<endl;
	printPopInfo(fp_out,fts,fse,fsg);
	for(int i=0; i<10; i++){
		//cout<<"Starting Generation "<<i<<endl;
		pop->propagate();
		//cout<<"Done Running Propagate!"<<endl;
		printPopInfo(fp_out,fts,fse,fsg);
		pop->updateAlleleFrequencies();
		//cout<<"Done Updating Allele Freqs!"<<endl;
		printPopInfo(fp_out,fts,fse,fsg);
		pop->computeMeanFitness();
		//cout<<"After "<<i<<" generations:"<<endl;
		printPopInfo(fp_out,fts,fse,fsg);
	}
	*/
}

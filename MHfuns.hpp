//
//  MHfuns.hpp
//  OUprocess
//
//  Created by Michael Warner on 5/18/17.
//  Copyright © 2017 Michael Warner. All rights reserved.
//

#ifndef MHfuns_hpp
#define MHfuns_hpp

#include <stdio.h>
#include <iostream>
#include <ctime>
#include <random>
#include <fstream>
#include <string>
using namespace std;

const int nParam = 4,nDataPoint = 5,nNode=2;
const double c1 = 0.282095; //sqrt(1/(2*pi))
void InitializeParameters();
void InitializeFile();
void runML();
void Randomize();
void GenerateData();
void GenProposal();
void TestProposal();
void PrintToFile();
void CalcPosterior();
void CalcPrior();
void CalcLikelihood();
double predDiff();
double determinantOfMatrix(double mat[nNode][nNode], int n);
void getCofactor(double mat[nNode][nNode], double temp[nNode][nNode], int p, int q, int n);
void inverse(double A[nNode][nNode], double inverse[nNode][nNode]);
void adjoint(double A[nNode][nNode],double adj[nNode][nNode]);
double dUnif(double start, double end, double value);
double dNorm(double mean, double sd, double value);
void CalcTrueVals();
void CalcEstimatedVars();
typedef mt19937 MyRNG;  // the Mersenne Twister with a popular choice of parameters
extern double prob1, prob2, prior, likelihood;
extern double TestData[nDataPoint][2];
extern double RealVal[nParam],Prop[nParam],CParam[nParam],stepSize[nParam];
extern int boots,BurnIn;
extern bool Burn,accept;
extern ifstream in;
extern ofstream out;
extern MyRNG rng;



#endif /* MHfuns_hpp */

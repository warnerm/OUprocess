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

const int nParam = 4,nDataPoint = 100,nNode=4;
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
void inverse(double A[nNode][nNode]);
void adjoint(double A[nNode][nNode]);
double dUnif(double start, double end, double value);
double dNorm(double mean, double sd, double value);
void CalcTrueVals();
void CalcEstimatedVars();
extern double prob1, prob2, prior, likelihood;
extern double TestData[nDataPoint][nNode];
extern double RealVal[nParam],Prop[nParam],CParam[nParam],stepSize[nParam];
extern int boots,BurnIn;
extern bool Burn,accept;
extern ifstream in;
extern ofstream out;
extern double adj[nNode][nNode],inv[nNode][nNode];
static std::random_device rd; // random device engine, usually based on /dev/random on UNIX-like systems
// initialize Mersennes' twister using rd to generate the seed
static std::mt19937 rng(rd());




#endif /* MHfuns_hpp */

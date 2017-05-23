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

const int nOptimal = 2,nDataPoint = 100,nTip=5;
const int nParam = nOptimal + 3;
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
void InitializeTipDist();
void CalcLikelihood();
void InitializeIndex();
double CalcVariance(double var,double Par[nParam],int branch);
double CalcExpr(double expr, double Par[nParam],int branch);
double predDiff();
double determinantOfMatrix(double mat[nTip][nTip], int n);
void getCofactor(double mat[nTip][nTip], double temp[nTip][nTip], int p, int q, int n);
void inverse(double A[nTip][nTip]);
void adjoint(double A[nTip][nTip]);
double dUnif(double start, double end, double value);
double dNorm(double mean, double sd, double value);
void CalcTrueVals();
void CalcEstimatedVars();
extern double prob1, prob2, prior, likelihood;
extern double TestData[nDataPoint][nTip];
extern double RealVal[nParam],Prop[nParam],CParam[nParam],stepSize[nParam];
extern int boots,BurnIn;
extern bool Burn,accept;
extern ifstream in;
extern ofstream out;
extern double adj[nTip][nTip],inv[nTip][nTip];
static std::random_device rd; // random device engine, usually based on /dev/random on UNIX-like systems
// initialize Mersennes' twister using rd to generate the seed
static std::mt19937 rng(rd());




#endif /* MHfuns_hpp */

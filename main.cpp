//
//  main.cpp
//  OUprocess
//
//  Created by Michael Warner on 5/18/17.
//  Copyright © 2017 Michael Warner. All rights reserved.
//

#include "MHfuns.hpp"
int boots = 100000,BurnIn = 10000;
bool Burn = true,accept;
ofstream out;
double prob1, prob2;
double TestData[nDataPoint][nNode];
//For the following, define individual variance, phenotypic drift, selection, and optimum expression
double RealVal[nParam] = {10,5,0.5,80},Prop[nParam],CParam[nParam],stepSize[nParam] = {0.2,0.2,0.05,1};
double AncestorExpr = 100, AncestorVar = 25;
double branchTimes[nNode] = {3,2};
double TrueNodeExpr[nNode];
double TrueNodeVar[nNode];
MyRNG rng;


int main(int argc, const char * argv[]) {
    CalcTrueVals();
//    for (int i = 0; i < nNode; i++){
//        cout << TrueNodeVar[i] << endl;
//        cout << TrueNodeExpr[i] << endl;
//    }
    GenerateData(); //While testing utility of approach
//    InitializeParameters();
//    InitializeFile();
//    for (int i = 0; i < nDataPoint; i++){
//        cout << TestData[i][0] << endl;
//        cout << TestData[i][1] << endl;
//    }
//    for (int i=0; i < BurnIn; i++){
//        runML();
//    }
//    Burn = false; //Begin keeping track of values
//    for (int i = 0; i < boots; i++)
//        runML();
    return 0;
}

//Calculate true expected values based on parameter values
void CalcTrueVals(){
    for (int i = 0; i < nNode; i++){
        TrueNodeExpr[i] = AncestorExpr*exp(-RealVal[2]*branchTimes[i]) + RealVal[3]*(1 - exp(-RealVal[2]*branchTimes[i]));
        TrueNodeVar[i] = (RealVal[1]/(2*RealVal[2]))*(1 - exp(-2*RealVal[2]*branchTimes[i])) + AncestorVar*exp(-2*RealVal[2]*branchTimes[i]);
    }
}

//Generate dummy data from random values centered around true values
void GenerateData(){
    for (int i=0 ; i<nNode; i++){
        double var = TrueNodeVar[i] + RealVal[0]; //Includes individual variation
        normal_distribution<double> dis(TrueNodeExpr[i],var); //Initialize standard deviation
        for (int j = 0; j < nDataPoint; j++){
            TestData[j][i] = dis(rng);
        }
    }
}

//Add header to output file
void InitializeFile(){
    out.open("Results.txt");
    out << "tau\tdrift\tselection\toptimal" << endl;
    out.close();
}

//Initialize parameters in middle of distribution
void InitializeParameters(){
    CParam[0] = 0;
    CParam[1] = 0;
    CParam[2] = 15;
}

//Calculate prior probability. Must edit this for each MH algorithm
double CalcPrior(double Par[nParam]){
    double tau = dUnif(0,30,Par[0]);
    double drift = dUnif(0,30,Par[1]);
    double selection = dUnif(0,2,Par[2]);
    double optimal = dUnif(0,1000,Par[3]);
    //If probability can't be calculated, pass back -1, which we'll recognize to reject
    if (tau == -1 || drift == -1 || selection == -1 || optimal == -1) return -1;
    double total = tau + drift + selection + optimal;
    return(total);
}

//Calculate likelihood
double CalcLikelihood(double Par[nParam]){
    double like = 0;
    for (int i = 0; i<nDataPoint; i++){
        double pred = TestData[i][0]*Par[0]+Par[1];
        like = like + dNorm(pred,Par[2],TestData[i][1]);
    }
    return like;
}

//Calculate likelihood of observing the data from a uniform distribution
double dUnif(double start, double end, double value){
    if (value<start || value>end){
        return -1;//Will interpret this value as impossible and skip
    } else {
        double prob = end - start + 1;
        prob = 1/prob;
        return log10(prob);
    }
}

//Calculate likelihood of observing data from Normal distribution
double dNorm(double mean, double sd, double value){
    double prob = (value - mean);
    prob = -prob*prob/(2*sd*sd);
    prob = exp(prob);
    prob = prob*c1/sd;
    return log10(prob);
}



















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
ifstream in;
ofstream out;
double prob1, prob2;
double TestData[nDataPoint][2];
double RealVal[nParam] = {1,5,5},Prop[nParam],CParam[nParam],stepSize[nParam] = {0.2,0.2,0.5};
MyRNG rng;


int main(int argc, const char * argv[]) {
    GenerateData(); //While testing utility of approach
    InitializeParameters();
    InitializeFile();
    for (int i = 0; i < nDataPoint; i++){
        cout << TestData[i][0] << endl;
        cout << TestData[i][1] << endl;
    }
    for (int i=0; i < BurnIn; i++){
        runML();
    }
    Burn = false; //Begin keeping track of values
    for (int i = 0; i < boots; i++)
        runML();
    return 0;
}

//Generate dummy data from random values centered around true values
void GenerateData(){
    normal_distribution<double> dis(0,RealVal[2]); //Initialize standard deviation
    for (int i=0 ; i<nDataPoint; i++){
        TestData[i][0] = i - nDataPoint/2;
        TestData[i][1] = TestData[i][0]*RealVal[0] + RealVal[1] + dis(rng);
    }
}

//Add header to output file
void InitializeFile(){
    out.open("Results.txt");
    out << "a\tb\tsig2\tAccept" << endl;
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
    double a = dUnif(-10,10,Par[0]);
    double b = dUnif(-10,10,Par[1]);
    double sd = dUnif(0,30,Par[2]);
    //If probability can't be calculated, pass back -1, which we'll recognize to reject
    if (a == -1 || sd == -1 || b == -1) return -1;
    double total = a + b + sd;
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



















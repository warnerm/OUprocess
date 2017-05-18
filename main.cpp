//
//  main.cpp
//  OUprocess
//
//  Created by Michael Warner on 5/18/17.
//  Copyright © 2017 Michael Warner. All rights reserved.
//

#include "MHfuns.hpp"
int boots = 100000,BurnIn = 0;
bool Burn = true,accept;
ofstream out;
double prob1, prob2,prior,likelihood;
double TestData[nDataPoint][nNode];
//For the following, define individual variance, phenotypic drift, selection, and optimum expression
double RealVal[nParam] = {10,5,0.5,80},Prop[nParam],CParam[nParam],stepSize[nParam] = {0.2,0.2,0.05,1};
double AncestorExpr = 100, AncestorVar = 25;
double branchTimes[nNode] = {3,2};
double TrueNodeExpr[nNode];
double TrueNodeVar[nNode];
double EstimatedExpr[nNode], EstimatedVar[nNode],Cov[nNode][nNode];
double adj[nNode][nNode],inv[nNode][nNode];
MyRNG rng;


int main(int argc, const char * argv[]) {
    CalcTrueVals();
//    for (int i = 0; i < nNode; i++){
//        cout << TrueNodeVar[i] << endl;
//        cout << TrueNodeExpr[i] << endl;
//    }
    GenerateData(); //While testing utility of approach
    InitializeParameters();
    InitializeFile();
//    for (int i = 0; i < nDataPoint; i++){
//        cout << TestData[i][0] << endl;
//        cout << TestData[i][1] << endl;
//    }
    for (int i=0; i < BurnIn; i++){
        runML();
    }
    Burn = false; //Begin keeping track of values
    for (int i = 0; i < boots; i++)
        runML();
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
    CParam[0] = 15;
    CParam[1] = 15;
    CParam[2] = 0.3;
    CParam[3] = 200;
    std::copy(Prop,Prop+nParam,CParam); //Necessary to calculate Posterior,
    CalcPrior();
    CalcLikelihood();
    CalcPosterior();
}

//Calculate prior probability. Must edit this for each MH algorithm
void CalcPrior(){
    double tau = dUnif(0,30,Prop[0]);
    double drift = dUnif(0,30,Prop[1]);
    double selection = dUnif(0,2,Prop[2]);
    double optimal = dUnif(0,1000,Prop[3]);
    //If probability can't be calculated, pass back -1, which we'll recognize to reject
    if (tau == -1 || drift == -1 || selection == -1 || optimal == -1) prior =  -1;
    else prior = tau + drift + selection + optimal;
}

//Calculate likelihood of multivariate normal distribution
void CalcLikelihood(){
    CalcEstimatedVars();
    likelihood = -(nNode/2)*log10(determinantOfMatrix(Cov,nNode)) -0.5*predDiff(); //predDiff calculates [x - E[x]]'Cov^-1[x - E[x]]
}

//predDiff calculates [x - E[x]]'Cov^-1[x - E[x]], where x is observed expression and E[x] is expected based on parameter values
double predDiff(){
    double temp[nNode][nDataPoint];
    inverse(Cov);
    for (int i = 0; i < nNode; i++){
        for (int j = 0; j < nDataPoint; j++){
            temp[i][j] = (TestData[i][j] - EstimatedExpr[i]); //x - X
        }
    }
    double total = 0;
    for (int j = 0; j < nDataPoint; j++){
        for (int i = 0; i < nNode; i++){
            double temp2 = 0;
            for (int k = 0; k < nNode; k++){
                temp2 = temp2 + temp[k][j]*inv[k][i];
            }
            total = total + temp2 * temp[i][j];
        }
    }
    return(total);
}

//Calculate the estimated mean expression and variance in expression from parameter values
void CalcEstimatedVars(){
    for (int i = 0; i < nNode; i++){
        EstimatedExpr[i] = AncestorExpr*exp(-Prop[2]*branchTimes[i]) + (1 - exp(-Prop[2]*branchTimes[i]))*Prop[3];
        EstimatedVar[i] = (Prop[1]/(2*Prop[2]))*(1 - exp(-2*Prop[2]*branchTimes[i])) + AncestorVar*exp(-2*Prop[2]*branchTimes[i]);
    }
    for (int i = 0; i < nNode; i++){
        for (int j = 0; j < nNode; j++){
            if ( i == j ) Cov[i][j] = EstimatedVar[i]; //Covariance with itself is variance
            else if (i > j) Cov[i][j] = Cov[j][i]; //Already calculated it, so save a bit of time
            else Cov[i][j] = AncestorVar*exp(-Prop[2]*(branchTimes[i]+branchTimes[j]));
        }
    }
    
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



















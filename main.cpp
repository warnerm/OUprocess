//
//  main.cpp
//  OUprocess
//
//  Created by Michael Warner on 5/18/17.
//  Copyright © 2017 Michael Warner. All rights reserved.
//

#include "MHfuns.hpp"
int boots = 50000000,BurnIn = 10000000;
const int nNode = nTip*2 - 1;
int ancestor[nNode];
bool Burn = true,accept;
ofstream out;
double prob1, prob2,prior,likelihood;
double TestData[nDataPoint][nTip];
//For the following, define individual variance, phenotypic drift, selection, and optimum expression
double RealVal[nParam] = {3,5,2,100,70},Prop[nParam],CParam[nParam],stepSize[nParam] = {0.1,0.1,0.1,0.5,0.5};
int optimalIndex[nNode - 1] = {3,3,3,3,3,3,3,4}; //Define optimal expression levels for different branches
double branchTimes[nNode - 1] = {0.5,0.4,2.9,3,0.2,0.2,0.5,2};
double TrueTipExpr[nNode];
double TrueTipVar[nNode];
double EstimatedExpr[nNode], EstimatedVar[nNode],Cov[nTip][nTip];
double adj[nTip][nTip],inv[nTip][nTip];
int MutAncestor[nTip][nTip];
double TipDist[nTip][nTip];

int main(int argc, const char * argv[]) {
    InitializeParameters();
    InitializeIndex();
    InitializeTipDist();
    CalcTrueVals();
    GenerateData(); //While testing utility of approach
    InitializeFile();
    for (int i = 0; i < nTip; i++){
        for (int j = 0; j < nTip; j++){
           
        }

    }
    for (int i=0; i < BurnIn; i++){
        runML();
//        cout << likelihood << endl;
//        cout << "pre" << predDiff() << endl;
//        for (int j = 0; j < nTip; j++){
//            cout << EstimatedVar[j] << endl;
//            cout << "expr" << EstimatedExpr[j] << endl;
//        }
    }
    Burn = false; //Begin keeping track of values
    for (int i = 0; i < boots; i++){
        runML();
    }
    return 0;
}

//Initialize key for which nodes are ancestor, matrix of divergence times
void InitializeIndex(){
    ancestor[0] = 0;
    for (int i = 1; i < (nTip - 1); i++){ //interior nodes
        ancestor[i] = i - 1;
    }
    for (int i = (nTip - 1); i < (nNode - 1); i++){ //Tips, before last one
        ancestor[i] = i - nTip + 1;
    }
    ancestor[nNode - 1] = nTip - 2;//Last tip shares ancestor with tip before it
}

//Initialize matrix of branch distances between tips
void InitializeTipDist(){
    for (int i = 0; i < nTip; i++){
        for (int j = 0; j < nTip; j++){
            if (i > j) MutAncestor[i][j] = MutAncestor[j][i];
            else MutAncestor[i][j] = i;
            if (i > j) TipDist[i][j] = TipDist[j][i];
            else {
                TipDist[i][j] = branchTimes[nTip + i - 1] + branchTimes[nTip + j - 1];
                int temp = ancestor[j + nTip - 1];
                int ancI = ancestor[i + nTip - 1];
                //Only j's ancestor can be greater because we restrict to cases where i < j
                if (temp > ancI){
                    TipDist[i][j] = TipDist[i][j] + branchTimes[temp]; //Adds branch length of interior node to nearest interior node
                    temp--;
                }
            }
        }
    }
}

//Calculate true expected values based on parameter values
void CalcTrueVals(){
    TrueTipExpr[0] = CParam[3]; //Set as optimum of 1st branch
    TrueTipVar[0] = 0;
    //After entering the ancestral variance, each successive node based on ancestor
    for (int i = 1; i < nNode; i++){
        TrueTipExpr[i] = CalcExpr(TrueTipExpr[ancestor[i]],RealVal,i);
        TrueTipVar[i] = CalcVariance(TrueTipVar[ancestor[i]],RealVal,i);
    }
}

//Calculate expected variance based on variance of immediately ancestral node
double CalcVariance(double var,double Par[nParam],int branch){
    double Variance = (Par[1]/(2*Par[2]))*(1 - exp(-2*Par[2]*branchTimes[branch])) + var*exp(-2*Par[2]*branchTimes[branch]);
    return(Variance);
}

//Calculate expected expression based on variance of immediately ancestral node
double CalcExpr(double expr, double Par[nParam],int branch){
    double Expression = expr*exp(-Par[2]*branchTimes[branch]) + Par[optimalIndex[branch]]*(1 - exp(-Par[2]*branchTimes[branch]));
    return(Expression);
}

//Generate dummy data from random values centered around true values
void GenerateData(){
    for (int i=0 ; i<nTip; i++){ //Don't generate values for interior nodes
        double var = TrueTipVar[i + nTip - 1] + RealVal[0]; //Includes individual variation
        normal_distribution<double> dis(TrueTipExpr[i + nTip - 1],var); //Initialize standard deviation
        for (int j = 0; j < nDataPoint; j++){
            TestData[j][i] = dis(rng);
        }
    }
}

//Add header to output file
void InitializeFile(){
    out.open("Results2.txt");
    out << "tau\tdrift\tselection\t";
    for (int i = 0; i < nOptimal; i++){
        out << "optimal" << i << "\t";
    }
    out << "accept" << endl;
    out.close();
}

//Initialize parameters in middle of distribution
void InitializeParameters(){
    for (int i = 0; i < nParam; i++){
        CParam[i] = 5;
    }
    CParam[3] = 50;
    CParam[4] = 100;
    std::copy(CParam,CParam+nParam,Prop); //Necessary to calculate Posterior,
    CalcPrior();
    CalcLikelihood();
    CalcPosterior();
    prob2 = prob1;
}

//Calculate prior probability. Must edit this for each MH algorithm
void CalcPrior(){
    prior = 0;
    for (int i = 0; i < nParam; i++){
        double prob = dUnif(0,1000,Prop[i]);
        if (prob == - 1){
            //If probability can't be calculated, pass back -1, which we'll recognize to reject
            prior = -1;
            break;
        } else {
            prior = prior + prob;
        }
    }
}

//Calculate likelihood of multivariate normal distribution
void CalcLikelihood(){
    CalcEstimatedVars();
    likelihood = -(nDataPoint/2)*log10(determinantOfMatrix(Cov,nTip)) -0.5*predDiff(); //predDiff calculates [x - E[x]]'Cov^-1[x - E[x]]
}

//predDiff calculates [x - E[x]]'Cov^-1[x - E[x]], where x is observed expression and E[x] is expected based on parameter values
double predDiff(){
    double temp[nTip];
    double total = 0;
    inverse(Cov);
    for (int j = 0; j < nDataPoint; j++){
        //Deviation of observed values from estimator
        for (int i = 0; i < nTip; i++){
            temp[i] = (TestData[j][i] - EstimatedExpr[nTip + i - 1]); //x - X
        }
        for (int i = 0; i < nTip; i++){
            for (int k = 0; k< nTip; k++){
                total = total + temp[i]*inv[i][k]*temp[k];
            }
        }
    }
    return(total);
}

//Calculate the estimated mean expression and variance in expression from parameter values
void CalcEstimatedVars(){
    EstimatedExpr[0] = Prop[3];
    EstimatedVar[0] = 0;
    for (int i = 1; i < nNode; i++){ //skip root of tree
        EstimatedExpr[i] = CalcExpr(EstimatedExpr[ancestor[i]], Prop, i);
        EstimatedVar[i] = CalcVariance(EstimatedVar[ancestor[i]], Prop, i);
    }
    for (int i = 0; i < nTip; i++){
        for (int j = 0; j < nTip; j++){
            if ( i == j ) Cov[i][j] = EstimatedVar[i + nTip - 1] + Prop[0]; //Covariance with itself is variance; add individual variance here
            else if (i > j) Cov[i][j] = Cov[j][i]; //Already calculated it, so save a bit of time
            else Cov[i][j] = EstimatedVar[MutAncestor[i][j]]*exp(-Prop[2]*TipDist[i][j]);
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



















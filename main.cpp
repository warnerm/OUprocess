//
//  main.cpp
//  OUprocess
//
//  Created by Michael Warner on 5/18/17.
//  Copyright © 2017 Michael Warner. All rights reserved.
//

#include "MHfuns.hpp"
int boots = 1000000,BurnIn = 100000;
const int nNode = nTip*2;
int ancestor[nNode];
bool Burn = true,accept;
ofstream out;
double prob1, prob2,prior,likelihood;
double TestData[nDataPoint][nTip];
//For the following, define individual variance, phenotypic drift, selection, and optimum expression
double RealVal[nParam] = {0,10,2,80},Prop[nParam],CParam[nParam],stepSize[nParam] = {0.5,0.5,0.25,1};
double AncestorExpr = 250, AncestorVar = 10;
double branchTimes[nNode] = {0,0.1,0.3,3,0.4,1};
double TrueTipExpr[nNode];
double TrueTipVar[nNode];
double EstimatedExpr[nNode], EstimatedVar[nNode],Cov[nTip][nTip];
double adj[nTip][nTip],inv[nTip][nTip];
int MutAncestor[nTip][nTip];
double TipDist[nTip][nTip];

int main(int argc, const char * argv[]) {
    TrueTipExpr[0] = EstimatedExpr[0] = AncestorExpr;
    TrueTipVar[0] = EstimatedVar[0] = AncestorVar;
    branchTimes[0] = 0;
    InitializeIndex();
    InitializeTipDist();
    CalcTrueVals();
    GenerateData(); //While testing utility of approach
    InitializeParameters();
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
    for (int i = 1; i < nTip; i++){ //interior nodes
        ancestor[i] = i - 1;
    }
    for (int i = nTip; i < (nNode - 1); i++){ //Tips, before last one
        ancestor[i] = i - (nTip - 1);
    }
    ancestor[nNode - 1] = nTip - 1;//Last tip shares ancestor with tip before it
}

//Initialize matrix of branch distances between tips
void InitializeTipDist(){
    for (int i = 0; i < nTip; i++){
        for (int j = 0; j < nTip; j++){
            if (i > j) MutAncestor[i][j] = MutAncestor[j][i];
            else MutAncestor[i][j] = i + 1;
            if (i > j) TipDist[i][j] = TipDist[j][i];
            else {
                TipDist[i][j] = branchTimes[nTip+i] + branchTimes[nTip+j];
                int temp = ancestor[j+nTip];
                int ancI = ancestor[i + nTip];
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
    double Expression = expr*exp(-Par[2]*branchTimes[branch]) + Par[3]*(1 - exp(-Par[2]*branchTimes[branch]));
    return(Expression);
}

//Generate dummy data from random values centered around true values
void GenerateData(){
    for (int i=0 ; i<nTip; i++){ //Don't generate values for interior nodes
        double var = TrueTipVar[i+nTip] + RealVal[0]; //Includes individual variation
        normal_distribution<double> dis(TrueTipExpr[i+nTip],var); //Initialize standard deviation
        for (int j = 0; j < nDataPoint; j++){
            TestData[j][i] = dis(rng);
        }
    }
}

//Add header to output file
void InitializeFile(){
    out.open("Results2.txt");
    out << "tau\tdrift\tselection\toptimal\taccept" << endl;
    out.close();
}

//Initialize parameters in middle of distribution
void InitializeParameters(){
    CParam[0] = 3;
    CParam[1] = 3;
    CParam[2] = 0.1;
    CParam[3] = 120;
    std::copy(CParam,CParam+nParam,Prop); //Necessary to calculate Posterior,
    CalcPrior();
    CalcLikelihood();
    CalcPosterior();
    prob2 = prob1;
}

//Calculate prior probability. Must edit this for each MH algorithm
void CalcPrior(){
    double tau = dUnif(0,30,Prop[0]);
    double drift = dUnif(0,1000,Prop[1]);
    double selection = dUnif(0,30,Prop[2]);
    double optimal = dUnif(0,1000,Prop[3]);
    //If probability can't be calculated, pass back -1, which we'll recognize to reject
    if (tau == -1 || drift == -1 || selection == -1 || optimal == -1) {
        prior =  -1;
    }
    else prior = tau + drift + selection + optimal;
}

//Calculate likelihood of multivariate normal distribution
void CalcLikelihood(){
    CalcEstimatedVars();
    likelihood = -(nDataPoint/2)*log10(determinantOfMatrix(Cov,nTip)) -0.5*predDiff(); //predDiff calculates [x - E[x]]'Cov^-1[x - E[x]]
    cout << "det" << log10(determinantOfMatrix(Cov,nTip)) << endl;
    cout << "pred" << predDiff() << endl;
}

//predDiff calculates [x - E[x]]'Cov^-1[x - E[x]], where x is observed expression and E[x] is expected based on parameter values
double predDiff(){
    double temp[nTip];
    double total = 0;
    inverse(Cov);
    for (int j = 0; j < nDataPoint; j++){
        //Deviation of observed values from estimator
        for (int i = 0; i < nTip; i++){
            temp[i] = (TestData[j][i] - EstimatedExpr[nTip + i]); //x - X
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
    for (int i = 1; i < nNode; i++){ //skip root of tree
        EstimatedExpr[i] = CalcExpr(EstimatedExpr[ancestor[i]], Prop, i);
        EstimatedVar[i] = CalcVariance(EstimatedVar[ancestor[i]], Prop, i);
    }
    for (int i = 0; i < nTip; i++){
        for (int j = 0; j < nTip; j++){
            if ( i == j ) Cov[i][j] = EstimatedVar[i + nTip]; //Covariance with itself is variance; add individual variance here
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



















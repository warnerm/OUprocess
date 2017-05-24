//
//  main.cpp
//  OUprocess
//
//  Created by Michael Warner on 5/18/17.
//  Copyright © 2017 Michael Warner. All rights reserved.
//

#include "MHfuns.hpp"
int boots = 1000000,BurnIn = 100000,nRun,nAccept,nValley=100,meanBranchLength=2;
double timeStep = 0.001;
const int nNode = nTip*2 - 1;
int ancestor[nNode];
bool Burn = true,accept;
ofstream out;
double prob1, prob2,prior,likelihood;
double TestData[nDataPoint][nTip];
//For the following, define individual variance, phenotypic drift, selection, and optimum expression
double RealVal[nParam],Prop[nParam],CParam[nParam],stepSize[nParam];
int optimalIndex[nNode - 1]; //Define optimal expression levels for different branches
double branchTimes[nNode - 1];
double EstimatedExpr[nNode], EstimatedVar[nNode],Cov[nTip][nTip],SimExpr[nNode];
double adj[nTip][nTip],inv[nTip][nTip];
int MutAncestor[nTip][nTip];
double TipDist[nTip][nTip];
double TrueTipExpr[nNode],TrueTipVar[nNode];

int main(int argc, const char * argv[]) {
    ConstructBranches();
    InitializeOptimal();
    InitializeIndex();
    InitializeTipDist();
    GenerateTrueVals();
    InitializeFile();
    SimulateData();
    BurnInML();
    for (nRun = 0; nRun < boots; nRun++){
        runML();
    }
    return 0;
}

//Construct phylogeny randomly, drawing branch times from a uniform distribution, mean = meanBranchLength
void ConstructBranches(){
    uniform_real_distribution<double> dis(0,2*meanBranchLength); //Initialize standard deviation
    for (int i = 0; i < (nNode - 1); i++){
        branchTimes[i] = dis(rng); //For now, make everything with the same optimal
    }
}

//Construct array of optimal expression
void InitializeOptimal(){
    for (int i = 0; i < (nNode - 1); i++){
        optimalIndex[i] = 3; //For now, make everything with the same optimal
    }
}

//Burn-in period
void BurnInML(){
    //Want to do burn in and make sure that we've left "valleys" of vanishingly low likelihood
    while (true){
        nAccept = 0;
        InitializeParameters();
        for (int i=0; i < BurnIn; i++){
            runML();
        }
        if (nAccept > nValley){
            break;
        }
    }
    Burn = false; //Begin keeping track of values
    for (int i = 0; i < nParam; i++){
        stepSize[i] = 0.25; //Decrease step size for ML estimation
    }
}

//Generate true values based on command line input, drawing parameters from distributions
void GenerateTrueVals(){
    uniform_real_distribution<double> dis(0,1000);
    for (int i = 0; i < nParam; i++){
        RealVal[i] = dis(rng);
        stepSize[i] = 2; //Start with big step sizes for the burn-in period
    }
}

//Simulate data under the O-U model with the given parameter states
void SimulateData(){
    SimExpr[0] = RealVal[3];
    for (int i = 1; i < nNode; i++){
        double expr = SimExpr[ancestor[i]];
        double time = 0;
        //Stochastic differential equation governing OU process
        while (time < branchTimes[i - 1]){
            expr = expr + RealVal[2]*(RealVal[optimalIndex[i - 1]] - expr)*timeStep + Drift(RealVal[1]);
            time = time + timeStep;
        }
        SimExpr[i] = expr;
    }
    GenerateSimData();
}

//Takes Simulated mean expression values and draws random values based on individual variance
void GenerateSimData(){
    for (int i = 0; i < nTip; i++){
        normal_distribution<double> dis(SimExpr[i+nTip - 1],RealVal[0]); //Initialize standard deviation
        for (int j = 0; j < nDataPoint; j++){
            TestData[j][i] = dis(rng);
            //cout << i << j << TestData[j][i] << endl;
        }
    }
}
//Calculates change per time-step due to drift
double Drift(double sig2){
    normal_distribution<double> dis(0,sig2*timeStep); //Initialize standard deviation
    return(dis(rng));
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

//Calculate expected variance based on variance of immediately ancestral node
double CalcVariance(double var,double Par[nParam],int branch){
    double Variance = (Par[1]/(2*Par[2]))*(1 - exp(-2*Par[2]*branchTimes[branch])) + var*exp(-2*Par[2]*branchTimes[branch]);
    return(Variance);
}

//Calculate expected expression based on variance of immediately ancestral node
double CalcExpr(double expr, double Par[nParam],int branch){
    double Expression = expr*exp(-Par[2]*branchTimes[branch - 1]) + Par[optimalIndex[branch - 1]]*(1 - exp(-Par[2]*branchTimes[branch - 1]));
    return(Expression);
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
    uniform_real_distribution<double> dis(0,1000);
    for (int i = 0; i < nParam; i++){
        CParam[i] = dis(rng);
    }
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



















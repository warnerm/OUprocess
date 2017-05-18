//
//  MHfuns.cpp
//  OUprocess
//
//  Created by Michael Warner on 5/18/17.
//  Copyright © 2017 Michael Warner. All rights reserved.
//

#include "MHfuns.hpp"
extern int boots;

//Use Maximum Likelihood to update parameters
void runML(){
    GenProposal();
    prob1 = Posterior(Prop);
    if (prob1 == -1){
        accept = false;
    } else {
        prob2 = Posterior(CParam);
        TestProposal();
        if (!Burn) PrintToFile(); //Don't keep burn-in values
    }
}

//Generate proposal parameters based on values, priors
void GenProposal(){
    for (int i = 0; i < nParam; i++){
        normal_distribution<double> dis(CParam[i],stepSize[i]);
        Prop[i] = dis(rng);
    }
}

//Calculate posterior probability
double Posterior(double Par[nParam]){
    double prior = CalcPrior(Par);
    if (prior == -1) return -1;
    double likelihood = CalcLikelihood(Par);
    return(prior + likelihood);
}

//Accept or reject proposal
void TestProposal(){
    uniform_real_distribution<double> dis(0,1);
    double rNum = dis(rng);
    prob1 = exp(prob1 - prob2);
    if (rNum < prob1){
        std::copy(Prop,Prop+nParam,CParam);
        accept = true;
    } else {
        accept = false;
    }
}

//Print current parameter states, whether proposal was accepted
void PrintToFile(){
    ofstream out;
    out.open("Results.txt",ios::app);
    for (int i = 0; i < (nParam); i++){
        out << CParam[i] << '\t';
    }
    out << accept << endl;
    out.close();
}

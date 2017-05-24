//
//  MHfuns.cpp
//  OUprocess
//
//  Created by Michael Warner on 5/18/17.
//  Copyright © 2017 Michael Warner. All rights reserved.
//

#include "MHfuns.hpp"

//Use Maximum Likelihood to update parameters
void runML(){
    GenProposal();
    CalcPosterior();
    if (prob1 == -1){
        accept = false;
    } else {
        TestProposal();
        if (!Burn & (nRun % 10000 == 0)) PrintToFile(); //Don't keep burn-in values
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
void CalcPosterior(){
    CalcPrior();
    if (prior == -1) prob1 = -1;
    else {
        CalcLikelihood();
        prob1 = prior + likelihood;
    }
}

//Accept or reject proposal
void TestProposal(){
    uniform_real_distribution<double> dis(0,1);
    double rNum = dis(rng);
    double prob = exp(prob1 - prob2);
    if (rNum < prob){
        std::copy(Prop,Prop+nParam,CParam);
        accept = true;
        nAccept++;
        prob2 = prob1; //Next round, the posterior for the current parameters is equal to the current proposal
    } else {
        accept = false;
    }
}

//Print current parameter states, whether proposal was accepted
void PrintToFile(){
    ofstream out;
    out.open(file.c_str(),ios::app);
    for (int i = 0; i < nParam; i++){
        out << CParam[i] << '\t';
    }
    out << accept << endl;
    out.close();
}

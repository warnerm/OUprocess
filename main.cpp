//
//  main.cpp
//  OUprocess
//
//  Created by Michael Warner on 5/18/17.
//  Copyright © 2017 Michael Warner. All rights reserved.
//

#include <iostream>
#include <ctime>
#include <random>
#include <fstream>
#include <string>
using namespace std;

void InitializeParameters();
void InitializeFile();
void runML();
void Randomize();
void GenerateData();
typedef mt19937 MyRNG;  // the Mersenne Twister with a popular choice of parameters
int boots = 1000,BurnIn = 100;
bool Burn = true;
struct Params {
    double a,b,sd;
};
ifstream in;
ofstream out;
Params RealVal,Prop,CParam;



int main(int argc, const char * argv[]) {
    GenerateData(); //While testing utility of approach
    InitializeParameters();
    InitializeFile();
    RealVal = {1,5,10};
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
    //Generate dummy data from random values centered around true values
}

//Add header to output file
void InitializeFile(){
    //Add header to output file
    out.open("Results.txt");
    out << "a\tb\tsig\n" << endl;
    out.close();
}

//Initialize parameters based on priors
void InitializeParameters(){
    //Initialize parameters
}

//Use Maximum Likelihood to update parameters
void runML(){
    //Run ML
}







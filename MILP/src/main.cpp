#include <iostream>
#include <fstream>
#include "modelAES.hpp"


using namespace std;

int main(int argc, char const *argv[]){

    unsigned R = stoi(argv[1]);	
    int version = stoi(argv[2]); 
    modelAES(R, version);
    return 0;
}

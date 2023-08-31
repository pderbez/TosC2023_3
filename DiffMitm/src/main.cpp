#include <iostream>
#include <fstream>
#include "modelAES.hpp"


using namespace std;

int main(int argc, char const *argv[]){

    unsigned Round = stoi(argv[1]);	
    
    unsigned min_rin = 0;
    unsigned max_rin = 3;
    
    unsigned min_rout = 0;
    unsigned max_rout = 3;
    
    if (argc > 2) {
    	min_rin = stoi(argv[2]);
    	max_rin = min_rin + 1;
    	
    	min_rout = stoi(argv[3]);
    	max_rout = min_rout + 1;
    }
    

 
    for(unsigned rin = min_rin; rin < max_rin; rin++){
        for(unsigned rout = min_rout; rout < max_rout; rout++){
            unsigned rmiddle = Round - rin - rout;

            cout << rin << " - " << rmiddle << " - " << rout << endl;
            modelAES(rin, rmiddle, rout, 128);
        }
    }
    return 0;
}

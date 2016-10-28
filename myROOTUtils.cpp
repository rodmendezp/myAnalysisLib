#include "myROOTUtils.h"

using namespace std;

Bool_t fexists(string fName){
    if(FILE *file = fopen(fName.c_str(), "r")){
        fclose(file);
        return true;
    }
    else
        return false;
}

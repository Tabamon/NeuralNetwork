#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

//|(0 0)|(1 0)|
//|-----+-----|
//|(0 1)|(1 1)|


int main(void){
int a;
double square [2][2];
double result[4];
std::ofstream fileOut;

fileOut.open("training_data.txt");

for (unsigned int i = 0; i<6000; i++){
    a=rand()%2;
    square[0][0]=1.0*a;
    a=rand()%2;
    square[0][1]=1.0*a;
    a=rand()%2;
    square[1][0]=1.0*a;
    a=rand()%2;
    square[1][1]=1.0*a;

    if (square[0][0]==1.0&&square[1][1]==1.0&&square[0][1]==0.0&&square[1][0]==0.0){
        result [0] = 1.0;
        result [1] = 0.0;
        result [2] = 0.0;
        result [3] = 0.0;
    }else if (square[0][0]==0.0&&square[1][1]==0.0&&square[0][1]==1.0&&square[1][0]==1.0){
        result [0] = 1.0;
        result [1] = 0.0;
        result [2] = 0.0;
        result [3] = 0.0;
    }else if (square[0][0]==1.0&&square[1][1]==0.0&&square[0][1]==1.0&&square[1][0]==0.0){
        result [0] = 0.0;
        result [1] = 1.0;
        result [2] = 0.0;
        result [3] = 0.0;
    }else if (square[0][0]==1.0&&square[1][1]==0.0&&square[0][1]==0.0&&square[1][0]==1.0){
        result [0] = 0.0;
        result [1] = 1.0;
        result [2] = 0.0;
        result [3] = 0.0;
    }else if (square[0][0]==0.0&&square[1][1]==1.0&&square[0][1]==1.0&&square[1][0]==0.0){
        result [0] = 0.0;
        result [1] = 0.0;
        result [2] = 1.0;
        result [3] = 0.0;
    }else if (square[0][0]==1.0&&square[1][1]==0.0&&square[0][1]==0.0&&square[1][0]==1.0){
        result [0] = 0.0;
        result [1] = 0.0;
        result [2] = 1.0;
        result [3] = 0.0;
    }else{
        result [0] = 0.0;
        result [1] = 0.0;
        result [2] = 0.0;
        result [3] = 1.0;
    }

    fileOut << square[0][0] << " " << square[1][0] << " " << square[0][1] << " " << square[1][1] << "\n";
    fileOut << result [0] << " " << result [1] << " " << result [2] << " " << result [3] << "\n";

    
}

}
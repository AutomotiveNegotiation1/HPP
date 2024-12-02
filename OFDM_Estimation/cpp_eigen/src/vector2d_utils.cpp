#include <vector2d_utils.h>
using namespace std; 

vector<vector<complex<double>>> transpose_conjugate(vector<vector<complex<double>>> orgV){
    
    int orgCol = orgV[0].size() ; 
    int orgRow = orgV.size(); 

    int newRow = orgCol; 
    int newCol = orgRow; 

    vector<vector<complex<double>>> transposeV(newRow, vector<complex<double>>(newCol,0));  

    for(int i=0; i<newRow; i++){
        for(int j=0; j<newCol; j++){
            transposeV[i][j] = std::conj(orgV[j][i]); 
        }
    }
    return transposeV; 
}


vector<vector<complex<double>>> matrix_multiplication(vector<vector<complex<double>>> vA, vector<vector<complex<double>>> vB){
    
    int A_Row = vA.size(); 
    int A_Col = vA[0].size() ; 

    int B_Row = vB.size(); 
    int B_Col = vB[0].size() ; 

    vector<vector<complex<double>>> retV(A_Row, vector<complex<double>>(B_Col,0)); 

    for(int i=0; i<A_Row; i++){
        for(int j=0; j<B_Col; j++){
            complex<double> temp = (0,0); 
            for(int k=0; k<A_Col; k++){
                temp += vA[i][k] * vB[k][j];
            }
            retV[i][j] = temp; 
        }    
    }

    return retV; 
}


void printout_vector(vector<vector<complex<double>>> V){
    for (int i=0; i<V.size(); i++){
        for(int j=0; j<V[i].size(); j++){
            cout << V[i][j] << " "; 
        }
        cout << std::endl; 
    }
    cout << std::endl; 
}

Eigen::MatrixXcd Vector2Eigen(vector<vector<complex<double>>> inputV){

    int Row = inputV.size(); 
    int Col = inputV[0].size() ; 

    Eigen::MatrixXcd retEigen(Row, Col);

    for (int i=0; i<Row; i++){
        for(int j=0; j<Col; j++){
            retEigen(i,j) = inputV[i][j]; 
        }
    }
    return retEigen; 
}
#include <vector> 
#include <complex> 
#include <eigen3/Eigen/Dense>
#include <iostream> 

using namespace std; 

vector<vector<complex<double>>> transpose_conjugate(vector<vector<complex<double>>> orgV); 
void printout_vector(vector<vector<complex<double>>> V);

Eigen::MatrixXcd Vector2Eigen(vector<vector<complex<double>>> inputV);


vector<vector<complex<double>>> matrix_multiplication(vector<vector<complex<double>>> vA, vector<vector<complex<double>>> vB);
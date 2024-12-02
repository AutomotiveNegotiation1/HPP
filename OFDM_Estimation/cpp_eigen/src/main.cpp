#include <vector> 
#include <complex> 
#include <iostream> 
#include <eigen3/Eigen/Dense>
#include <iomanip> 
#include <load_data.h> 
#include <vector2d_utils.h>
#include <ESPRIT_algorithm.h>

using namespace std; 

int main(){

    //Load Data
    vector<vector<complex<double>>> Azi_cov = load_AziCov(); 
    vector<vector<complex<double>>> Ele_cov = load_EleCov(); 

    //Input check 
    std::cout << "input Azi_cov" << std::endl; 
    printout_vector(Azi_cov);
    
    // Vector to Eigen 
    Eigen::MatrixXcd Azi_cov_eigen = Vector2Eigen(Azi_cov);
    
    // Apply SVD using Eigen Library 
    Eigen::JacobiSVD<Eigen::MatrixXcd> svd(Azi_cov_eigen, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::MatrixXcd V = svd.matrixV(); 

    // Apply ESPRIT Algorithm 

    double EstAzi_esprit = get_EstAzi_esprit(V); 

    // Print the answer
    std::cout << "EstAzi_esprit: " << EstAzi_esprit << std::endl;    

    return 0; 

}
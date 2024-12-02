#include <ESPRIT_algorithm.h> 


double get_EstAzi_esprit(Eigen::MatrixXcd V){

    int row_num = V.rows();

    Eigen::MatrixXcd U1_azi = V.topRows(row_num-1); // 3 = (Row num) -1 
    Eigen::MatrixXcd U2_azi = V.bottomRows(row_num-1); // 3 = (Row num) -1 

    Eigen::MatrixXcd temp = U1_azi.colPivHouseholderQr().solve(U2_azi); 
    Eigen::VectorXcd eivals= temp.eigenvalues();

    vector<double> eival_args ;  

    for (int i=0; i<row_num; i++){
        std::complex tmp = eivals(i); 
        
        eival_args.push_back(std::arg(tmp)); 
    }

    double EstAzi_esprit = (-1) * std::asin(eival_args[row_num-1] / M_PI); // 3 : last number of Row (0,1,2,3)

    return EstAzi_esprit;
}

        
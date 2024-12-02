#include <load_data.h>
#include <vector2d_utils.h>

vector<vector<complex<double>>> load_AziData(){
    
    vector<complex<double>> V;
    //test 1
    
    V.push_back(std::complex(-0.189451219294928 , 0.099964611969799));
    V.push_back(std::complex(0.054021408912110, -0.207283321851416));
    V.push_back(std::complex(0.116561748086865, 0.179716574142569));
    V.push_back(std::complex(-0.211294689121474, -0.035202876127765));
    
    // Test 2 
    /*
    V.push_back(std::complex(-0.164692234606777, 0.136971369266842)); 
    V.push_back(std::complex(0.044645204341920, -0.209502968637381)); 
    V.push_back(std::complex(0.094542544956705, 0.192214451412262)); 
    V.push_back(std::complex(-0.193197110088571, -0.092517916059711)); 
    */

    //Test 3 
    /*
    V.push_back(std::complex(-0.161682595127447, 0.140511304070033)); 
    V.push_back(std::complex(0.030922869633668, -0.211963355965663)); 
    V.push_back(std::complex(0.114800365717460, 0.180846797510337)); 
    V.push_back(std::complex(-0.204971951936703, -0.062218864164056)); 
    */
    
    vector<vector<complex<double>>> retV;
    retV.push_back(V); 
   
    return retV; 
}

vector<vector<complex<double>>> load_EleData(){
    
    vector<vector<complex<double>>> retV;
    
    //test1 
    
    vector<complex<double>> tempV = {std::complex(-0.189451219294928, 0.099964611969799)};
    retV.push_back(tempV);
    tempV = {std::complex(-0.212358339349166, -0.028082447319843)};
    retV.push_back(tempV);
    tempV = {std::complex(-0.156954587870231, -0.145773610386425)};
    retV.push_back(tempV);
    tempV = {std::complex(-0.043671072380530, -0.209708191484600)};
    

    //Test 2 
    /*
    vector<complex<double>> tempV = {std::complex(-0.164692234606777, 0.136971369266842)};
    retV.push_back(tempV);
    tempV = {std::complex(-0.205566560266340, 0.060225222613724)};
    retV.push_back(tempV);
    tempV = {std::complex(-0.212568512433695, -0.026444577144160)};
    retV.push_back(tempV);
    tempV = {std::complex(-0.184544339551780, -0.108756953239802)};
    */
    
    //Test 3
    /*
    vector<complex<double>> tempV = {std::complex(-0.161682595127447, 0.140511304070033)};
    retV.push_back(tempV);
    tempV = {std::complex(-0.202761198719023, 0.069083894162331)};
    retV.push_back(tempV);
    tempV = {std::complex(-0.213838247380748, -0.012565512155795)};
    retV.push_back(tempV);
    tempV = {std::complex(-0.193274725926358, -0.092355662829620)};
    */
    retV.push_back(tempV);
    
    return retV; 
}


vector<vector<complex<double>>> load_AziCov(){
    
    vector<vector<complex<double>>> Azi_data = load_AziData(); // 1  by 4
    vector<vector<complex<double>>> Azi_t = transpose_conjugate(Azi_data);
    vector<vector<complex<double>>> Azi_cov = matrix_multiplication(Azi_t, Azi_data); 
    
    return Azi_cov; 
} 

vector<vector<complex<double>>> load_EleCov(){
    
    vector<vector<complex<double>>> Ele_data = load_EleData(); // 1  by 4
    vector<vector<complex<double>>> Ele_t = transpose_conjugate(Ele_data);
    vector<vector<complex<double>>> Ele_cov = matrix_multiplication(Ele_data, Ele_t); 
    
    return Ele_cov; 
} 

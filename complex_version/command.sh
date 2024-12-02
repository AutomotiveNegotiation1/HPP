
rm -rf *.o 
g++ -c jacobi_cplx_evd.o jacobi_cplx_evd.cpp 
g++ -c mat_utils_cplx.o mat_utils_cplx.cpp 
g++ -c jacobi_cplx_svd.o jacobi_cplx_svd.cpp 
g++ -c jacobi_cplx_svd_test.o jacobi_cplx_svd_test.cpp 
g++ jacobi_cplx_evd.o jacobi_cplx_svd_test.o mat_utils_cplx.o jacobi_cplx_svd.o -o output 
./output 

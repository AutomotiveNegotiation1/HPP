# rm -rf *.o 
# g++ -Wall -Wextra -c jacobi_evd.o jacobi_evd.c
# g++ -Wall -Wextra -c jacobi_svd.o jacobi_svd.c
# g++ -Wall -Wextra -c mat_utils.o mat_utils.c
# g++ -Wall -Wextra -c jacobi_svd_test.o jacobi_svd_test.cpp

# #g++ -Wall -Wextra -c jacobi_cplx_evd.o jacobi_cplx_evd.cpp

# g++ -Wall -Wextra jacobi_svd_test.o jacobi_evd.o jacobi_svd.o mat_utils.o -o output
# ./output


rm -rf *.o 
g++ -c jacobi_cplx_evd.o jacobi_cplx_evd.cpp 
g++ -c mat_utils_cplx.o mat_utils_cplx.cpp 
g++ -c jacobi_cplx_svd.o jacobi_cplx_svd.cpp 
g++ -c jacobi_cplx_svd_test.o jacobi_cplx_svd_test.cpp 
g++ jacobi_cplx_evd.o jacobi_cplx_svd_test.o mat_utils_cplx.o jacobi_cplx_svd.o -o output 
./output 

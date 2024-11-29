g++ -Wall -Wextra -c jacobi_evd.o jacobi_evd.c
g++ -Wall -Wextra -c jacobi_svd.o jacobi_svd.c
g++ -Wall -Wextra -c mat_utils.o mat_utils.c
g++ -Wall -Wextra -c jacobi_svd_test.o jacobi_svd_test.cpp
g++ -Wall -Wextra jacobi_svd_test.o jacobi_evd.o jacobi_svd.o mat_utils.o -o output
./output

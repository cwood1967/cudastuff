extern "C"
#include <math.h>

__global__ void test(int n, float *a, float*b, float *res)
{
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i<n)
    {
        res[i] = exp(-a[i]*b[i]);
    }
}        
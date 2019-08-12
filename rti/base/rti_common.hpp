#ifndef RTI_COMMON_HPP
#define RTI_COMMON_HPP

/// \file
///
/// A header including CUDA related headers and functions

#if defined(__CUDACC__)

#   include <cuda_runtime.h>
#   include <cublas.h>
#   include <curand_kernel.h>
//#include "helper/helper_cuda.h"
//#include "helper/helper_math.h"
#   include <curand.h>
#   include <stdio.h>

#   define CUDA_HOST_DEVICE __host__ __device__
#   define CUDA_HOST        __host__ 
#   define CUDA_DEVICE      __device__

//#   define rti_sqrt(__CUDACC__)
//
//sqrtf : for float, sqrtg: for double
//T norm()const{return sqrtf(x*x+y*y+z*z);}
//#else
//T norm()const{return std::sqrt(x*x+y*y+z*z);}
//#endif


#else

#   define CUDA_HOST_DEVICE 
#   define CUDA_HOST        
#   define CUDA_DEVICE      


#endif

#endif


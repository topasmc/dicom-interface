#ifndef RTI_CONST1D_H
#define RTI_CONST1D_H

/// \file
/// 
/// Distribution functions (meta-header file for all distributions)
 
#include <random>
#include <functional>
#include <queue>
#include <iostream>
#include <fstream>
#include <random>
#include <chrono>
#include <functional>
#include <array>

#include <rti/base/rti_vec.hpp>
#include <rti/base/rti_matrix.hpp>
#include <rti/base/distributions/rti_pdfMd.hpp>


namespace rti{

/// \class const_1d
///
/// 1-dimensional const pdf.
/// \tparam T type of return value
/// \note sigam values are ignored
template<typename T>
class const_1d : public pdf_Md<T,1> {

public:
    
    /// Constructor
    CUDA_HOST_DEVICE
    const_1d(
        std::array<T,1>& m, 
        std::array<T,1>& s) 
        : pdf_Md<T,1>(m,s)
    {;}
    
    /// Constructor
    CUDA_HOST_DEVICE
    const_1d(
        const std::array<T,1>& m, 
        const std::array<T,1> &s)
        : pdf_Md<T,1>(m,s) 
    {;}

    /// Returns mean_ 
    CUDA_HOST_DEVICE
    virtual 
    std::array<T,1>
    operator()(void){
        return pdf_Md<T,1>::mean_;
    };

};

}
#endif
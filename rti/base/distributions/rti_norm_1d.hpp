
#ifndef RTI_NORM_1D_H
#define RTI_NORM_1D_H

/// \file
/// 
/// Distribution functions for 1D normal distribution


#include <rti/base/distributions/rti_pdfMd.hpp>

namespace rti{
/// \class norm_1d
///
/// 1-dimensional normal pdf.
/// \tparam T type of return value
/// \note in future, the randome engine will be relocated.
template<typename T>
class norm_1d : public pdf_Md<T,1> {
public:

    std::default_random_engine  gen_ ; ///< C++ random engine
    std::normal_distribution<T> func_; ///< C++ distribution function

    /// Constructor to copy mean and sigma and initialize normal_distributions
    CUDA_HOST_DEVICE
    norm_1d(
        std::array<T,1>& m_,
        std::array<T,1>& s_)
    : pdf_Md<T,1>(m_,s_)
    {
        #if !defined(__CUDACC__)
        func_ = std::normal_distribution<T>(m_[0], s_[0]);
        #endif
    }
    
    /// Constructor to copy mean and sigma
    CUDA_HOST_DEVICE
    norm_1d(
        const std::array<T,1>& m_,
        const std::array<T,1> &s_)
    : pdf_Md<T,1>(m_,s_)
    {
        #if !defined(__CUDACC__)
        func_ = std::normal_distribution<T>(m_[0], s_[0]);
        #endif
    }

    /// Returns value sampled from normal distribution
    CUDA_HOST_DEVICE
    virtual 
    std::array<T,1>
    operator()(void){
        return {func_(gen_)};
    };

};

}

#endif

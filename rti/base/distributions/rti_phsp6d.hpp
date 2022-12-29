#ifndef RTI_PHSP6D_H
#define RTI_PHSP6D_H

/// \file
/// 
/// Distribution functions (meta-header file for all distributions)

#include <rti/base/distributions/rti_pdfMd.hpp>

namespace rti{
/// \class phsp_6d
///
/// 6-dimensional uniform pdf for phase-space variables for a spot
//  phase-space variables are position (x,y,z) and direction (x',y',z')
/// \tparam T type of return value
template<typename T>
class phsp_6d : public pdf_Md<T,6> {

/// correlations for x-x' and y-y' respectively
std::array<T,2> rho_; ///< For X,Y

public:
    
    /// Random engine and distribution function
    std::default_random_engine  gen_ ;
    std::normal_distribution<T> func_;

    /// Constructor to initializes mean, sigma, rho, and random engine
    /// \param m[0,1,2]: mean spot-position  of x, y, z
    /// \param m[3,4,5]: mean spot-direction of x', y', z'.
    /// \param s[0,1,2]: std  spot-position of x,y,z -> spot-size
    /// \param s[3,4,5]: std  spot-direction of x,y,z. s[5] is ignored but calculated internally
    /// \param r[0] : x, xp correlation
    /// \param r[1] : y, yp correlation
    /// \note seed setup needs to be done by public method
    /// gen_.seed(from outside, topas or UI); //Ideally set from TOPAS?
    CUDA_HOST_DEVICE
    phsp_6d(
        std::array<T,6>& m_,
        std::array<T,6>& s_,
        std::array<T,2>& r_)
    : pdf_Md<T,6>(m_,s_), rho_(r_)
    {
        #if !defined(__CUDACC__)
        gen_.seed(std::chrono::system_clock::now().time_since_epoch().count());
        func_ = std::normal_distribution<T>(0, 1); 
        #endif
    }
    
    /// Constructor to initializes mean, sigma, rho, and random engine
    CUDA_HOST_DEVICE
    phsp_6d(
        const std::array<T,6>& m_,
        const std::array<T,6>& s_,
        const std::array<T,6>& r_)
    : pdf_Md<T,6>(m_,s_), rho_(r_)
    {
        #if !defined(__CUDACC__)
        gen_.seed(std::chrono::system_clock::now().time_since_epoch().count());
        func_ = std::normal_distribution<T>(0, 1); 
        #endif
    }

    /// Sample 6 phase-space variables and returns
    CUDA_HOST_DEVICE
    virtual 
    std::array<T,6>
    operator()(void)
    {
        std::array<T,6> phsp = pdf_Md<T,6>::mean_;
        T Ux = func_(gen_); T Vx = func_(gen_);
        T Uy = func_(gen_); T Vy = func_(gen_);
        T Uz = func_(gen_); //T Vz = func_(gen_);
        phsp[0] += pdf_Md<T,6>::sigma_[0]* Ux ; 
        phsp[1] += pdf_Md<T,6>::sigma_[1]* Uy ; 
        phsp[2] += pdf_Md<T,6>::sigma_[2]* Uz ; 

        phsp[3] += pdf_Md<T,6>::sigma_[3] * (rho_[0] * Ux + Vx * std::sqrt(1.0-rho_[0]*rho_[0]));
        phsp[4] += pdf_Md<T,6>::sigma_[4] * (rho_[1] * Uy + Vy * std::sqrt(1.0-rho_[1]*rho_[1]));
        phsp[5]  =  -1.0*std::sqrt( 1.0 - phsp[3]*phsp[3]- phsp[4]*phsp[4] );

        return phsp;
    };

};

}
#endif

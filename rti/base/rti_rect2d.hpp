#ifndef RTI_RECT2D_H
#define RTI_RECT2D_H

/// \file
///
/// 2D rectangle grid

#include <cmath>
#include <array>
#include <valarray>
#include <sstream>
#include <fstream>
#include <map>

#if defined(__CUDACC__)
#include <cuda_runtime.h>
#include <cublas.h>
#include <curand_kernel.h>
//#include "helper/helper_cuda.h"
//#include "helper/helper_math.h"
#include <curand.h>
#include <stdio.h>
#endif

#include <rti/base/rti_vec.hpp>
#include <rti/base/rti_line.hpp>

namespace rti{

/// \class rect2d
/// \tparam T for grid coordinate (float or double) a type of x, y, and z positions
template<typename T>
class rect2d {

public:
    /// Enumerate
    /// Following lines can be accessible as
    /// rti::rect2d mycell(...);
    /// mycell[XP]; will return const l11to10
    /// mycell[XP].p1 ; will return first point
    /// mycell[XP].p2 ; will return first point
    typedef enum 
    {   //corner a,b,c,d 
        XP=0,  ///< X plus 
        YP=1,  ///< Y plus
        XM=2,  ///< X minus
        YM=3,  ///< Y minus
        NOSIDE=4 
    } cell_side;

    rti::line2d<T> l00to01 ; ///< xminus
    rti::line2d<T> l01to11 ; ///< yplus
    rti::line2d<T> l11to10 ; ///< xplus
    rti::line2d<T> l10to00 ; ///< yminus

    rti::vec2<T> center   ; 
    rti::vec2<T> dxy      ;
    rti::vec2<T> d05_xy   ;

    CUDA_HOST_DEVICE
    rect2d(
        rti::vec2<T> c, 
        rti::vec2<T> d)
    : center(c) , 
      dxy(d), 
      d05_xy(d*0.5)
    {
        update_lines();
    }

    CUDA_HOST_DEVICE
    rect2d(const rect2d& ref ) 
    :   l00to01(ref.l00to01), 
        l01to11(ref.l01to11), 
        l11to10(ref.l11to10),
        l10to00(ref.l10to00),
        center(ref.center),
        dxy(ref.dxy),
        d05_xy(ref.d05_xy)
    {;}

    CUDA_HOST_DEVICE
    ~rect2d(){;}

    CUDA_HOST_DEVICE
    void
    change_position( rti::vec2<T>& c )
    {
        center = c;
        update_lines();
    }

    CUDA_HOST_DEVICE
    void
    change_size( rti::vec2<T>& d )
    {
        dxy = d;
        update_lines();
    }

    CUDA_HOST_DEVICE
    void
    update_lines( void )
    {

        l00to01.p1.x = center.x - d05_xy.x; l00to01.p1.y = center.y - d05_xy.y;
        l00to01.p2.x = center.x - d05_xy.x; l00to01.p2.y = center.y + d05_xy.y;

        l01to11.p1.x = center.x - d05_xy.x; l01to11.p1.y = center.y + d05_xy.y;
        l01to11.p2.x = center.x + d05_xy.x; l01to11.p2.y = center.y + d05_xy.y;
        
        l11to10.p1.x = center.x + d05_xy.x; l11to10.p1.y = center.y + d05_xy.y;
        l11to10.p2.x = center.x + d05_xy.x; l11to10.p2.y = center.y - d05_xy.y;

        l10to00.p1.x = center.x + d05_xy.x; l10to00.p1.y = center.y - d05_xy.y;
        l10to00.p2.x = center.x - d05_xy.x; l10to00.p2.y = center.y - d05_xy.y;

    }
    

    CUDA_HOST_DEVICE
    T
    intersect_alpha(
        rti::line2d<T>& line,
        cell_side which_side) const
    {
        return (*this)[which_side].intersect_alpha(line);
    }
    
    CUDA_HOST_DEVICE
    std::array<T, 4>
    intersect_alpha(
        rti::line2d<T>& line)
    const
    {
        std::array<T, 4> v;
        v[XP] = this->intersect_alpha(line, XP);
        v[XM] = this->intersect_alpha(line, XM);
        v[YP] = this->intersect_alpha(line, YP);
        v[YM] = this->intersect_alpha(line, YM);
        return v;
    }

    CUDA_HOST_DEVICE
    const rti::line2d<T>&
    operator[](cell_side which_side)
    const
    {
        //check which_side > NOSIDE
        switch( which_side ){
        case XM:
            return l00to01;
        case XP:
            return l11to10;
        case YM:
            return l10to00;
        case YP:
            return l01to11;
        default:
            throw; //should be an error
        }
    }

};

}

#endif
#ifndef RTI_UTILS_H
#define RTI_UTILS_H


/// \file
/// Header file containing general functions useful for several classes
/// \see http://www.martinbroadhurst.com/how-to-trim-a-stdstring.html
/// \see https://stackoverflow.com/questions/1798112/removing-leading-and-trailing-spaces-from-a-string


#include <string>
#include <map> 
#include <array>
#include <vector>
#include <tuple>

namespace rti{

    ///< Remove white-space on right
    inline std::string trim_right_copy(const std::string& s,const std::string& delimiters = " \f\n\r\t\v\0\\" )
    {
        return s.substr( 0, s.find_last_not_of( delimiters ) + 1 );
    }

    ///< Remove white-space on left
    inline std::string trim_left_copy(const std::string& s,const std::string& delimiters = " \f\n\r\t\v\0\\" )
    {
    return s.substr( s.find_first_not_of( delimiters ) );
    }

    ///< Remove white-space on left and right
    inline std::string trim_copy(const std::string& s,const std::string& delimiters = " \f\n\r\t\v\0\\" )
    {
    return trim_left_copy( trim_right_copy( s, delimiters ), delimiters );
    }

    /// Interpolates values in map
    /// \tparam T is type of values 
    /// \tparam S is size of columns
    /// \param db map of dataset
    /// \param x position value where interpolated value is calculated at
    /// \param y_col column index of map for y-values to be used for interpolation
    template<typename T, size_t S>
    inline 
    T interp_linear(
        const std::map<T, std::array<T,S>>& db,
        const T x,
        const size_t y_col=0)
    {
        auto it_up = db.upper_bound(x);
        if(it_up == db.end()){
            return ((--it_up)->second)[y_col]; 
        }
        if(it_up == db.begin()){
            return (it_up->second)[y_col]; 
        }
        T x1 = it_up->first; 
        T y1 = (it_up->second)[y_col];
        auto it_down = --it_up;
        T x0 = it_down->first; T y0 = (it_down->second)[y_col];
        return y0 + (x-x0)*(y1-y0)/(x1 - x0);
    }

    /// Interpolates values in vector
    /// \tparam T is type of values 
    /// \tparam S is size of columns
    /// \param db is vector for dataset
    /// \param x position value where interpolated value is calculated at
    /// \param x_col column index of map for x-values to be used for interpolation
    /// \param y_col column index of map for y-values to be used for interpolation
    template<typename T, size_t S>
    inline 
    T interp_linear(
        const std::vector<std::array<T,S>>& db,
        const T x,
        const size_t x_col =0,
        const size_t y_col=1)
    {
        if( x <= db[0][x_col] ) return db[0][y_col];
        
        for(size_t i = 1 ; i < db.size() - 2 ; ++i){
            if( x <= db[i][x_col] ){

                T x0 = db[i-1][x_col];
                T x1 = db[i][x_col];
                T y0 = db[i-1][y_col];
                T y1 = db[i][y_col];
                
                return y0 + (x-x0)*(y1-y0)/(x1 - x0);

            }
        }
        
        return db[db.size()-1][y_col];

    }



    
    /// Interpolate table lambda function
    /// \param vector_X The array of x coordinates
    /// \param vector_Y The array of y coordinates
    /// \param x the ordinate to evaluate
    /// \param npoints the number of coordinates in the table
    /// \return the y value corresponding to the x ordinate
    /// \note from gpmc code
    inline 
    float 
    TableInterpolation(
        float* const vector_X, 
        float* const vector_Y, 
        const float x, 
        const int npoints)
    {
        float result;
        int order = 4; // order of the poly
        // Allocate enough space for any table we'd like to read.
        float lambda[npoints]; // lambda[npoints] error: variable length arrays are not supported in OpenCL
        // check order of interpolation
        if (order > npoints) order = npoints;
        // if x is ouside the vector_X[] interval
        if (x <= vector_X[0]) return result = vector_Y[0];
        if (x >= vector_X[npoints-1]) return result = vector_Y[npoints-1];
        // loop to find j so that x[j-1] < x < x[j]
        int j=0;
        while (j < npoints)
        {
            if (vector_X[j] >= x) break;
            j++;
        }
        // shift j to correspond to (npoint-1)th interpolation
        j = j - order/2;
        // if j is ouside of the range [0, ... npoints-1]
        if (j < 0) j=0;
        if (j+order > npoints ) j=npoints-order;
        result = 0.0;
        for (int is = j; is < j+order; is++)
        {
            lambda[is] = 1.0;
            for (int il=j; il < j+order; il++)
            {
                if(il != is) lambda[is] = lambda[is]*(x-vector_X[il])/(vector_X[is]-vector_X[il]);
            }
            result += vector_Y[is]*lambda[is];
        }
        return result;
    }
}

#endif

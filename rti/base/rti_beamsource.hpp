#ifndef RTI_BEAMSOURCE_H
#define RTI_BEAMSOURCE_H

/// \file
///
/// A beamsource is a collection of beamlets and provides an interface for sampling.

#include <rti/base/rti_beamlet.hpp>

namespace rti{

/// \class beamsource
///
/// \tparam T for types for return values by the distributions
/// \note 
template <typename T> 
class beamsource {
public:

    /// A beamlet container. 
    /// Each element is a tuple of beamlet, number of histories, and accumlated histories
    std::vector<std::tuple<rti::beamlet<T>, size_t, size_t>> beamlets_;
    
    /// A lookup table to map a history to beamlet id.
    /// note: To run GPU, this std::map needs to be replaced GPU compatible data structure.
    std::map<std::size_t, std::uint32_t>  cdf2beamlet_; 

    /// A lookup table to map a time to beamlet id.
    //std::map<std::double_t, std::uint32_t> deliverytime2beamlet;

    /// Default constructor
    beamsource(){beamlets_.clear();}

    /// Add a beamlet to internal containers
    /// \param b a beamlet
    /// \param h number of histories
    /// \param p coordinate transformation
    void 
    append_beamlet(
        rti::beamlet<T> b, 
        size_t h,
        rti::coordinate_transform<T> p)
    {
        b.set_coordinate_transform(p);
        this->append_beamlet(b,h);
    }

    /// Add a beamlet to internal containers.
    /// \param b a beamlet
    /// \param h number of histories
    void 
    append_beamlet(
        rti::beamlet<T> b, 
        size_t h)
    {
        size_t acc = total_histories() + h ; 
        cdf2beamlet_.insert( std::make_pair(acc, this->total_beamlets()) ); 
        beamlets_.push_back( std::make_tuple(b,h,acc) );
    }

    /// Returns size of beamlets.
    /// \return total number of beamlets
    std::size_t
    total_beamlets()
    {return beamlets_.size();}

    /// Returns total number of histories
    /// \return total number of histories
    std::size_t
    total_histories()
    {    
        return this->total_beamlets() == 0 ? 
            0 : std::get<2>(beamlets_.back());
    }

    /// Returns a tuple for given beamlet id.
    /// \return a tuple of beamlet, histories, and accumulated histories
    const 
    std::tuple<rti::beamlet<T>, size_t, size_t>& 
    operator[]
    (unsigned int i)
    {return beamlets_[i];}

    /// Returns a beamlet of a history
    /// \return a beamlet reference (const)
    const
    rti::beamlet<T>&
    operator()
    (size_t h)
    {
        size_t beamlet_id = cdf2beamlet_.upper_bound(h)->second;
        return std::get<0>(beamlets_[beamlet_id]);
    }


};

}

#endif

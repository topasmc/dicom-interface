#ifndef RTI_VERTEX_HPP
#define RTI_VERTEX_HPP

#include <rti/base/rti_vec.hpp>

namespace rti{

///< Particle properties at a position
///< Physics will propose next vertex
///< Geometry will propose next vertex
///< Step limiter will propose next vertex
///< Magnetic field will propose next vertex
///< Vertex doesn't include particle type
template<typename T>
struct vertex_t
{
	 T          ke   ;  //< kinetic energy
	 vec3<T>    pos  ;  //< position
	 vec3<T>    dir  ;  //< direction

	 CUDA_HOST_DEVICE
	 vertex_t<T>&
	 operator=(const vertex_t<T>& rhs)
	 {
		 ke  = rhs.ke ;
		 pos = rhs.pos;
		 dir = rhs.dir;
		 return *this;
	 }
};

}
#endif

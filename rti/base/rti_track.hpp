#ifndef RTI_TRACK_HPP
#define RTI_TRACK_HPP

#include <rti/base/rti_vertex.hpp>
#include <rti/base/rti_node.hpp>

namespace rti{
	
///< Particle type
typedef
enum{
 PHOTON   = 0,
 ELECTRON = 1,
 PROTON   = 2,
 NEUTRON  = 3
} particle_t;
//PDG[PHOTON] =
//PDG[PROTON] = 2212 //

///< Process type
typedef
enum{
	 BEGIN    ,
	 MAX_STEP ,
	 BOUNDARY ,
	 CSDA     ,
	 D_ION    ,
	 PP_E     ,
	 PO_E     ,
	 PO_I
} process_t;
///< status of particle
///< CREATED : just created
///< ALIVE   : under tracking
///< STOPPED : no further tracking is required due to limitation, e.g., cut
///< KILLED  : no further tracking is required due to energy 0 or exit world boundary?
typedef
enum{
	 CREATED    = 0,  ///< created 
	 //	 PRE_STEP   = 1,  ///< tracking. next vertex is calculated
	 //POST_STEP  = 2,
	 STOPPED    = 3  ///< stopped by physics process
} status_t;

	
///< Track class
/// pointer to geometry where current track is placed
/// pre-/post- vertex points: vtx0 vtx1   
   
template<typename R>
class track_t {
public:
	status_t       status   ; ///< particle status	
	process_t      process  ; ///< id of physics process that limit the step,
	                          ///< -1: geometry, 0: CSDA, 1: delta, 2: pp-e, 3: po-e, 4: po-i 
	bool           primary  ; 
	particle_t     particle ; ///< particle type, 

	vertex_t<R>    vtx0     ; ///< vertex pre
	vertex_t<R>    vtx1     ; ///< vertex post
	
	R              dE       ; ///< total energy deposit between vtx0 to vtx1

	node_t<R>*     c_node = nullptr ;  //current node
	
	///< auxiliary information for current node
	uint32_t        cnb  =0 ;  ///copy number of geometry of a node, 0 is default
	vec3<ijk_t>     cell ;
	cell_side       side = rti::NONE_XYZ_PLANE  ;

	///< Defaut constructor
	CUDA_HOST_DEVICE
	track_t() : status(CREATED), process(BEGIN), primary(true), dE(0)
	{;}

	///< Constructor
	CUDA_HOST_DEVICE
	track_t(const vertex_t<R>& v) :  status(CREATED), process(BEGIN), primary(true), dE(0)
	{ vtx0 = v;
	  vtx1 = v;}

	///< Constructor
	CUDA_HOST_DEVICE
	track_t
	(status_t  s,
	 process_t p,
	 bool   is_p,
	 particle_t t,
	 vertex_t<R> v0,
	 vertex_t<R> v1,
	 const R& dE)
		: status(s), process(p), primary(is_p),
		  particle(t), vtx0(v0), vtx1(v1), dE(0)
	{;}

	///< copy constructor
	CUDA_HOST_DEVICE
	track_t
	(const track_t& rhs)
	{
		status   = rhs.status   ;
		process  = rhs.process  ;
		particle = rhs.particle ;
		primary  = rhs.primary  ;
		vtx0     = rhs.vtx0     ;
		vtx1     = rhs.vtx1     ;
		dE       = rhs.dE       ;
		c_node   = rhs.c_node   ;
		cnb      = rhs.cnb      ;
		cell     = rhs.cell     ;
		side     = rhs.side     ;
	}

	///< Destructor
	CUDA_HOST_DEVICE
	~track_t()
	{;}


	///< Deposit energy 
	CUDA_HOST_DEVICE
	void
	deposit(R e)
	{
		dE += e ;
	}

	CUDA_HOST_DEVICE
	bool
	is_stopped()
	{ return status == STOPPED ;}

	///< Update vertex point for given R
	CUDA_HOST_DEVICE
	void
	shorten_step
	(R ratio) //0 < ratio < 1
	{
		vtx1.pos  = vtx0.pos + (vtx1.pos - vtx0.pos) * ratio ;
	}
	
	
	///< 
	CUDA_HOST_DEVICE
	void
	update_post_vertex_direction
	(const R& theta,
	 const R& phi,
	 const vec3<R>& dir_z = vec3<R>(0,0,-1))
	{

		rti::mat3x3<R> m_local(0, theta, phi);
		rti::vec3<R>   d_local = m_local * dir_z ; d_local.normalize();
		rti::mat3x3<R> m_global(dir_z, vtx1.dir);
		vtx1.dir  = m_global * d_local ;
		vtx1.dir.normalize();
	}

	///< called by CSDA.
	///< no needs to get called by nuclear interactions
	CUDA_HOST_DEVICE
	void
	update_post_vertex_position
	(const R& len)
	{
		vtx1.pos  = vtx0.pos + vtx0.dir * len ;
	}

	///< 
	CUDA_HOST_DEVICE
	void
	update_post_vertex_energy
	(const R& e)
	{
		vtx1.ke  -= e;
	}
	
	///< Proceed a step
	///< vtx0 = vtx1
	CUDA_HOST_DEVICE
	void
	move()
	{
		vtx0 = vtx1 ;
		dE   = 0    ;
	}

	///< Change particle status
	CUDA_HOST_DEVICE
	void
	stop()
	{
		//this->move(); //don't put move here
		status = STOPPED;
	}

	///< Prints out track values CPU only as std libary is used
	/*
	CUDA_HOST
	void
	dump(){
		std::cout
			<< "Trak ID: " << id
			<< ", Particle: " << particle
			<< ", Status: " << status
			<< ", Vertex (Energy, pos, dir): "
			<< vtx0.ke <<", ("
			<< vtx0.pos.x <<", " << vtx0.pos.y <<", "<< vtx0.pos.z <<") ,"
			<< vtx0.dir.x <<", " << vtx0.dir.y <<", "<< vtx0.dir.z
			<<") to "
			<< vtx1.ke <<", ("
			<< vtx1.pos.x <<", " << vtx1.pos.y <<", "<< vtx1.pos.z <<") ,"
			<< vtx1.dir.x <<", " << vtx1.dir.y <<", "<< vtx1.dir.z <<") \n";			
	}
	*/
	
};

	
}
#endif

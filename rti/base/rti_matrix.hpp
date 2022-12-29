#ifndef RTI_MATRIX_H
#define RTI_MATRIX_H

/// \file
///
/// RT-Ion matrix 3x3 and 4x4 


#include <iostream>
#include <cmath>
#include <array>

#include <rti/base/rti_vec.hpp>

namespace rti{

/// \class mat3x3
/// 3x3 matrix
/// \tparam T type of matrix elements
template<typename T>
class mat3x3 {
public:
    ///< angle in radian
    T x;
    T y;
    T z;
    //matrix element
    T xx;
    T xy;
    T xz;
    T yx;
    T yy;
    T yz;
    T zx;
    T zy;
    T zz;
    
    CUDA_HOST_DEVICE
    mat3x3(): 
    x(0),y(0),z(0),
    xx(1.0),xy(0),xz(0),
    yx(0),yy(1.0),yz(0),
    zx(0),zy(0),zz(1.0)
    {;}

    CUDA_HOST_DEVICE
    mat3x3(T xx_, T xy_, T xz_,T yx_, T yy_, T yz_,T zx_, T zy_, T zz_ ):
    x(0),y(0),z(0),
    xx(xx_),xy(xy_),xz(xz_),
    yx(yx_),yy(yy_),yz(yz_),
    zx(zx_),zy(zy_),zz(zz_)
    {;}

    CUDA_HOST_DEVICE
    mat3x3(const mat3x3& ref) :
    x (ref.x ),
    y (ref.y ),
    z (ref.z ),
    xx(ref.xx),
    xy(ref.xy),
    xz(ref.xz),
    yx(ref.yx),
    yy(ref.yy),
    yz(ref.yz),
    zx(ref.zx),
    zy(ref.zy),
    zz(ref.zz)
    {;}

    CUDA_HOST_DEVICE
    mat3x3(T a, T b, T c):
        x(a),y(b),z(c),
        xx(1.0),xy(0),xz(0),
        yx(0),yy(1.0),yz(0),
        zx(0),zy(0),zz(1.0){
	//this routine may need a sequence of rotation,
	//e.g., x->y->z, y->x->z, etc. a total of 12
        if (x != 0) this->rotate_x(x);
        if (y != 0) this->rotate_y(y);
        if (z != 0) this->rotate_z(z);
    }

    CUDA_HOST_DEVICE
    mat3x3(std::array<T,3>& abc): 
        x(abc[0]),y(abc[1]),z(abc[2]),
        xx(1.0),xy(0),xz(0),
        yx(0),yy(1.0),yz(0),
        zx(0),zy(0),zz(1.0)
    {
        if (x != 0) this->rotate_x(x);
        if (y != 0) this->rotate_y(y);
        if (z != 0) this->rotate_z(z);
    }

    /// Assignment operator
    CUDA_HOST_DEVICE
    mat3x3<T>&
    operator=(const mat3x3<T>& ref)
    {
        x  = ref.x ;
        y  = ref.y;
        z  = ref.z;
        xx = ref.xx;
        xy = ref.xy;
        xz = ref.xz;
        yx = ref.yx;
        yy = ref.yy;
        yz = ref.yz;
        zx = ref.zx;
        zy = ref.zy;
        zz = ref.zz;
        return *this;
    }

    //Calculate rotation matrix from two vectors, from (f) and to (t).
    //the matrix rotates f and aligns t.
    // f and t are normalized one
	/*
    CUDA_HOST_DEVICE
    mat3x3
	(rti::vec3<T>& f, 
	 rti::vec3<T>& t)
    {
		f.normalize();
		t.normalize();
		printf("here\n");
		rti::vec3<T> v = f.cross(t);
		rti::vec3<T> u  = v / v.norm();
		T c = f.dot(t);
		T h = 1-c/1-c*c;

		xx = c + h * v.x * v.x   ;
		xy = h * v.x * v.y - v.z ;
		xz = h * v.x * v.z + v.y ;
		yx = h * v.x * v.y + v.z ;
		yy = c + h * v.y * v.y   ;
		yz = h * v.y * v.z - v.x ;
		zx = h * v.x * v.z - v.y ;
		zy = h * v.y * v.z + v.x ;
		zz = c + h * v.z * v.z   ;
		//calculat euler angle x->y->z sequence
    }
	*/
	
    //expecting normalized vectors, f and t
    CUDA_HOST_DEVICE
    mat3x3
	(const rti::vec3<T>& f, 
	 const rti::vec3<T>& t)
    {
		///< a matrix aligns vector (f) to vector (t)
		///< by Moller & Hughe, 1999
		rti::vec3<T> v = f.cross(t);
		T c = 1 ;
		T h = 0 ;
		if ( v.norm() <= 0.0001 ) {
			v.x = 0 ;
			v.y = 0 ;
			v.z = 0 ;
		}else{
			v.normalize();
			c = f.dot(t);
			h = (1.0-c)/(1.0-c*c);
		}
				
		xx = c + h * v.x * v.x   ;
		xy = h * v.x * v.y - v.z ;
		xz = h * v.x * v.z + v.y ;
		yx = h * v.x * v.y + v.z ;
		yy = c + h * v.y * v.y   ;
		yz = h * v.y * v.z - v.x ;
		zx = h * v.x * v.z - v.y ;
		zy = h * v.y * v.z + v.x ;
		zz = c + h * v.z * v.z   ;

    }

    CUDA_HOST_DEVICE
    ~mat3x3(){;}

    
    CUDA_HOST
    rti::vec3<T>
    euler_xyz(bool y_is_2nd_quad=false){
	//R_z(phi)*R_y(theta)*R_x(psi)
	//psi_th_phi.x = psi
	//psi_th_phi.y = th
	//psi_th_phi.z = phi
	rti::vec3<T> psi_th_phi;
	
	if( zx > -1.0 && zx < 1.0 ){ //zx neq -1 or 1
	    T psi[2], th[2], phi[2];
	    th[0] = -1.0 * std::asin(zx);
	    th[1] = M_PI - th[0];
	    psi[0] = std::atan2(zy/std::cos(th[0]), zz/std::cos(th[0]));
		psi[1] = std::atan2(zy/std::cos(th[1]), zz/std::cos(th[1]));
	    phi[0] = std::atan2(yx/std::cos(th[0]), xx/std::cos(th[0]));
		phi[1] = std::atan2(yx/std::cos(th[1]), xx/std::cos(th[1]));

	    if ( th[0] > 0 && th[0] < M_PI*0.5 ){
		psi_th_phi.x = psi[0];
	        psi_th_phi.y =  th[0];
	        psi_th_phi.z = phi[0];
	    }else{
		psi_th_phi.x = psi[1];
	        psi_th_phi.y =  th[1];
	        psi_th_phi.z = phi[1];
	    }
	    //
	    //psi_th_phi.x = psi[y_is_2nd_quad];
	    //psi_th_phi.y =  th[y_is_2nd_quad];
	    //psi_th_phi.z = phi[y_is_2nd_quad];
	}else{
	    std::cout<<"here\n";//never get here
	    psi_th_phi.z = 0.0 ; 
	    if ( zx == -1.0  ){//zx = -1
		psi_th_phi.x = psi_th_phi.z + std::atan2(xy, xz);
		psi_th_phi.y = M_PI * 0.5 ; 
	    }else{
		psi_th_phi.x = -1.0 * psi_th_phi.z + std::atan2(-1.0*xy, -1.0*xz);
		psi_th_phi.y = -1.0 * M_PI * 0.5 ;
	    }
	}
	
        #if defined(__CUDACC__)
		
	#else
	
	#endif
	return psi_th_phi;
    }
    
    CUDA_HOST_DEVICE
    mat3x3& 
    rotate_x(T a){
        x = a;
        #if defined(__CUDACC__)
        T c1 = cosf(x);
        T s1 = sinf(x);
        #else
        T c1 = std::cos(x);
        T s1 = std::sin(x);
        #endif
        T x1 = yx, y1 = yy, z1 = yz; 
        yx = c1*x1 - s1*zx;
        yy = c1*y1 - s1*zy;
        yz = c1*z1 - s1*zz;
        zx = s1*x1 + c1*zx;
        zy = s1*y1 + c1*zy;
        zz = s1*z1 + c1*zz;
        return *this;
    }

    CUDA_HOST_DEVICE
    mat3x3&
    rotate_y(T a){
        y = a;
        #if defined(__CUDACC__)
        T c1 = cosf(y);
        T s1 = sinf(y);
        #else
        T c1 = std::cos(y);
        T s1 = std::sin(y);
        #endif
        T x1 = zx, y1 = zy, z1 = zz; 
        zx = c1*x1 - s1*xx;
        zy = c1*y1 - s1*xy;
        zz = c1*z1 - s1*xz;
        xx = s1*x1 + c1*xx;
        xy = s1*y1 + c1*xy;
        xz = s1*z1 + c1*xz;
        return *this;
    }
    
    CUDA_HOST_DEVICE
    mat3x3&
    rotate_z(T a){
        z = a;
        #if defined(__CUDACC__)
        T c1 = cosf(z);
        T s1 = sinf(z);
        #else
        T c1 = std::cos(z);
        T s1 = std::sin(z);
        #endif
        T x1 = xx, y1 = xy, z1 = xz; 
        xx = c1*x1 - s1*yx;
        xy = c1*y1 - s1*yy;
        xz = c1*z1 - s1*yz;
        yx = s1*x1 + c1*yx;
        yy = s1*y1 + c1*yy;
        yz = s1*z1 + c1*yz;
        return *this;
    }

    CUDA_HOST_DEVICE
    std::array<T,3>
    operator * (const std::array<T,3>& r) const {
        return std::array<T,3>({ 
            xx * r[0] + xy * r[1] + xz * r[2],
            yx * r[0] + yy * r[1] + yz * r[2],
            zx * r[0] + zy * r[1] + zz * r[2]});
    }

    CUDA_HOST_DEVICE
    rti::vec3<T>
    operator * (const rti::vec3<T>& r) const {
        return rti::vec3<T>(
            xx * r.x + xy * r.y + xz * r.z,
            yx * r.x + yy * r.y + yz * r.z,
            zx * r.x + zy * r.y + zz * r.z);
    }

    CUDA_HOST_DEVICE
    mat3x3 inverse() const {        
        return rti::mat3x3<T>(xx, yx, zx, xy, yy, zy, xz, yz, zz );
    }

    /*
    mat3x3<T> 
    operator* (const mat3x3<T>& r){
        return mat3x3<T>(
            xx*r.xx + xy*r.yx + xz*r.zx,
            xx*r.xy + xy*r.yy + xz*r.zy,
            xx*r.xz + xy*r.yz + xz*r.zz,
            yx*r.xx + yy*r.yx + yz*r.zx,
            yx*r.xy + yy*r.yy + yz*r.zy,
            yx*r.xz + yy*r.yz + yz*r.zz,
            zx*r.xx + zy*r.yx + zz*r.zx,
            zx*r.xy + zy*r.yy + zz*r.zy,
            zx*r.xz + zy*r.yz + zz*r.zz );
    }*/

    CUDA_HOST_DEVICE
    mat3x3<T>
    operator* (const mat3x3<T>& r) const {
        return mat3x3<T>(
            xx*r.xx + xy*r.yx + xz*r.zx,
            xx*r.xy + xy*r.yy + xz*r.zy,
            xx*r.xz + xy*r.yz + xz*r.zz,
            yx*r.xx + yy*r.yx + yz*r.zx,
            yx*r.xy + yy*r.yy + yz*r.zy,
            yx*r.xz + yy*r.yz + yz*r.zz,
            zx*r.xx + zy*r.yx + zz*r.zx,
            zx*r.xy + zy*r.yy + zz*r.zy,
            zx*r.xz + zy*r.yz + zz*r.zz );
    }

    CUDA_HOST_DEVICE
    void dump() const {
    #if defined(__CUDACC__)
        printf("xx,xy,xz %f, %f, %f\n", xx,xy, xz);
        printf("yx,yy,yz %f, %f, %f\n", yx,yy, yz);
        printf("zx,zy,zz %f, %f, %f\n", zx,zy, zz);
    #else
        std::cout<<"xx,xy,xz "<< xx <<" " << xy <<" " << xz << std::endl;
        std::cout<<"yx,yy,yz "<< yx <<" " << yy <<" " << yz << std::endl;
        std::cout<<"zx,zy,zz "<< zx <<" " << zy <<" " << zz << std::endl;
    #endif
    }
};

/// \class mat4x4
/// 4x4 matrix
/// \tparam T type of matrix elements
template<typename T>
class mat4x4 {
public:
    //matrix elements
    T xx;
    T xy;
    T xz;
    T xs;
    T yx;
    T yy;
    T yz;
    T ys;
    T zx;
    T zy;
    T zz;
    T zs;
    T sx;
    T sy;
    T sz;
    T ss;
    
    CUDA_HOST_DEVICE
    mat4x4(): 
    xx(1.0),xy(0),xz(0),xs(0),
    yx(0),yy(1.0),yz(0),ys(0),
    zx(0),zy(0),zz(1.0),zs(0),
    sx(0),sy(0),sz(0.0),ss(1.0)
    {;}

    CUDA_HOST_DEVICE
    mat4x4(
        T xx_, T xy_, T xz_, T xs_,
        T yx_, T yy_, T yz_, T ys_,
        T zx_, T zy_, T zz_, T zs_,
        T sx_, T sy_, T sz_, T ss_
    ): 
    xx(xx_),xy(xy_),xz(xz_),xs(xs_),
    yx(yx_),yy(yy_),yz(yz_),ys(ys_),
    zx(zx_),zy(zy_),zz(zz_),zs(zs_),
    sx(sx_),sy(sy_),sz(sz_),ss(ss_)
    {;}

    CUDA_HOST_DEVICE
    mat4x4(const mat4x4& ref)
    {
        xx = ref.xx;
        xy = ref.xy;
        xz = ref.xz;
        xs = ref.xs;

        yx = ref.yx;
        yy = ref.yy;
        yz = ref.yz;
        ys = ref.ys;
        
        zx = ref.zx;
        zy = ref.zy;
        zz = ref.zz;
        zs = ref.zs;

        sx = ref.sx;
        sy = ref.sy;
        sz = ref.sz;
        ss = ref.ss;
    }

    CUDA_HOST_DEVICE
    mat4x4(const T* a)
    {
        xx = a[0]; xy = a[1]; xz = a[2]; xs = a[3];
        yx = a[4]; yy = a[5]; yz = a[6]; ys = a[7];
        zx = a[8]; zy = a[9]; zz = a[10]; zs = a[11];
        sx = a[12]; sy = a[13]; sz = a[14];ss = a[15];
    }

    CUDA_HOST_DEVICE
    mat4x4(
        const mat3x3<T>& rot, 
        const vec3<T>& tra)
    {
        xx = rot.xx;
        xy = rot.xy;
        xz = rot.xz;
        xs = tra.x;

        yx = rot.yx;
        yy = rot.yy;
        yz = rot.yz;
        ys = tra.y;
        
        zx = rot.zx;
        zy = rot.zy;
        zz = rot.zz;
        zs = tra.z;

        sx = 0.0;
        sy = 0.0;
        sz = 0.0;
        ss = 1.0;
    }

    CUDA_HOST_DEVICE
    mat4x4(const mat3x3<T>& rot)
    {
        xx = rot.xx;
        xy = rot.xy;
        xz = rot.xz;
        xs = 0.0;
        
        yx = rot.yx;
        yy = rot.yy;
        yz = rot.yz;
        ys = 0.0;
        
        zx = rot.zx;
        zy = rot.zy;
        zz = rot.zz;
        zs = 0.0;

        sx = 0.0;
        sy = 0.0;
        sz = 0.0;
        ss = 1.0;
    }

    CUDA_HOST_DEVICE
    mat4x4(const vec3<T>& ref)
    {
        xx = 1.0;
        xy = 0.0;
        xz = 0.0;
        xs = ref.x;

        yx = 0.0;
        yy = 1.0;
        yz = 0.0;
        ys = ref.y;
        
        zx = 0.0;
        zy = 0.0;
        zz = 1.0;
        zs = ref.z;

        sx = 0.0;
        sy = 0.0;
        sz = 0.0;
        ss = 1.0;
    }


    CUDA_HOST_DEVICE
    ~mat4x4(){;}

    CUDA_HOST_DEVICE
    std::array<T,4>
    operator * (const std::array<T,4>& r) const 
    {
        return std::array<T,4>({ 
            xx * r[0] + xy * r[1] + xz * r[2] + xs*r[3],
            yx * r[0] + yy * r[1] + yz * r[2] + ys*r[3],
            zx * r[0] + zy * r[1] + zz * r[2] + zs*r[3],
            sx * r[0] + sy * r[1] + sz * r[2] + ss*r[3]});
    }

    CUDA_HOST_DEVICE
    rti::vec4<T>
    operator * (const rti::vec4<T>& r) const 
    {
        return rti::vec4<T>(
            xx * r.x + xy * r.y + xz * r.z + xs * r.s,
            yx * r.x + yy * r.y + yz * r.z + ys * r.s,
            zx * r.x + zy * r.y + zz * r.z + zs * r.s,
            sx * r.x + sy * r.y + sz * r.z + ss * r.s);
    }

    CUDA_HOST_DEVICE
    rti::vec3<T>
    operator * (const rti::vec3<T>& r) const 
    {
        return rti::vec3<T>(
            xx * r.x + xy * r.y + xz * r.z + xs ,
            yx * r.x + yy * r.y + yz * r.z + ys ,
            zx * r.x + zy * r.y + zz * r.z + zs );
    }

    /*
    mat4x4<T> 
    operator* (const mat3x3<T>& r){
        return mat3x3<T>(
            xx*r.xx + xy*r.yx + xz*r.zx,
            xx*r.xy + xy*r.yy + xz*r.zy,
            xx*r.xz + xy*r.yz + xz*r.zz,
            yx*r.xx + yy*r.yx + yz*r.zx,
            yx*r.xy + yy*r.yy + yz*r.zy,
            yx*r.xz + yy*r.yz + yz*r.zz,
            zx*r.xx + zy*r.yx + zz*r.zx,
            zx*r.xy + zy*r.yy + zz*r.zy,
            zx*r.xz + zy*r.yz + zz*r.zz );
    }*/

    CUDA_HOST_DEVICE
    mat4x4<T>
    operator* (const mat4x4<T>& r) const 
    {
        return mat4x4<T>(
            xx*r.xx + xy*r.yx + xz*r.zx + xs*r.sx,
            xx*r.xy + xy*r.yy + xz*r.zy + xs*r.sy,
            xx*r.xz + xy*r.yz + xz*r.zz + xs*r.sz,
            xx*r.xs + xy*r.ys + xz*r.zs + xs*r.ss,
            
            yx*r.xx + yy*r.yx + yz*r.zx + ys*r.sx,
            yx*r.xy + yy*r.yy + yz*r.zy + ys*r.sy,
            yx*r.xz + yy*r.yz + yz*r.zz + ys*r.sz,
            yx*r.xs + yy*r.ys + yz*r.zs + ys*r.ss,
            
            zx*r.xx + zy*r.yx + zz*r.zx + zs*r.sx,
            zx*r.xy + zy*r.yy + zz*r.zy + zs*r.sy,
            zx*r.xz + zy*r.yz + zz*r.zz + zs*r.sz,
            zx*r.xs + zy*r.ys + zz*r.zs + zs*r.ss,
            
            sx*r.xx + sy*r.yx + sz*r.zx + ss*r.sx,
            sx*r.xy + sy*r.yy + sz*r.zy + ss*r.sy,
            sx*r.xz + sy*r.yz + sz*r.zz + ss*r.sz,
            sx*r.xs + sy*r.ys + sz*r.zs + ss*r.ss
            );
    }

    /*
    CUDA_HOST_DEVICE
    mat4x4 inverse() const {        
        T det = xx*(yy*zz*ss + yx*zs*sy + )
        //return rti::mat3x3<T>(xx, yx, zx, xy, yy, zy, xz, yz, zz );
    }
    */

    CUDA_HOST_DEVICE
    void dump() const 
    {
    #if defined(__CUDACC__)
        printf("xx,xy,xz,xs %f, %f, %f, %f\n", xx,xy, xz, xs);
        printf("yx,yy,yz,ys %f, %f, %f, %f\n", yx,yy, yz, ys);
        printf("zx,zy,zz,zs %f, %f, %f, %f\n", zx,zy, zz, zs);
        printf("sx,sy,sz,ss %f, %f, %f, %f\n", sx,sy, sz, ss);
    #else
        std::cout<<"xx,xy,xz,xs "<< xx <<" " << xy <<" " << xz << " " << xs << std::endl;
        std::cout<<"yx,yy,yz,ys "<< yx <<" " << yy <<" " << yz << " " << ys << std::endl;
        std::cout<<"zx,zy,zz,zs "<< zx <<" " << zy <<" " << zz << " " << zs << std::endl;
        std::cout<<"zx,zy,zz,ss "<< sx <<" " << sy <<" " << sz << " " << ss << std::endl;
    #endif
    }
};



}
#endif

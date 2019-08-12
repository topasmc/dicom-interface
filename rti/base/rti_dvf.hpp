#ifndef RTI_DVF_H
#define RTI_DVF_H

/// \file
///
/// Deformation Vector field

#include <vtkMetaImageReader.h>
#include <vtkSmartPointer.h>

//#include <vtkGDCMImageReader2.h>
#include <rti/base/rti_rect3d.hpp>

namespace rti{

/// \class dvf
/// this dvf class is a grid with vec3 pixel type.
/// \tparam S for type of vector elements, e.g. vec3<double> or vec3<float>
/// \tparam R for grid coordinate (float or double) a type of x, y, and z positions
/// \see http://dicom.nema.org/MEDICAL/dicom/2017b/output/chtml/part03/sect_C.20.3.html
template<typename S, typename R> 
class dvf : public rect3d<rti::vec3<S>, R> {

private:
    
    /// \note this is a temporal solution to have two reading capabilities
    bool is_dicom_; 

    /// DVF file name
    char* dvf_file_;

    /// DVF grid resolution.
    R dx_;
    R dy_;
    R dz_;

public: 

    /// Default constructor
    CUDA_HOST
    dvf(){;}

    
    /// Construct a rectlinear grid from DVF (deformation vector field) DICOM file
    /// \param f DVF file name.
    /// \param is_print flag to print DVF information
    CUDA_HOST
    dvf(std::string f, 
        bool is_print=false)
    {

        dvf_file_ = new char[f.length()+1];
        strcpy(dvf_file_, f.c_str());
        
        is_dicom_ = true;
        std::string fix = f.substr( f.size()- 3, 3 );
        //std::cout<< "file: " << f << ", extension: " << fix << " , compare: " << fix.compare("mha") << std::endl;
        if( !fix.compare("mhd") || 
            !fix.compare("mha") || 
            !fix.compare("MHD") || 
            !fix.compare("MHA") 
         ){
            is_dicom_ = false;
        }

        is_dicom_ ? this->read_dicom_reg() : this->read_mhd();

        
        std::cout<<"DVF: nx,ny,nz: (" 
                 << rect3d<rti::vec3<S>,R>::dim_.x << ", " 
                 << rect3d<rti::vec3<S>,R>::dim_.y <<", " 
                 << rect3d<rti::vec3<S>,R>::dim_.z <<")" << std::endl;
        std::cout<<"DVF (dx,dy,dz): (" << dx_ << ", " << dy_ <<", " << dz_ <<")" << std::endl;
        std::cout<<"DVF (x,y,z): (" << rect3d<rti::vec3<S>,R>::x_[0] << ", " 
                                    << rect3d<rti::vec3<S>,R>::y_[0] << ", " 
                                    << rect3d<rti::vec3<S>,R>::z_[0] <<")" << std::endl;
        this->flip_xyz_if_any();
    }

    /// Reads vector fields from MetaImage file
    CUDA_HOST
    virtual 
    void read_mhd()
    {
        assert("not yet implemented");
        
        vtkSmartPointer<vtkMetaImageReader> reader 
            = vtkSmartPointer<vtkMetaImageReader>::New();
        reader->SetFileName(dvf_file_);
        reader->Update();
        
        int ext[6];
        reader->GetDataExtent(ext);

        rect3d<rti::vec3<S>,R>::dim_.x = ext[1] - ext[0] + 1; 
        rect3d<rti::vec3<S>,R>::dim_.y = ext[3] - ext[2] + 1;
        rect3d<rti::vec3<S>,R>::dim_.z = ext[5] - ext[4] + 1;

        double* space  = reader->GetDataSpacing();

        dx_ = space[0]; dy_ = space[1]; dz_ = space[2];

        double* origin = reader->GetDataOrigin();
        rect3d<rti::vec3<S>,R>::x_ = new R[rect3d<rti::vec3<S>,R>::dim_.x]; 
        rect3d<rti::vec3<S>,R>::y_ = new R[rect3d<rti::vec3<S>,R>::dim_.y];
        rect3d<rti::vec3<S>,R>::z_ = new R[rect3d<rti::vec3<S>,R>::dim_.z];

        for(size_t i=0; i < rect3d<rti::vec3<S>,R>::dim_.x ; ++i){
            rect3d<rti::vec3<S>,R>::x_[i] = i*dx_ + origin[0];
        }
        
        for(size_t i=0; i < rect3d<rti::vec3<S>,R>::dim_.y ; ++i){
            rect3d<rti::vec3<S>,R>::y_[i] = i*dy_ + origin[1]; 
        } 

        for(size_t i=0; i < rect3d<rti::vec3<S>,R>::dim_.z ; ++i){
            rect3d<rti::vec3<S>,R>::z_[i] = i*dz_ + origin[2];
        }
        
    }

    /// Reads vector fields from DICOM-REG file
    CUDA_HOST
    virtual 
    void read_dicom_reg(){

        gdcm::Reader reader; reader.SetFileName(dvf_file_);
        if(!reader.Read()) throw std::runtime_error("Invalid DICOM file is given to interface.");

        ///< Check the media storage
        const gdcm::DataSet& header = reader.GetFile().GetHeader();
        gdcm::Attribute<0x0002,0x0002> tag4ms;
        tag4ms.SetFromDataElement(header.GetDataElement(tag4ms.GetTag()));
        std::string ms = tag4ms.GetValue();
        if (ms.compare("1.2.840.10008.5.1.4.1.1.66.3")){
            throw std::runtime_error("It's not MediaStorage:Deformable Spation Registration Storage");
        }

        ///< Read position, resolution, voxels
        const gdcm::DataSet& ds = reader.GetFile().GetDataSet();

        ///< Access Deformable Registration Sequence
        gdcm::Tag t_drs(0x0064,0x0002);
        if(!ds.FindDataElement(t_drs)){
            std::runtime_error("Can't find Deformable Registration Sequence from "+ std::string(dvf_file_));
        }

        const gdcm::DataElement& drs = ds.GetDataElement(t_drs);
        gdcm::SmartPointer<gdcm::SequenceOfItems> p_drs = drs.GetValueAsSQ();

        if( p_drs->GetNumberOfItems() == 0){
            std::runtime_error("Can't find Deformable Registration Sequence item");
        }

        ///< Deformable Registration 2
        /// \note MGH DVF has DVF information in sequence item 2
        /// Not sure this is general. 
        gdcm::Item& item2 = p_drs->GetItem(2); //
        gdcm::DataSet& ds_r2 = item2.GetNestedDataSet();
        
        ///< Deformable Dose Grid: dimension, data 2
        ///tag for Deformable registration grid sequence
        gdcm::Tag t_drgs(0x0064,0x0005);     
        const gdcm::DataElement& drgs = ds_r2.GetDataElement(t_drgs);
        gdcm::SmartPointer<gdcm::SequenceOfItems> p_drgs = drgs.GetValueAsSQ();
        
        if( p_drgs->GetNumberOfItems() == 0){
            std::runtime_error("Can't find Deformable Registration Grid Sequence");
        }

        gdcm::Item& p_drgs_i1 = p_drgs->GetItem(1);
        const gdcm::DataSet& ds_drgs = p_drgs_i1.GetNestedDataSet();
        
        gdcm::Attribute<0x0020,0x0032> a_pos;
        gdcm::Attribute<0x0020,0x0037> a_ori;
        gdcm::Attribute<0x0064,0x0007> a_nxyz;
        gdcm::Attribute<0x0064,0x0008> a_dxyz;

        a_pos.SetFromDataElement(ds_drgs.GetDataElement(a_pos.GetTag())); //[XYZ]_start
        a_ori.SetFromDataElement(ds_drgs.GetDataElement(a_ori.GetTag())); //[XYZ]_row/column 
        a_nxyz.SetFromDataElement(ds_drgs.GetDataElement(a_nxyz.GetTag())); //[]
        a_dxyz.SetFromDataElement(ds_drgs.GetDataElement(a_dxyz.GetTag()));

        ///
        /// Get Vector information described in DICOM. C.20.3.1.1
        /// [X Y Z]_column
        /// [X Y Z]_depth
        /// [X Y Z]_R
        //auto xyz_mm  = a_pos.GetValues();
        //* [X Y Z]_start
        rti::vec3<R> xyz_mm( a_pos.GetValues()[0], a_pos.GetValues()[1], a_pos.GetValues()[2]);
        
        rti::vec4<R> xyz0( a_pos.GetValues()[0], a_pos.GetValues()[1], a_pos.GetValues()[2], 1.0);

        rti::vec3<R> row(a_ori.GetValues()[0], a_ori.GetValues()[1], a_ori.GetValues()[2]);
        rti::vec3<R> col(a_ori.GetValues()[3], a_ori.GetValues()[4], a_ori.GetValues()[5]);

        rti::vec3<uint> nxyz(a_nxyz.GetValues());

        rti::vec3<R> dxyz(a_dxyz.GetValues()[0], a_dxyz.GetValues()[1], a_dxyz.GetValues()[2]);

        rti::vec3<R> depth = row.cross(col);

        rti::mat4x4<R> init_m(
            row.x * dxyz.x, col.x * dxyz.y, depth.x * dxyz.z, xyz_mm.x, 
            row.y * dxyz.x, col.y * dxyz.y, depth.y * dxyz.z, xyz_mm.y, 
            row.z * dxyz.x, col.z * dxyz.y, depth.z * dxyz.z, xyz_mm.z,
            0.0, 0.0, 0.0, 1.0
        );
        
        ///< Pre Deformable registration matrix sequence 
        /// tag for Pre Deformable registration Matrix sequence
        gdcm::Tag t_pre_dmrs(0x0064,0x000f); 
        rti::mat4x4<R> pre_m = this->transformation_matrix(ds_r2,t_pre_dmrs);

        /// tag for Post Deformable registration Matrix sequence        
        gdcm::Tag t_post_dmrs(0x0064,0x0010); 
        rti::mat4x4<R> post_m = this->transformation_matrix(ds_r2, t_post_dmrs);

        rti::mat4x4<R> matrix_final = post_m*pre_m*init_m;
        
        rect3d<rti::vec3<S>,R>::dim_.x = nxyz.x; 
        rect3d<rti::vec3<S>,R>::dim_.y = nxyz.y;
        rect3d<rti::vec3<S>,R>::dim_.z = nxyz.z;
        dx_ = dxyz.x; dy_ = dxyz.y; dz_ = dxyz.z;

        rect3d<rti::vec3<S>,R>::x_ = new R[rect3d<rti::vec3<S>,R>::dim_.x]; 
        rect3d<rti::vec3<S>,R>::y_ = new R[rect3d<rti::vec3<S>,R>::dim_.y];
        rect3d<rti::vec3<S>,R>::z_ = new R[rect3d<rti::vec3<S>,R>::dim_.z];

        /// Initialize x_,y_,z_ vector
        /// we assumed that dose grid has same voxel size, 
        /// especially we assumed that Z framebuffer is same
        /// Direction cosine is assumed to be (1/-1,0,0,0,1/-1,0,0,0,1/-1)

        for(size_t i=0; i < rect3d<rti::vec3<S>,R>::dim_.x ; ++i){
            rect3d<rti::vec3<S>,R>::x_[i] = i*matrix_final.xx + matrix_final.xs;
        }
        
        for(size_t i=0; i < rect3d<rti::vec3<S>,R>::dim_.y ; ++i){
            rect3d<rti::vec3<S>,R>::y_[i] = i*matrix_final.yy + matrix_final.ys; 
        } 

        for(size_t i=0; i < rect3d<rti::vec3<S>,R>::dim_.z ; ++i){
            rect3d<rti::vec3<S>,R>::z_[i] = i*matrix_final.zz + matrix_final.zs;
        }

    }

    ///Loads DVFs
    CUDA_HOST
    virtual 
    void
    load_data()
    {
        size_t nb_voxels 
            = rect3d<rti::vec3<S>,R>::dim_.x 
            * rect3d<rti::vec3<S>,R>::dim_.y 
            * rect3d<rti::vec3<S>,R>::dim_.z;
                          
        rect3d<rti::vec3<S>,R>::data_.resize(nb_voxels);

        is_dicom_ ? this->load_dicom() : this->load_meta();

        this->flip_data();

    }

    ///Loads DVFs from DICOM
    CUDA_HOST
    virtual
    void 
    load_dicom()
    {
        gdcm::Reader reader; reader.SetFileName(dvf_file_);
        if(!reader.Read()) throw std::runtime_error("Invalid DICOM file is given to interface.");

        /**< Read position, resolution, voxels */
        const gdcm::DataSet& ds = reader.GetFile().GetDataSet();

        //tag for Deformable registration sequence
        gdcm::Tag t_drs(0x0064,0x0002); 
        const gdcm::DataElement& drs = ds.GetDataElement(t_drs);
        gdcm::SmartPointer<gdcm::SequenceOfItems> p_drs = drs.GetValueAsSQ();
        
        //!!!! temporal, MGH DVF has DVF information in sequence item 2
        //Not sure this is general 
        gdcm::Item& item2 = p_drs->GetItem(2); //!!!!!
        /**< Deformable Registration 2 */
        gdcm::DataSet& ds_r2 = item2.GetNestedDataSet();
        
        /**< Deformable Dose Grid: dimension, data 2 */
        //tag for Deformable registration grid sequence
        gdcm::Tag t_drgs(0x0064,0x0005);     
        const gdcm::DataElement& drgs = ds_r2.GetDataElement(t_drgs);
        gdcm::SmartPointer<gdcm::SequenceOfItems> p_drgs = drgs.GetValueAsSQ();
        
        gdcm::Item& p_drgs_i1 = p_drgs->GetItem(1);
        const gdcm::DataSet& ds_drgs = p_drgs_i1.GetNestedDataSet();
        
        gdcm::Tag pixeldata(0x0064,0x0009);
        gdcm::Attribute<0x0064,0x0009> dvf_data;

        const gdcm::ByteValue* bv = ds_drgs.GetDataElement(pixeldata).GetByteValue();

        size_t nb_voxels 
            = rect3d<rti::vec3<S>,R>::dim_.x 
            * rect3d<rti::vec3<S>,R>::dim_.y 
            * rect3d<rti::vec3<S>,R>::dim_.z;

        std::valarray<float> tmp_data(nb_voxels*3);
        bv->GetBuffer( (char*) &tmp_data[0], nb_voxels*3*sizeof(float) );
        
        //http://dicom.nema.org/MEDICAL/dicom/2017b/output/chtml/part17/chapter_P.html
        //8.9944\-6.24271\0.00387573

        this->transform_to_type(tmp_data);     
    }

    ///Loads DVFs from MetaImage
    CUDA_HOST 
    virtual 
    void
    load_meta()
    {
        //not yet implemented
        vtkSmartPointer<vtkMetaImageReader> reader 
            = vtkSmartPointer<vtkMetaImageReader>::New();
        reader->SetFileName(dvf_file_);
        reader->Update();

        reader->GetDataScalarType();

        void* p = reader->GetMemoryBuffer(); 
        std::cout<< std::sizeof(*p) << std::endl;
        //GetBuffer((char*) &rect3d<int16_t,R>::data_[i*nb_voxels_2d]);
        //reader->GetBu


    }

    /// Finds x index
    /// \param x  x-position
    /// \return x-index
    /// Index-finding is simplified cause pixel spacing is same along x
    CUDA_HOST_DEVICE
    inline virtual 
    size_t
    find_c000_x_index( const R& x)
    {   
        assert( (x >= rect3d<rti::vec3<S>,R>::x_[0]) && (x < rect3d<rti::vec3<S>,R>::x_[rect3d<vec3<S>,R>::dim_.x-1]) );
        return floor(abs(x-rect3d<rti::vec3<S>,R>::x_[0])/dx_);
    }

    /// Finds y index
    /// \param y  y-position
    /// \return y-index
    CUDA_HOST_DEVICE
    inline virtual size_t
    find_c000_y_index(const R& y)
    {
        assert( (y >= rect3d<rti::vec3<S>,R>::y_[0]) && (y < rect3d<rti::vec3<S>,R>::y_[rect3d<rti::vec3<S>,R>::dim_.y-1]) );
        
        return floor( abs(y- rect3d<rti::vec3<S>,R>::y_[0])/dy_);
    }

    /// Finds z index
    /// \param z  z-position
    /// \return z-index
    CUDA_HOST_DEVICE
    inline virtual size_t
    find_c000_z_index(const R& z)
    {
        assert( (z >= rect3d<rti::vec3<S>,R>::z_[0]) && (z < rect3d<rti::vec3<S>,R>::z_[rect3d<rti::vec3<S>,R>::dim_.z-1]) );
        return floor(abs(z-rect3d<rti::vec3<S>,R>::z_[0])/dz_);
    }

    /// Write data into file
    CUDA_HOST
    virtual void
    write_data(const std::string filename)
    {

	    std::ofstream file1( filename, std::ios::out | std::ofstream::binary);
        size_t nb_voxels =  rect3d<rti::vec3<S>,R>::dim_.x 
                          * rect3d<rti::vec3<S>,R>::dim_.y 
                          * rect3d<rti::vec3<S>,R>::dim_.z;
                          
        std::valarray<S> tmp0(nb_voxels*3); //temporal copy object
        for(size_t i=0; i < nb_voxels ; ++i){
            tmp0[3*i]   = rect3d<rti::vec3<S>,R>::data_[i].x;
            tmp0[3*i+1] = rect3d<rti::vec3<S>,R>::data_[i].y;
            tmp0[3*i+2] = rect3d<rti::vec3<S>,R>::data_[i].z;
        }        

        /*
        size_t pixel_id = 0 ;
        std::cout << rect3d<rti::vec3<S>,R>::data_[pixel_id].x <<", "
                  << rect3d<rti::vec3<S>,R>::data_[pixel_id].y <<", "
                  << rect3d<rti::vec3<S>,R>::data_[pixel_id].z <<"  "
                  <<std::endl;

        std::cout << tmp0[3*pixel_id+0] <<", "
                  << tmp0[3*pixel_id+1] <<", "
                  << tmp0[3*pixel_id+2] <<"  "
                  <<std::endl;
        std::cout << "Max: " << tmp0.max() <<", Min:" << tmp0.min() <<"  "
                  <<std::endl;
        */

        file1.write(reinterpret_cast<const char *>(&tmp0[0]), tmp0.size() * sizeof(S) );
        file1.close();

    }

private:
    
    /// Convert integer type PixelData to float 
    /// by using std::transform
    CUDA_HOST
    void
    transform_to_type(std::valarray<float>& in)
    {

        float nan_voxels = std::count_if(begin(in), end(in),[](float i){return std::isnan(i);});
        std::cout<<"number of NaN values:" << nan_voxels << std::endl;

        size_t nb_voxels = rect3d<rti::vec3<S>,R>::dim_.x 
                          * rect3d<rti::vec3<S>,R>::dim_.y 
                          * rect3d<rti::vec3<S>,R>::dim_.z;

        for(size_t i=0; i < nb_voxels ; ++i){
            rect3d<rti::vec3<S>,R>::data_[i] = rti::vec3<S>( in[3*i],in[3*i+1],in[3*i+2] );
        }
    }


    /// Get deformation Pre/Post matrix from Deformable Registration Grid Sequence
    /// As these tags are not available in gdcm dictionary, gdcm's approach was applied.
    CUDA_HOST
    rti::mat4x4<R>
    transformation_matrix(
        const gdcm::DataSet& ds,
        gdcm::Tag& pre_or_post)
    {

        const gdcm::DataElement& dspdrms = ds.GetDataElement(pre_or_post);
        gdcm::SmartPointer<gdcm::SequenceOfItems> p_dspdrms = dspdrms.GetValueAsSQ();
        
        if( p_dspdrms->GetNumberOfItems() == 0){
            std::runtime_error("Can't find Deformable Registration matrix sequence");
        }

        gdcm::Item& p_dspdrm_i1 = p_dspdrms->GetItem(1);
        const gdcm::DataSet& ds_pdmr = p_dspdrm_i1.GetNestedDataSet();
        
        gdcm::Attribute<0x0070,0x030c> a_mtype;  ///< matrix type ~ RIGID
        gdcm::Attribute<0x3006,0x00c6> a_matrix; ///< Pre-matrix 16 elements for 4x4 matrix
        a_mtype.SetFromDataElement(ds_pdmr.GetDataElement(a_mtype.GetTag()));
        a_matrix.SetFromDataElement(ds_pdmr.GetDataElement(a_matrix.GetTag()));
        const double* m = a_matrix.GetValues();
        rti::mat4x4<R> txm( m[0],  m[1], m[2], m[3],
                            m[4],  m[5], m[6], m[7],
                            m[8],  m[9], m[10],m[11],
                            m[12], m[13],m[14],m[15]);
        return txm;
    }
};

}

#endif
#ifndef RTI_RTDOSE_H
#define RTI_RTDOSE_H

/// \file
///
/// DICOM-RT dose

#include <rti/base/rti_rect3d.hpp>

namespace rti{

/// \class rtdose 
/// rtdose has a float as pixel type.
/// \tparam R for grid coordinate (float or double) a type of x, y, and z positions
template<typename R>
class rtdose : public rect3d<float,R> {
protected:

    /// path to rtdose file
    char* rtdose_file;

    /// Dimension of rtdose cube
    rti::vec3<R> lxyz;

public: 
    /// Default constructor
    CUDA_HOST
    rtdose(){;}

    /// Constructs a rectlinear grid from array of x/y/z with their size
    /// \param f rtdose file name.
    CUDA_HOST
    rtdose(std::string f)
    {
        rtdose_file = new char[f.length()+1];
        strcpy(rtdose_file, f.c_str());

        gdcm::ImageReader reader;
        reader.SetFileName(rtdose_file);
        reader.Read();
        const gdcm::Image& img  = reader.GetImage();

        /// We assumed that dose grid has same voxel size, 
        /// especially we assumed that Z framebuffer is same
        /// Direction cosine is assumed to be (1,0,0,0,1,0)
        auto xyz_mm    = img.GetOrigin()  ;
        auto dxyz_mm   = img.GetSpacing() ;
        auto nxyz      = img.GetDimensions();
        //std::cout<< img << std::endl;

        lxyz.x = dxyz_mm[0] * nxyz[0];
        lxyz.y = dxyz_mm[1] * nxyz[1];
        lxyz.z = dxyz_mm[2] * nxyz[2];
        
        rect3d<float,R>::dim_.x = nxyz[0]; 
        rect3d<float,R>::dim_.y = nxyz[1];
        rect3d<float,R>::dim_.z = nxyz[2];

        /// Initializes x_,y_,z_ vector
        rect3d<float,R>::x_ = new R[rect3d<float,R>::dim_.x]; 
        rect3d<float,R>::y_ = new R[rect3d<float,R>::dim_.y];
        rect3d<float,R>::z_ = new R[rect3d<float,R>::dim_.z];

        for(size_t i=0; i < rect3d<float,R>::dim_.x ; ++i) rect3d<float,R>::x_[i] = xyz_mm[0] + i*dxyz_mm[0];
        for(size_t i=0; i < rect3d<float,R>::dim_.y ; ++i) rect3d<float,R>::y_[i] = xyz_mm[1] + i*dxyz_mm[1];

        const gdcm::DataSet& ds = reader.GetFile().GetDataSet();
        
        //gdcm::Attribute<0x0028,0x0009> frame_increment_ptr ;
        //frame_increment_ptr.SetFromDataElement(ds.GetDataElement(frame_increment_ptr.GetTag()));
        //assert(frame_increment_ptr.GetValue() == gdcm::Tag(0x3004,0x000c) && "GridFrameOffsetVector is missing.");

        gdcm::Attribute<0x3004,0x000c> gridframe_offsetvector ;
        gridframe_offsetvector.SetFromDataElement(ds.GetDataElement(gridframe_offsetvector.GetTag()));
        auto offset_vector = gridframe_offsetvector.GetValues();
        auto offset_size   = gridframe_offsetvector.GetNumberOfValues();
        
        for(size_t i=0; i < offset_size ; ++i) rect3d<float,R>::z_[i] = xyz_mm[2] + offset_vector[i];

        std::cout<<"RTDOSE (nx,ny,nz): (" << nxyz[0] << ", " << nxyz[1] <<", " << nxyz[2] <<")" << std::endl;
        std::cout<<"RTDOSE (dx,dy,dz): (" << dxyz_mm[0] << ", " << dxyz_mm[1] <<", " << dxyz_mm[2] <<")" << std::endl;
        std::cout<<"RTDOSE (x,y,z): (" << rect3d<float,R>::x_[0] << ", " << rect3d<float,R>::y_[0] <<", " << rect3d<float,R>::z_[0] <<")" << std::endl;
    }


    
    /// Initializes data by reading PixelData in RTDOSE
    /// This operation needs to be done in CPU side
    CUDA_HOST
    virtual void
    load_data()
    {
        size_t nb_voxels = rect3d<float,R>::dim_.x * rect3d<float,R>::dim_.y * rect3d<float,R>::dim_.z;
        rect3d<float,R>::data_.resize(nb_voxels);

        gdcm::ImageReader reader;
        reader.SetFileName(rtdose_file);
        reader.Read();
        const gdcm::Image& img  = reader.GetImage();
        float intercept = float(img.GetIntercept());
        float slope     = float(img.GetSlope());

        gdcm::PixelFormat pixeltype = img.GetPixelFormat();
        
        switch( pixeltype )
        {
            case gdcm::PixelFormat::INT8:
                {
                std::valarray<int8_t> d(nb_voxels) ;
                img.GetBuffer((char*) &d[0]);
                this->transform_to_float(d, intercept, slope);
                }
                break;
            case gdcm::PixelFormat::UINT8:
                {
                std::valarray<uint8_t> d(nb_voxels) ;
                img.GetBuffer((char*) &d[0]);
                this->transform_to_float(d, intercept, slope);
                }
                break;
            case gdcm::PixelFormat::INT16:
                {
                std::valarray<int16_t> d(nb_voxels) ;
                img.GetBuffer((char*) &d[0]);
                this->transform_to_float(d, intercept, slope);
                }
                break;
            case gdcm::PixelFormat::UINT16:
                {
                std::valarray<uint16_t> d(nb_voxels) ;
                img.GetBuffer((char*) &d[0]);
                this->transform_to_float(d, intercept, slope);
                }
                break;
            case gdcm::PixelFormat::INT32:
                {
                std::valarray<int32_t> d(nb_voxels) ;
                img.GetBuffer((char*) &d[0]);
                this->transform_to_float(d, intercept, slope);
                }
                break;
            case gdcm::PixelFormat::UINT32:
                {
                //Only this Pixelformat was tested.
                //couldn't put buffer to rect3d img.GetBuffer( (char*) &rect3d<float,R>::data_[0] );
                std::valarray<uint32_t> d(nb_voxels) ;
                img.GetBuffer((char*) &d[0]); 
                this->transform_to_float(d, intercept, slope);
                }
                break;
            default:
                assert(0);
        }
        
    }


    
    /// Converts integer type PixelData to float 
    /// by using std::transform
    CUDA_HOST
    template<class S>
    void
    transform_to_float(
        std::valarray<S>& in, 
        float inter, 
        float slope)
    {
        std::transform(begin(in), end(in),             //from
                       begin(rect3d<float,R>::data_),  //to
                       [inter, slope](S a) { return (float(a))*slope + inter; }
                       );
    }

    /// Returns the size of rtdose box
    CUDA_HOST
    rti::vec3<R>
    get_dosegrid_size(){
        return lxyz;
    }

};



}

#endif
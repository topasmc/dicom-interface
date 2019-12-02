cmake_minimum_required (VERSION 3.10)

project (dicom-interface)

add_library(dicom_interface INTERFACE)

target_compile_features(dicom_interface INTERFACE cxx_std_11)
target_include_directories(dicom_interface INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR}
)

# This fixes strange whitespace errors when trying to compile against GDCM
set(GDCM_USE_VTK OFF)

# We need at least GDCM version 2.6
find_package(GDCM 2.6 REQUIRED)

target_link_libraries(dicom_interface INTERFACE gdcmMSFF)

# This make GDCM header files visible to rti
target_include_directories(dicom_interface INTERFACE ${GDCM_INCLUDE_DIRS})

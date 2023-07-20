cmake_minimum_required (VERSION 2.8.11)
project (session)

#find_package (GDCM REQUIRED)
#include (${GDCM_USE_FILE})

set(CMAKE_CXX_STANDARD 11) 
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON) 

#file(GLOB gdcmlibs "${GDCM_LIBRARY_DIRS}/*.so*" "${GDCM_LIBRARY_DIRS}/*.a" "${GDCM_LIBRARY_DIRS}/*.dylib")
#file(GLOB excludegdcmlibs "${GDCM_LIBRARY_DIRS}/libgdcmMEXD.*" "${GDCM_LIBRARY_DIRS}/libsocketxx.*")

find_library (PATHS /opt/sw/topas/3.9.0/Frameworks/ )

if(APPLE)
	find_library(COREFOUNDATION_LIBRARY CoreFoundation )
endif()

include_directories (
    ${PROJECT_SOURCE_DIR}
    /opt/sw/dicom-interface/git.public.20221230/gdcm-2.6
    ../../../
)

link_directories (
    ${PROJECT_SOURCE_DIR}
    /opt/sw/topas/3.9.0/Frameworks/ 
)

add_executable(treatment_session treatment_session.cpp)
add_executable(smc_linescan smc_linescan.cpp)

#install(PROGRAMS treatment_session DESTINATION .)

set (GDCM_Libs 
        gdcmMSFF 
        gdcmDICT 
        gdcmIOD 
        gdcmexpat 
        gdcmDSED 
        gdcmCommon 
        gdcmzlib 
        gdcmjpeg8 
        gdcmjpeg12 
        gdcmjpeg16 
        gdcmopenjpeg 
        gdcmcharls 
        gdcmuuid)  


target_link_libraries(treatment_session
    ${COREFOUNDATION_LIBRARY}
    ${GDCM_Libs}
)

target_link_libraries(smc_linescan
    ${COREFOUNDATION_LIBRARY}
    ${GDCM_Libs}
)



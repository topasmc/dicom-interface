cmake_minimum_required (VERSION 2.8.11)
project (session)

find_package (GDCM REQUIRED)
include (${GDCM_USE_FILE})

set(CMAKE_CXX_STANDARD 11) 
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON) 

file(GLOB gdcmlibs "${GDCM_LIBRARY_DIRS}/*.so*" "${GDCM_LIBRARY_DIRS}/*.a" "${GDCM_LIBRARY_DIRS}/*.dylib")
file(GLOB excludegdcmlibs "${GDCM_LIBRARY_DIRS}/libgdcmMEXD.*" "${GDCM_LIBRARY_DIRS}/libsocketxx.*")

if(APPLE)
	find_library(COREFOUNDATION_LIBRARY CoreFoundation )
endif()


include_directories (
    ${PROJECT_SOURCE_DIR},
    ../../../
)

add_executable(treatment_session treatment_session.cpp)

install(PROGRAMS treatment_session DESTINATION .)

link_directories (${GDCM_DIR} )

target_link_libraries(treatment_session
    ${COREFOUNDATION_LIBRARY}
    gdcmCommon
    gdcmDSED     
    gdcmMEXD     
    gdcmcharls   
    gdcmjpeg12   
    gdcmjpeg8    
    gdcmuuid     
    socketxx
    gdcmDICT     
    gdcmIOD      
    gdcmMSFF     
    gdcmexpat    
    gdcmjpeg16   
    gdcmzlib   
)



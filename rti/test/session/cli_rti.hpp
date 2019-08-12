#ifndef CLI_RTI_HPP
#define CLI_RTI_HPP


#include <array>
#include <iostream>
#include <valarray>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <fstream>
#include <limits>

#include <sys/stat.h>
//#include <errno.h>
#include <fcntl.h>
#include <unistd.h>
#include <cstdlib>
#include <sys/mman.h>


namespace rti{

/**
 *  A command line interface for Unit tests
 *  cli is a top virtual class
 *  1. reading RT-Ion plan file and its type (plan or log).
 *  2. creating a geometries, patient, dosegrid, and beamline components of a machine
 *  3. creating beam sources of a machine
 */
class cli{

protected:
    std::map<const std::string, std::vector<std::string> > parameters; 

public:
    cli(){;};
    ~cli(){;}

	void 
	read(int argc, char** argv){

		std::cout<<"# of arguments: "<< argc<<std::endl;
	    if (argc ==1) {print_help(argv[0]); exit(1);};

        for(int i=1; i < argc ; ++i){
			//Find a parameter
		    auto it = parameters.find(argv[i]);
		    if( it != parameters.end()){
		    	//Accumulate argv until it meets next option		
		    	int j = i+1;
		    	do{
		    		if(std::string(argv[j]).compare(0,2,"--") == 0) break;
		    		it->second.push_back(argv[j])	;
		    		j++;
		    	}while( j < argc );
		    	//Print out options 
		    	std::cout<<it->first<<" : " ;
		    	for(auto parm:it->second) std::cout<< parm <<" ";
		    	std::cout<<std::endl;
		    }
	    }
	}

    const
    std::vector<std::string>
    operator[](
        const std::string& t)
    {
        return parameters[t];
    }

    virtual void print_help(char* s){;}
	
};



/**
 *  A command line interface for Unit tests of patient geometry
 *
 */
class cli_patient : public cli{

public:
    cli_patient(){
		parameters = 
		std::map<const std::string, std::vector<std::string> > 
		({ 
			//option name,  {parameters}
			{"--ctdir"        ,{} }
		});
	}

    ~cli_patient(){;}

	virtual void 
	print_help(
		char* s
	)
	{
        std::cout<<"Usage:   "<<s<<" [-option] [argument]"<<std::endl;
        std::cout<<"options:  "<<"--ctdir CT_DIRECTORY "<<std::endl;
    }
};

/**
 *  A command line interface for Unit tests of dose grid
 *  it reads in RTDOSE and returns the dimension and pixel-pitch
 *  or 
 *  it reads in CT-images and returns the dosegrid to cover patient image
 */
class cli_dosegrid : public cli{

public:
    cli_dosegrid(){
		parameters = 
		std::map<const std::string, std::vector<std::string> > 
		({ 
			//option name,  {parameters}
			{"--ctdir"    ,{} },
			{"--dosefile" ,{} },
			{"--pxyz"     ,{} }
		});
	}

    ~cli_dosegrid(){;}

	virtual void 
	print_help(
		char* s
	)
	{
        std::cout<<"Usage:   "<<s<<" [-option] [argument]"<<std::endl;
        std::cout<<"options:  "<<"--ctdir CT_DIRECTORY "<<std::endl;
		std::cout<<"          "<<"--dosefile rtdosefile "<<std::endl;
		std::cout<<"          "<<"--pxyz x,y,z "<<std::endl;
    }
};

/**
 *  A command line interface of unit tests for DVF (deformation vector field)
 *  it reads in RTDOSE, CT, DVF
 *  and performs
 *  1. interpolation mathing RTDOSE to RTDOSE in CT grid 
 *  2. perform dose warping with linear transportation
 */
class cli_dvf : public cli{

public:
    cli_dvf(){
		parameters = 
		std::map<const std::string, std::vector<std::string> > 
		({ 
			//option name,  {parameters}
			{"--ctdir"    ,{} },
			{"--rtdose" ,{} },
			{"--dvf"   ,{} },
			{"--output" ,{}}
		});
	}

    ~cli_dvf(){;}

	virtual void 
	print_help(
		char* s
	)
	{
        std::cout<<"Usage:   "<<s<<" [-option] [argument]"<<std::endl;
        std::cout<<"options:  "<<"--ctdir CT_DIRECTORY "<<std::endl;
		std::cout<<"          "<<"--rtdose rtdosefile "<<std::endl;
		std::cout<<"          "<<"--dvf    REG-DICOM "<<std::endl;
		std::cout<<"          "<<"--output  output file (binary) "<<std::endl;
    }
};

/**
 *  A command line interface for Unit tests of basic dicom reading capabilities
 *
 */
class cli_rti_read : public cli{

public:
    cli_rti_read(){
		parameters = 
			std::map<const std::string, std::vector<std::string>> 
			({ 
				//option name,  {parameters}
				{"--rti"        ,{} },  //RTIP or RTIBTR
				{"--bname"      ,{} },  //beam name
				{"--bnumber"    ,{} },  //beam number
				{"--machine"    ,{} }  //specify machine name to use
			});
	}
    ~cli_rti_read(){;}

	virtual void 
	print_help(
		char* s
	)
	{
        std::cout<<"Usage:   "<<s<<" [-option] [argument]"<<std::endl;
        std::cout<<"options:  "<<"--rti RTIP or RTIBTR "<<std::endl;  //mandatory
		std::cout<<"          "<<"--machine name"<<std::endl;         //default '' 
        std::cout<<"          "<<"--bname number"<<std::endl;         //default 0
        std::cout<<"          "<<"--bnumber beamname"<<std::endl;     //no default
    }
    
};

/**
 *  A command line interface for Unit tests of beam generation
 *
 */
class cli_beam_read : public cli{

public:
    cli_beam_read(){
		parameters = std::map<const std::string, std::vector<std::string> >({ 
		//option name,  {parameters}
		{"--rti"        ,{} },  //RTIP or RTIBTR
		{"--machine"    ,{} },  //treatment machine: Institution:Tx_machine_name
		{"--mc_code"    ,{} },  //mc_code, e.g., topas:3.0.1 (default), gpmc
		{"--bname"      ,{} },  //beam name
		{"--bnumber"    ,{} },  //beam number
		{"--spots"       ,{} }, //spot-id, obsolute, will be replaced by beamlets
		{"--beamlets"       ,{} }, //beamlet-id
		{"--pph"        ,{} },  //particles_per_history
		{"--pxyz"       ,{} },  //position xyz w.r.t mother coordinate system
		{"--rxyz"       ,{} },  //rotation xyz w.r.t mother coordinate system 
		{"--sid"        ,{} },  //beam generation distance from isocenter.
		{"--output"     ,{} }   //output file name, for empty, such as "", it printout to console
		});
	}

    ~cli_beam_read(){;}

	virtual void 
	print_help(
		char* s
	)
	{
        std::cout<<"Usage:   "<<s<<" [-option] [argument]"<<std::endl;
        std::cout<<"options:  "<<"--rti RTIP or RTIBTR "<<std::endl;
        std::cout<<"          "<<"--machine machine_name"<<std::endl; //default, automatic search
		std::cout<<"          "<<"--mc_code  mccode"<<std::endl; //default topas
        std::cout<<"          "<<"--bname number"<<std::endl;         //default 0
        std::cout<<"          "<<"--bnumber beamname"<<std::endl;     //no default
		std::cout<<"          "<<"--spots i for single, i j for range, -1 for all."<<std::endl;  //e.g) i:i-th spot, i-j:i-th to j-th, -1: all
    	std::cout<<"          "<<"--beamlets i for single, i j for range, no options for all."<<std::endl;  //e.g) i:i-th spot, i-j:i-th to j-th, -1: all
		std::cout<<"          "<<"--pph scale #-1 is 1 history/spot"<<std::endl;            //default -1 -> 1 histories 
        std::cout<<"          "<<"--pxyz x y z"<<std::endl;           //position of phase-space coordinates system
        std::cout<<"          "<<"--rxyz x y z"<<std::endl;           //rotation of phase-space coordinates system
		std::cout<<"          "<<"--sid d"<<std::endl;           //distance between beam and isocenter in z direction
        std::cout<<"          "<<"--output "<<std::endl;              //output file name
    }
};

}

#endif
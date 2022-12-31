## DICOM RT-Ion interface for MC simulation 
***

<!-- 
title: "Note: DICOM interface for MC simulation"
author: Jungwook Shin
date: July 1, 2019
-->

> *<span style="color:red">Current version is BETA as of Aug14, 2019</span>*

DICOM RT-Ion interface, shortly **RTI**, is a library to convert treatment planning information in DICOM format into Monte Carlo components on-the-fly such as geometries and beam source. The **RTI** allows integration of DICOM-RT Ion interface within a MC engine to provide more reliable and consistent performance than with a script-based processing to interpret DICOM information.
For the Monte Carlo simulation, [TOPAS extensions](https://topas.readthedocs.io/en/latest/extension-docs/intro.html#) of geometry and particle source to connect **RTI** are provided.

DICOM-RT objects utilized in this interface are
- Ion Plan (**RTIP**)
- Ion Beams Treatment Record (**RTIBTR**) 
- CT image 
- RT-Dose
> note: RT Plan is not supported (2019. July 22).

Beamline geometries, gantry and patient coordinate systems, and fluence map are determined from RTIP and/or RTIBTR.  CT image and Dose are used only for patient dose calculations.

## Getting started

### Requirements

- **[GDCM](http://gdcm.sourceforge.net)** (tested versions are 2.4 and 2.6.) library (header and objects) is required.
  tested versions are 2.4 and 2.6.8
  version 3.0 doesn't work with RTI (thank to Dohyeon).

### Monte Carlo 

- **[TOPAS](http://www.topasmc.org)** are tested for versions 3.1, 3.2, 3.5, 3.7, and 3.9

### Tested operating systems
  - Mac OS X
    - 10.14.5, clang++ (LLVM version 10.0.1, /usr/bin/c++) 
    - 12.6.2,  clang++ (LLVM version 14.0.0, /usr/bin/c++) with topas 3.9
    - cmake version 3.24
  - Linux
    - Red Hat Enterprise Linux Server release 6.7 g++ (4.9.0)
    - Ubuntu 19.10
  
### Installation

As **RTI** is a header only library, copy or download the files to somewhere in your system.
```bash
$ cd /<your_sw_path>/
$ git clone https://github.com/topasmc/dicom-interface rti.git
```

If you are compiling dicom-interface with topas releases ( not from source), you need to use same headers of gdcm that topas used.
```bash
$ cd /<your_sw_path>/rti.git
$ tar -zxf gdcm-2.6.include.tar.gz
$ ls gdcm-2.6
```

Then, go to the directory where your topas is installed. 

```bash
$ cd /<your_topas_path>/
```

Add two lines in CMakeList.txt (in TOPAS) as following:
```cmake
include_directories (
    ${CMAKE_BINARY_DIR}
    ${CMAKE_BINARY_DIR}/extensions
    ${PROJECT_SOURCE_DIR}/Geant4Headers
    ${PROJECT_SOURCE_DIR}/include
    /<your_sw_path>/rti.git/
    /<your_sw_path>/rti.git/gdcm-2.6
)
```


```bash
$ cmake -DTOPAS_EXTENSIONS_DIR=/<your_sw_path>/rti.git/rti/topas/rtion .
$ make -j4
```

### Run a test

```bash
$ cd /your_sw_path/rti.git/topas/tutorial/
$ topas beam_view.txt
```

For more detail, please see this project's wiki page [https://github.com/topasmc/dicom-interface/wiki](https://github.com/topasmc/dicom-interface/wiki).

## Authors

- **Shin, Jungwook**: main developments and documentations. 
- **Hueso-Gonzalez, Fernando**: documentations and discussion
- **Edmunds, David**: GPU (cuda) implementation
- **Kooy, Hanne, M.** : Advisor
- **Paganetti, Harald** : Advisor
- **Clasie, Benjamin M.** : Advisor
  
## Acknowledgements

This project is initiated by Massachusetts General Hospital and continuously supported by a NIH grant, TOPAS project (U24).

## How to cite

Please cite this paper: https://doi.org/10.1016/j.ejmp.2020.04.018

## License (TBD):

This project is licensed under - see the LICENSE.md file for details.

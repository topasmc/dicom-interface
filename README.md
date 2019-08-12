## DICOM RT-Ion interface for MC simulation
***

<!-- 
title: "Note: DICOM interface for MC simulation"
author: Jungwook Shin
date: July 1, 2019
-->

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

- **[TOPAS](http://www.topasmc.org)** (version 3.2) is required and assumed pre-installed.
- **[GDCM](http://gdcm.sourceforge.net)** (version > 2.6.8) library (header and objects) is required.

### Tested operating systems
  - Mac OS X
    - 10.14.5, clang++ (LLVM version 10.0.1, /usr/bin/c++)
    - cmake version 3.x
  - Red Hat Enterprise Linux Server release 6.7
    - g++ (4.9.0)

### Installation

As **RTI** is a header only library, copy or download the files to somewhere in your system.
```bash
$ cd /<your_sw_path>/
$ git clone https://github.com/mghro/rtion-interface-mc.git rti.git
```

Then, go to the directory where your topas is installed. 

```bash
$ cd /<your_topas_path>/
```

Add two lines in CMakeList.txt as following:
```cmake
include_directories (
    ${CMAKE_BINARY_DIR}
    ${CMAKE_BINARY_DIR}/extensions
    ${PROJECT_SOURCE_DIR}/Geant4Headers
    ${PROJECT_SOURCE_DIR}/include
    /<your_sw_path>/rti.git/
    /<your_gdcm_include_path>/gdcm-<version>/
)
```


```bash
$ cmake -DTOPAS_EXTENSIONS_DIR=/<your_sw_path>/rti.git/topas/rtion .
$ make -j4
```

### Run a test

```bash
$ cd /your_sw_path/rti.git/test/topas/
$ topas rti_beam.txt
```

For more detail, please see this project's wiki page [https://github.com/mghro/rtion-interface-mc/wiki](https://github.com/mghro/rtion-interface-mc/wiki).

## Authors

- **Shin, Jungwook**: main developments and documentations. 
- **Hueso-Gonzalez, Fernando**: documentations and discussion
- **Edmunds, David**: GPU (cuda) implementation
- **Kooy, Hanne, M.** : Advisor
- **Paganetti, Harald** : Advisor
- **Clasie, Benjamin M.** : Advior
  
## Acknowledgements

This project is initiated by Massachusetts General Hospital and continuously supported by a NIH grant, TOPAS project (U24).

## License (TBD):

This project is licensed under - see the LICENSE.md file for details.

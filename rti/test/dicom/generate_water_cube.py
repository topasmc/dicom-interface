#!/usr/bin/env python3
#__author__ = 'Fernando Hueso-Gonzalez'
# python3 generate_water_cube.py --outdir /tmp/watercube50/

import argparse
import sys
import os

from pydicom import uid, write_file, dcmread
from pydicom.dataset import Dataset, FileDataset
from pydicom.sequence import Sequence

from scipy.interpolate import *
from math import sqrt, log
import numpy as np

def main(dir, args):

    parser = argparse.ArgumentParser(description='Generate an ideal CT of a water cube and a dummy RTplan with certain energy layers, defined in script. Every energy layer has a single spot at (0,0)')

    parser.add_argument('--outdir', dest='outdir'   , type=str, required=True )
    parser.add_argument('--inDate', dest='inDate'   , type=str, required=False, default='20191011')
    parser.add_argument('--stDate', dest='stDate'   , type=str, required=False, default='20191011')
    parser.add_argument('--seDate', dest='seDate'   , type=str, required=False, default='20191011')
    parser.add_argument('--acDate', dest='acDate'   , type=str, required=False, default='20191011')
    parser.add_argument('--coDate', dest='coDate'   , type=str, required=False, default='20191011')
    parser.add_argument('--inTime', dest='inTime'   , type=str, required=False, default='150000')
    parser.add_argument('--stTime', dest='stTime'   , type=str, required=False, default='150000')
    parser.add_argument('--seTime', dest='seTime'   , type=str, required=False, default='150000')
    parser.add_argument('--acTime', dest='acTime'   , type=str, required=False, default='150000')
    parser.add_argument('--coTime', dest='coTime'   , type=str, required=False, default='150000')
    parser.add_argument('--access', dest='access'   , type=str, required=False, default='0')
    parser.add_argument('--manuf' , dest='manuf'    , type=str, required=False, default='MGH Physics Research')
    parser.add_argument('--station',dest='station'  , type=str, required=False, default='Nashua')
    parser.add_argument('--institution', dest='institution', type=str, required=False, default='rbe')
    parser.add_argument('--instaddr', dest='instaddr', type=str, required=False, default='30 Fruit St, Boston MA 02114, USA')
    parser.add_argument('--physician', dest='physician', type=str, required=False, default='')
    parser.add_argument('--stDesc', dest='stDesc'  , type=str, required=False, default='For commissioning')
    parser.add_argument('--seDesc', dest='seDesc'  , type=str, required=False, default='Water cube 50x50x50cm3')
    parser.add_argument('--operat', dest='operat'  , type=str, required=False, default='FHG')
    parser.add_argument('--manmod', dest='manmod'  , type=str, required=False, default='2019')
    parser.add_argument('--patname',dest='patname' , type=str, required=False, default='Cube^Water')
    parser.add_argument('--patid' , dest='patid'   , type=str, required=False, default='WaterCube50')
    parser.add_argument('--birth' , dest='birth'   , type=str, required=False, default='')
    parser.add_argument('--age'   , dest='age'     , type=str, required=False, default='')
    parser.add_argument('--sex'   , dest='sex'     , type=str, required=False, default='')
    parser.add_argument('--anon'  , dest='anon'    , type=str, required=False, default='NO')
    parser.add_argument('--sw'    , dest='sw'      , type=str, required=False, default='pydicom1.2.2')
    parser.add_argument('--lastCal',dest='lastCal' , type=str, required=False, default='')
    parser.add_argument('--stuid' , dest='stuid'   , type=str, required=False, default=uid.generate_uid())
    parser.add_argument('--seuid' , dest='seuid'   , type=str, required=False, default=uid.generate_uid())
    parser.add_argument('--rseuid' , dest='rseuid' , type=str, required=False, default=uid.generate_uid())
    parser.add_argument('--stid'  , dest='stid'    , type=str, required=False, default='1')
    parser.add_argument('--sen'   , dest='sen'     , type=str, required=False, default='1')
    parser.add_argument('--acqn'  , dest='acqn'    , type=str, required=False, default='1')
    parser.add_argument('--forUID', dest='forUID'  , type=str, required=False, default=uid.generate_uid())
    parser.add_argument('--RTlabel',dest='RTlabel' , type=str, required=False, default='Test')
    parser.add_argument('--RTname', dest='RTname'  , type=str, required=False, default='Pristine test')
    parser.add_argument('--RTdesc', dest='RTdesc'  , type=str, required=False, default='Dummy energies')
    parser.add_argument('--RTdate', dest='RTdate'  , type=str, required=False, default='20191011')
    parser.add_argument('--RTtime', dest='RTtime'  , type=str, required=False, default='150000')
    parser.add_argument('--RTgeom', dest='RTgeom'  , type=str, required=False, default='PATIENT')
    parser.add_argument('--sadX'  , dest='sadX'    , type=float, required=False, default=2000.0)
    parser.add_argument('--sadY'  , dest='sadY'    , type=float, required=False, default=1800.0)
    parser.add_argument('--nBlocks',dest='nBlocks' , type=int, required=False, default=0)
    parser.add_argument('--isoX'  , dest='isoX'    , type=float, required=False, default=0.0)
    parser.add_argument('--isoY'  , dest='isoY'    , type=float, required=False, default=0.0)
    parser.add_argument('--isoZ'  , dest='isoZ'    , type=float, required=False, default=0.0)
    parser.add_argument('--snout' , dest='snout'   , type=str, required=False, default='')
    parser.add_argument('--nRS'   , dest='nRS'     , type=int, required=False, default=0)
    parser.add_argument('--rsID'  , dest='rsID'    , type=str, required=False, default='')
    parser.add_argument('--snPos' , dest='snPos'   , type=float, required=False, default=200)
    parser.add_argument('--seX'   , dest='seX'     , type=float, required=False, default=0.0)
    parser.add_argument('--seY'   , dest='seY'     , type=float, required=False, default=0.0)
    parser.add_argument('--seZ'   , dest='seZ'     , type=float, required=False, default=0.0)
    parser.add_argument('--ssd'   , dest='ssd'     , type=float, required=False, default=2000.0)
    parser.add_argument('--beam'        , dest='beam'       , type=str  , required=False, nargs='?', default='G000'  )
    parser.add_argument('--machine'     , dest='machine'    , type=str  , required=False, nargs='?', default='1.1')
    parser.add_argument('--dosimeter'   , dest='dosimeter'  , type=str  , required=False, nargs='?', default='NP'    )
    parser.add_argument('--tuneid'      , dest='tuneid'     , type=str  , required=False, nargs='?', default='Tune'  )
    parser.add_argument('--gangle'      , dest='gangle'     , type=int  , required=False, nargs='?', default=0       )
    # Define geometry of water cube
    parser.add_argument('--rows'    , dest='rows'    , type=int  , required=False, nargs='?', default=512)
    parser.add_argument('--columns' , dest='columns' , type=int  , required=False, nargs='?', default=512)
    parser.add_argument('--slices'  , dest='slices'  , type=int  , required=False, nargs='?', default=512)
    parser.add_argument('--wrows'   , dest='wrows'   , type=float, required=False, nargs='?', default=1  )# mm
    parser.add_argument('--wcolumns', dest='wcolumns', type=float, required=False, nargs='?', default=1  )# mm
    parser.add_argument('--wslices' , dest='wslices' , type=float, required=False, nargs='?', default=1  )# mm
    parser.add_argument('--margin'  , dest='margin'  , type=int  , required=False, nargs='?', default=6  )# pixels of air around water cube, so that water cube is 50cm * 50cm * 50cm
    # Proton energies and their sigma
    parser.add_argument('--energies', dest='energies', type=float, required=False, nargs='?', default=[90 + i*5 for i in range(0,29)])# MeV
    parser.add_argument('--lat_sigma',dest='lat_sigma',type=float, required=False, nargs='?', default=[6 for i in range(0,29)])# mm
    parser.add_argument('--pPerSpot' ,dest='pPerSpot', type=float, required=False, nargs='?', default=[1e9 for i in range(0,29)])# protons per spot
   
    args = parser.parse_args(args)
    rows = args.rows
    columns = args.columns
    slices = args.slices
    wrows = args.wrows
    wcolumns = args.wcolumns
    wslices = args.wslices
    margin = args.margin
    energies = args.energies
    lat_sigma = args.lat_sigma
    pPerSpot = args.pPerSpot
    
    # Check output dir
    if not os.path.exists(args.outdir):
        print('Creating output folder',args.outdir)
        os.mkdir(args.outdir)
        os.mkdir(os.path.join(args.outdir,'ct'))
    else:
        print('Writing to preexisting output folder',args.outdir)
        if not os.path.exists(os.path.join(args.outdir,'ct')):
            print('Creating output subfolder',os.path.join(args.outdir,'ct'))
            os.mkdir(os.path.join(args.outdir,'ct'))

    cube = np.full((slices, rows, columns), -1000, dtype='int16') # -1000 HU as initialization value (air)
    cube[margin:-margin,margin:-margin,margin:-margin] = 0 # 0 HU cube (water)

    # Generate CT image
    for sl in range(slices):
        # ~ break
        output_file = os.path.join(args.outdir,'ct',str(sl+1)+'.dcm')
        image = cube[sl]
        print(output_file,image.shape)

        meta = Dataset()
        meta.MediaStorageSOPClassUID = '1.2.840.10008.5.1.4.1.1.2'#'CT Image Storage'
        meta.MediaStorageSOPInstanceUID = uid.generate_uid()
        meta.ImplementationClassUID = uid.PYDICOM_IMPLEMENTATION_UID
        meta.TransferSyntaxUID = uid.ImplicitVRLittleEndian
        ds = FileDataset(output_file, {}, file_meta=meta, preamble=b"\0" * 128)

        ds.SpecificCharacterSet     = 'ISO_IR 100'
        ds.ImageType                = ['DERIVED','SECONDARY','AXIAL']
        ds.InstanceCreationDate     = args.inDate
        ds.InstanceCreationTime     = args.inTime
        ds.SOPClassUID              = ds.file_meta.MediaStorageSOPClassUID
        ds.SOPInstanceUID           = ds.file_meta.MediaStorageSOPInstanceUID
        ds.StudyDate                = args.stDate
        ds.SeriesDate               = args.seDate
        ds.AcquisitionDate          = args.acDate
        ds.ContentDate              = args.coDate
        ds.StudyTime                = args.stTime
        ds.SeriesTime               = args.seTime
        ds.AcquisitionTime          = args.acTime
        ds.ContentTime              = args.coTime
        ds.AccessionNumber          = args.access
        ds.Modality                 = 'CT'
        ds.Manufacturer             = args.manuf
        ds.InstitutionName          = args.institution
        ds.InstitutionAddress       = args.instaddr
        ds.ReferringPhysicianName   = args.physician
        ds.StationName              = args.station
        ds.StudyDescription         = args.stDesc
        ds.SeriesDescription        = args.seDesc
        ds.OperatorsName            = args.operat
        ds.ManufacturerModelName    = args.manmod
        ds.PatientName              = args.patname
        ds.PatientID                = args.patid
        ds.PatientBirthDate         = args.birth
        ds.PatientSex               = args.sex
        ds.PatientAge               = args.age
        ds.PatientIdentityRemoved   = args.anon
        ds.AdditionalPatientHistory = 'Pseudo-CT'
        if args.anon == "YES":
            ds.DeidentificationMethod   = "Manual"
        ds.SliceThickness           = wslices
        # ~ ds.FocalSpots
        # ~ ds.KVP
        # ~ ds.DataCollectionDiameter
        # ~ ds.ReconstructionDiameter
        ds.SoftwareVersions         = args.sw
        # ~ ds.DistanceSourceToDetector
        # ~ ds.DistanceSourceToPatient
        # ~ ds.GantryDetectorTilt
        # ~ ds.ExposureTime
        # ~ ds.XRayTubeCurrent
        # ~ ds.RotationDirection
        # ~ ds.ConvolutionKerne
        # ~ ds.FilterType
        ds.ProtocolName             = 'RESEARCH'
        # ~ ds.ScanOptions = 'AXIAL MODE'
        ds.DateOfLastCalibration    = args.lastCal
        ds.PatientPosition          = 'HFS'
        ds.StudyInstanceUID         = args.stuid
        ds.SeriesInstanceUID        = args.seuid
        ds.StudyID                  = args.stid
        ds.SeriesNumber             = args.sen
        ds.AcquisitionNumber        = args.acqn
        ds.InstanceNumber           = sl + 1
        # Zero (origin) is at surface of water on anterior side, and centered in the other axes. Image position refers to center point of corner voxel.
        ds.ImagePositionPatient     = [(-columns/2+0.5)*wcolumns,(-margin+0.5)*wrows, (slices/2-0.5-sl)*wslices]
        ds.ImageOrientationPatient  = ['1', '0', '0', '0', '1', '0']
        ds.FrameOfReferenceUID      = args.forUID
        ds.PositionReferenceIndicator = ''#'OM'
        ds.SliceLocation            = ds.ImagePositionPatient[2]

        ds.SamplesPerPixel          = 1
        ds.PhotometricInterpretation= 'MONOCHROME2'
        ds.Rows                     = rows
        ds.Columns                  = columns
        ds.PixelSpacing             = [wcolumns, wrows]
        ds.BitsAllocated            = 16
        ds.BitsStored               = 16
        ds.HighBit                  = 15
        ds.PixelRepresentation      = 1
        ds.SmallestImagePixelValue  = np.amin(image).tobytes()
        ds.LargestImagePixelValue   = np.amax(image).tobytes()
        ds.WindowCenter             = 0
        ds.WindowWidth              = 1000
        ds.RescaleIntercept         = 0
        ds.RescaleSlope             = 1
        ds.RescaleType              = "HU"
        ds.PixelData                = image.tobytes()
        # ~ ds.PixelPaddingValue        = -2000

        ds.save_as(output_file, False)

    # Generate RTplan with certain energies
    meta = Dataset()
    meta.MediaStorageSOPClassUID    = '1.2.840.10008.5.1.4.1.1.481.8' #RT Ion Plan Storage
    meta.MediaStorageSOPInstanceUID = uid.generate_uid()
    meta.ImplementationClassUID     = uid.PYDICOM_IMPLEMENTATION_UID
    meta.TransferSyntaxUID = uid.ImplicitVRLittleEndian
    output_file = os.path.join(args.outdir,'rtplan.dcm')
    ds = FileDataset(output_file, {}, file_meta=meta, preamble=b"\0" * 128)

    ds.SpecificCharacterSet     = 'ISO_IR 100'
    ds.InstanceCreationDate     = args.inDate
    ds.InstanceCreationTime     = args.inTime
    ds.SOPClassUID              = ds.file_meta.MediaStorageSOPClassUID
    ds.SOPInstanceUID           = ds.file_meta.MediaStorageSOPInstanceUID
    ds.StudyDate                = args.stDate
    ds.SeriesDate               = args.seDate
    ds.StudyTime                = args.stTime
    ds.SeriesTime               = args.seTime
    ds.AccessionNumber          = ''
    ds.Modality                 = 'RTPLAN'
    ds.Manufacturer             = args.manuf
    ds.InstitutionName          = args.institution
    ds.ReferringPhysicianName   = args.physician
    ds.StudyDescription         = args.stDesc
    ds.SeriesDescription        = args.seDesc
    ds.OperatorsName            = args.operat
    ds.ManufacturerModelName    = args.manmod
    ds.PatientName              = args.patname
    ds.PatientID                = args.patid
    ds.PatientBirthDate         = args.birth
    ds.PatientSex               = args.sex
    ds.PatientAge               = args.age
    ds.PatientIdentityRemoved   = args.anon
    if args.anon=="YES":
        ds.DeidentificationMethod   = "Manual"
    ds.SoftwareVersions         = args.sw
    ds.DateOfLastCalibration    = args.lastCal
    ds.StudyInstanceUID         = args.stuid
    ds.SeriesInstanceUID        = args.rseuid
    ds.StudyID                  = args.stid
    ds.SeriesNumber             = args.sen
    ds.InstanceNumber           = args.acqn
    ds.FrameOfReferenceUID      = args.forUID
    ds.PositionReferenceIndicator = ''
    ds.RTPlanLabel              = args.RTlabel
    ds.RTPlanName               = args.RTname
    ds.RTPlanDescription        = args.RTdesc
    ds.RTPlanDate               = args.RTdate
    ds.RTPlanTime               = args.RTtime
    ds.RTPlanGeometry           = args.RTgeom

    ds.FractionGroupSequence = Sequence()
    dsfx = Dataset()
    dsfx.FractionGroupNumber      = 1
    dsfx.FractionGroupDescription = ''
    dsfx.NumberOfFractionsPlanned = 1
    dsfx.NumberOfBeams            = 1
    dsfx.NumberOfBrachyApplicationSetups = 0
    dsfx.ReferencedBeamSequence = Sequence()
    dsfx_b = Dataset()
    dsfx_b.BeamDoseSpecificationPoint = [args.isoX,args.isoY,args.isoZ]
    dsfx_b.BeamDose = 1 #dummy
    dsfx_b.BeamMeterset = float(sum(pPerSpot))#pPerSpot protons per spot
    dsfx_b.ReferencedBeamNumber = 1
    dsfx.ReferencedBeamSequence.append(dsfx_b)
    ds.FractionGroupSequence.append(dsfx)

    ds.PatientSetupSequence = Sequence()
    pss = Dataset()
    pss.PatientPosition = 'HFS'
    pss.PatientSetupNumber = 1
    pss.PatientSetupLabel = 'Standard'
    ds.PatientSetupSequence.append(pss)

    ds.IonBeamSequence = Sequence()
    be = Dataset()
    ds.IonBeamSequence.append(be)
    be.BeamName = args.beam
    be.IonControlPointSequence = Sequence()
    be.TreatmentMachineName = args.machine
    be.InstitutionName = args.institution
    be.PrimaryDosimeterUnit = args.dosimeter
    be.Manufacturer = args.manuf
    be.InstitutionName  = args.institution
    be.ManufacturerModelName  = args.machine
    be.InstitutionAddress = args.instaddr
    be.TreatmentMachineName   = args.machine
    be.PrimaryDosimeterUnit   = args.dosimeter
    be.BeamNumber             = '1'
    be.BeamName               = args.beam
    be.BeamDescription        = 'Gantry from top'
    be.BeamType               = 'STATIC'
    be.RadiationType          = 'PROTON'
    be.TreatmentDeliveryType  = 'TREATMENT'
    be.NumberOfWedges         = 0
    be.NumberOfCompensators   = 0
    be.NumberOfBoli           = 0
    be.NumberOfBlocks         = args.nBlocks
    be.FinalCumulativeMetersetWeight = int(sum(pPerSpot)) # pPerSpot protons per spot (energy)
    be.NumberOfControlPoints  = 2*len(energies)
    be.ScanMode                   = 'MODULATED'
    be.VirtualSourceAxisDistances = [args.sadX, args.sadY]
    if not args.snout == '':
        be.SnoutSequence = Sequence()
        sds = Dataset()
        sds.AccessoryCode = args.snout
        sds.SnoutID = args.snout
        be.SnoutSequence.append(sds)
    be.NumberOfRangeShifters   = args.nRS
    if args.nRS == 1:
        be.RangeShifterSequence = Sequence()
        rsds = Dataset()
        rsds.AccessoryCode = 'Undefined Accessory Code'
        rsds.RangeShifterNumber = 1
        rsds.RangeShifterID = args.rsID
        rsds.RangeShifterType = 'BINARY'
        be.RangeShifterSequence.append(rsds)
    be.NumberOfLateralSpreadingDevices = 0
    be.NumberOfRangeModulators = 0
    be.PatientSupportType = 'TABLE'
    cweight = 0
    for i,energy in enumerate(energies[::-1]):
        for j in range(2):
            icpoi = Dataset()
            icpoi.NominalBeamEnergyUnit = 'MEV'
            icpoi.ControlPointIndex = i*2 + j
            icpoi.NominalBeamEnergy = str(energy)
            if j == 0:
                icpoi.GantryAngle = args.gangle
                icpoi.GantryRotationDirection = 'NONE'
                icpoi.BeamLimitingDeviceAngle = 0
                icpoi.BeamLimitingDeviceRotationDirection = 'NONE'
                icpoi.PatientSupportAngle = 0
                icpoi.PatientSupportRotationDirection = 'NONE'
                icpoi.TableTopVerticalPosition     = 0
                icpoi.TableTopLongitudinalPosition = 0
                icpoi.TableTopLateralPosition      = 0
                icpoi.IsocenterPosition            = [args.isoX,args.isoY,args.isoZ]
                icpoi.SurfaceEntryPoint            = [args.seX,args.seY,args.seZ]
                icpoi.SourceToSurfaceDistance      = args.ssd
            icpoi.CumulativeMetersetWeight = cweight
            if j == 0:
                icpoi.TableTopPitchAngle = 0
                icpoi.TableTopPitchRotationDirection = 'NONE'
                icpoi.TableTopRollAngle  = 0
                icpoi.TableTopRollRotationDirection = 'NONE'
                icpoi.GantryPitchAngle = 0.0
                icpoi.GantryPitchRotationDirection = 'NONE'
                icpoi.SnoutPosition      = args.snPos
            icpoi.ScanSpotTuneID = args.tuneid
            icpoi.NumberOfScanSpotPositions = 1
            icpoi.ScanSpotPositionMap = [0.0, 0.0]
            icpoi.ScanSpotMetersetWeights = pPerSpot[len(energies)-i-1] if j == 0 else 0
            cweight += icpoi.ScanSpotMetersetWeights
            icpoi.ScanningSpotSize = [lat_sigma[len(energies)-i-1],lat_sigma[len(energies)-i-1]]
            icpoi.NumberOfPaintings = 1
            be.IonControlPointSequence.append(icpoi)
    be.ReferencedPatientSetupNumber = 1
    be.ReferencedToleranceTableNumber = 0

    ds.ApprovalStatus = 'UNAPPROVED'
    # ~ ds.ReviewDate     = args.inDate
    # ~ ds.ReviewTime     = args.inTime
    # ~ ds.ReviewerName   = 'You'

    ds.save_as(output_file, False)
    print('Done, saved to',output_file)

if __name__ == "__main__":
    main(sys.argv[0], sys.argv[1:])

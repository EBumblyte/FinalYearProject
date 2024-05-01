"""
Created on Sat Apr 13 14:43:19 2024

@author: Emilija Bumblyte(20431382)

Description:
    
    
"""

from astropy import units as u
import numpy as np
import ccdproc
from astropy.io import fits
from astropy.nddata import CCDData


bias1_0 = fits.getdata('E:/T120/Bias/CCD Image 2379.fits')
bias2 = fits.getdata('E:/T120/Bias/CCD Image 2380.fits')
bias3 = fits.getdata('E:/T120/Bias/CCD Image 2381.fits')
bias4 = fits.getdata('E:/T120/Bias/CCD Image 2382.fits')
bias5 = fits.getdata('E:/T120/Bias/CCD Image 2383.fits')
bias1 = CCDData(bias1_0, unit='adu')
bias2 = CCDData(bias2, unit='adu')
bias3 = CCDData(bias3, unit='adu')
bias4 = CCDData(bias4, unit='adu')
bias5 = CCDData(bias5, unit='adu')
MasterBias = ccdproc.combine([bias1, bias2, bias3, bias4, bias5], 'E:/T120/Bias/MasterBias3.fit', 'median')



flat1 = fits.getdata('E:/T120/Flat/Hafilter/CCD Image 2410.fits')
flat2 = fits.getdata('E:/T120/Flat/Hafilter/CCD Image 2411.fits')
flat1 = CCDData(flat1, unit ='adu')
flat2 = CCDData(flat2, unit ='adu')
MasterFlat = ccdproc.combine([flat1, flat2],'E:/T120/Flat/Hafilter/MasterFlatcopy3.fit', 'median')

"""
dark1 = fits.getdata('E:/T120/Dark/CCD Image 2441.fits')
dark1 = CCDData(dark1, unit = 'adu')
"""

sci1 = fits.getdata('E:/T120/images/CCD Image 2376.fits')
sci1 = CCDData(sci1, unit = 'adu')


#BiasSubDark = ccdproc.subtract_bias(dark1, MasterBias)
BiasSubFlat = ccdproc.subtract_bias(MasterFlat, MasterBias)
BiasSubSci1 = ccdproc.subtract_bias(sci1, MasterBias)

#DarkSubFlat = ccdproc.subtract_dark(BiasSubFlat, BiasSubDark,
#                                    dark_exposure = (3600.0 * u.second),
#                                    data_exposure = (40.0 * u.second),
#                                    scale = True)


#DarkSubSci1 = ccdproc.subtract_dark(BiasSubSci1, BiasSubDark,
#                                    dark_exposure = (3600.0 * u.second),
#                                    data_exposure = (35.0 * u.second),
#                                    scale = True)


FinalSci1 = ccdproc.flat_correct(BiasSubSci1, BiasSubFlat)

FinalSci1 = np.asarray(FinalSci1)
#FinalSci1 = (FinalSci1/2.0)
FinalSci1 = CCDData(FinalSci1, unit = 'adu')

FinalSci1 = ccdproc.cosmicray_lacosmic(FinalSci1, sigclip = 3.0)
FinalSci1.write('E:/T120/images/XPersei_final3.fits')
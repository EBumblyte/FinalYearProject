"""
Created on Wed Oct 25 15:14:33 2023

@author: Emilija Bumblyte(20431382)

Description:
    
    
"""

from astropy.io import fits
from astropy import units as u
from matplotlib import pyplot as plt
from astropy.visualization import quantity_support
from specutils import Spectrum1D
quantity_support() 
import numpy as np


filename = 'E:/SpectraXPersei/ADP_2021-04-14T141901_831.fits'
#filename = 'E:/SpectraXPersei/Standard/ADP.2020-06-19T09 50 48.327.fits'
#be star from another instrument 550 to 900nm
filename2 = 'E:/SpectraXPersei/Standard/ADP.2021-04-14T14 00_52.807.fits'
filename3 = 'E:/SpectraXPersei/Standard/ADP.2016-09-27T07 02 46.529.fits'
#filename3 = 'E:/SpectraXPersei/Standard/ADP.2020-06-19T09 50 48.356.fits'


with fits.open(filename) as f:  

    specdata = f[1].data  

with fits.open(filename2) as f2:  

    specdata2 = f2[1].data      
    
 
with fits.open(filename3) as f3:  

    specdata3 = f3[1].data      





lamb = specdata['WAVE_AIR'][0] * u.AA 

flux = specdata['FLUX'][0]  * u.Unit('erg cm-2 s-1 AA-1') 
s_n = specdata['SNR'][0] 
err = specdata['ERR'][0]

spec = Spectrum1D(spectral_axis=lamb, flux=flux, radial_velocity = -50*u.Unit("km/s")) 



lamb2 = specdata2['WAVE'][0] * u.AA 

flux2 = specdata2['FLUX'][0] /2.58173e-9 * u.Unit('erg cm-2 s-1 AA-1') 

spec2 = Spectrum1D(spectral_axis=lamb2, flux=flux2, radial_velocity = 18.5*u.Unit("km/s")) 



lamb3 = specdata3['WAVE'][0] * u.AA 

flux3 = specdata3['FLUX'][0]  *3 * u.Unit('erg cm-2 s-1 AA-1') 

spec3 = Spectrum1D(spectral_axis=lamb3, flux=flux3, radial_velocity = 32*u.Unit("km/s")) 




fig, ax = plt.subplots()

ax1= plt.subplot(211)
ax1.plot(spec.spectral_axis, spec.flux)
ax1.plot(spec2.spectral_axis, spec2.flux)
#ax1.plot(spec3.spectral_axis, spec3.flux)

ax2 = plt.subplot(212)
ax2.plot(spec.spectral_axis, spec.flux)
ax2.plot(spec3.spectral_axis, spec3.flux)


#plot of modeled spectra and observed spectra
fig, axs = plt.subplots(1, 1, figsize=(8, 6))
plt.plot(spec.spectral_axis, spec.flux)

plt.xlabel('Wavelength ($\AA$)', fontsize = 20)
plt.ylabel('Flux ($erg \ cm^{-2} \  s^{-1} \  \AA^{-1}$)', fontsize = 20)
plt.title('Spectra of X Persei', fontsize = 20)
plt.xticks(size = 20)
plt.yticks(size = 20)

axs.grid(linestyle='--', alpha = 0.6)
fig, axs = plt.subplots(1, 1, figsize=(8, 6))
plt.plot(spec.spectral_axis, err)

plt.tight_layout()
plt.show()

"""

plt.plot(spec.spectral_axis, spec.flux)
plt.plot(spec2.spectral_axis, spec2.flux)
plt.plot(spec3.spectral_axis, spec3.flux)
"""
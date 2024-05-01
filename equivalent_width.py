"""
Created on Tue Feb  6 17:20:34 2024

@author: Emilija Bumblyte(20431382)

Description:
    
    
"""
from astropy.io import fits
from astropy import units as u
from matplotlib import pyplot as plt
from astropy.visualization import quantity_support
from specutils.spectra import SpectralRegion
from specutils import analysis
from specutils.manipulation import extract_region
from specutils import Spectrum1D


quantity_support() 
filename = 'E:/SpectraXPersei/ADP_2021-04-14T141901_831.fits'

with fits.open(filename) as f:  
    specdata = f[1].data  
    
lamb = specdata['WAVE'][0]  * u.AA 
flux = specdata['FLUX'][0] * u.Unit('erg cm-2 s-1 AA-1') 
s_n = specdata['SNR'][0] 

spec = Spectrum1D(spectral_axis=lamb, flux=flux, radial_velocity = -50*u.Unit("km/s")) 

reg = SpectralRegion(6556*u.AA ,6572*u.AA)

Eq_wd = analysis.equivalent_width(spec, continuum=2.88e-12 , regions=reg)
print(Eq_wd)


red_spec = extract_region(spec, reg)

area_flux = red_spec.flux
area_wave = red_spec.wavelength


area = ((area_flux[0])) + (area_flux[1460])
area_c = 2* 2.88e-12
j = 0
dx = (area_wave[0] - area_wave[len(area_wave)-2])/(len(area_wave)-1)

for i in range((len(area_wave))):
    j = j + 1
    if i > 0 and i < 1460:
       
       if j % 2 == 0:
           area = area + (2 * ( area_flux[i]))
           area_c += 2* 2.88e-12
            
       else:
           area = area +(4 * (area_flux[i]))
           area_c += 4*2.88e-12
            
            
area = (dx/3) * area /  u.Unit('erg cm-2 s-1 AA-1') 
area_c = (dx/3) * area_c
       

print(area/(2.88e-12) - area_c/2.88e-12)


fig, ax = plt.subplots()
plt.plot(spec.spectral_axis, spec.flux)
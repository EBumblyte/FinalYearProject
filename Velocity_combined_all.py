"""
Created on Mon Nov 27 14:15:07 2023

@author: Emilija Bumblyte(20431382)

Description:
    
    
"""
from astropy.io import fits
from astropy import units as u
import numpy as np
from matplotlib import pyplot as plt
from astropy.visualization import quantity_support
from specutils import Spectrum1D
quantity_support() 


filename = 'E:/SpectraXPersei/ADP_2021-04-14T141901_831.fits'

with fits.open(filename) as f:  
    specdata = f[1].data  
    

#---------------------------------H alpha--------------------------------------

lamb_Ha = specdata['WAVE'][0] * u.AA 

flux_Ha_norm = (specdata['FLUX'][0] -2.963e-12) /1.3208e-11 * u.Unit('erg cm-2 s-1 AA-1') 
flux_Ha = (specdata['FLUX'][0] -2.963e-12)  * u.Unit('erg cm-2 s-1 AA-1') 

spec_Ha = Spectrum1D(spectral_axis=lamb_Ha, flux=flux_Ha,
                  radial_velocity = -50*u.Unit('km/s'),
                  velocity_convention = 'optical',
                  rest_value = 6562.81 * u.AA )


#----------------------------------H beta--------------------------------------

lamb_Hg = specdata['WAVE'][0] * u.AA 

flux_Hg_norm = (specdata['FLUX'][0] - 5.957e-12) /2.236e-12 * u.Unit('erg cm-2 s-1 AA-1') 
flux_Hg = (specdata['FLUX'][0] - 5.957e-12) * u.Unit('erg cm-2 s-1 AA-1') 

spec_Hg = Spectrum1D(spectral_axis=lamb_Hg, flux=flux_Hg,
                  radial_velocity = -50*u.Unit('km/s'),
                  velocity_convention = 'optical',
                  rest_value = 4340.47 * u.AA )


#---------------------------------H gamma--------------------------------------

lamb_Hb = specdata['WAVE'][0] * u.AA 

flux_Hb_norm = (specdata['FLUX'][0] -5.035e-12) /5.412e-12 * u.Unit('erg cm-2 s-1 AA-1') 
flux_Hb = (specdata['FLUX'][0] -5.035e-12)  * u.Unit('erg cm-2 s-1 AA-1') 

spec_Hb = Spectrum1D(spectral_axis=lamb_Hb, flux=flux_Hb,
                  radial_velocity = -50*u.Unit('km/s'),
                  velocity_convention = 'optical',
                  rest_value = 4861.34 * u.AA )


#----------------------------------Helium--------------------------------------

lamb_He = specdata['WAVE'][0] * u.AA 

flux_He_norm = (specdata['FLUX'][0] -3.802e-12) / 1.54699e-12 * u.Unit('erg cm-2 s-1 AA-1') 
flux_He = (specdata['FLUX'][0] -3.802e-12) * u.Unit('erg cm-2 s-1 AA-1') 

spec_He = Spectrum1D(spectral_axis=lamb_He, flux=flux_He,
                  radial_velocity = -50*u.Unit('km/s'),
                  velocity_convention = 'optical',
                  rest_value = 5875.618 * u.AA )


#----------------------------------0 line--------------------------------------

y_line1 = np.linspace(0, 1.65e-11, 100)
y_line2 = np.linspace(0, 1.05, 100)
x_line=[]
for i in range(len(y_line2)):
    x_line.append(0)


#------------------------------Plot Stuff--------------------------------------

fig, ax = plt.subplots()

ax1= plt.subplot(121)
ax1.set_xlim(-350*u.Unit('km/s'), 350*u.Unit('km/s'))
ax1.set_ylim(0, 1.4e-11)
plt.grid(True, linestyle = '--')
ax1.plot(x_line, y_line1, color = '#808080', linestyle = '--')
ax1.plot(spec_Ha.velocity, spec_Ha.flux, label = 'H\u03B1')
ax1.plot(spec_Hb.velocity, spec_Hb.flux, label = 'H\u03B2')
ax1.plot(spec_Hg.velocity, spec_Hg.flux,label = 'H\u03B3')
ax1.plot(spec_He.velocity, spec_He.flux, label = 'He')
plt.xlabel('Velocity (km/s)', fontsize = 15)
plt.ylabel('Flux ($erg \ cm^{-2} \  s^{-1} \  \AA^{-1}$)', fontsize = 15)
plt.legend(fontsize=15)


ax2= plt.subplot(122)

ax2.set_xlim(-350*u.Unit('km/s'), 350*u.Unit('km/s'))
ax2.set_ylim(0,1.05)
plt.grid(True, linestyle = '--')
ax2.plot(x_line, y_line2, color = '#808080', linestyle = '--')
plt.xlabel('Velocity (km/s)', fontsize = 15)
plt.ylabel('Normalised Flux', fontsize=15)


#---------------------------Visibility of Lines--------------------------------

line1, = ax2.plot(spec_Ha.velocity, flux_Ha_norm, label = 'H\u03B1')
line2, = ax2.plot(spec_Hb.velocity, flux_Hb_norm, label = 'H\u03B2')
line3, = ax2.plot(spec_Hg.velocity, flux_Hg_norm, label = 'H\u03B3')
line4, = ax2.plot(spec_He.velocity, flux_He_norm, label = 'He')

leg = ax2.legend(loc = 'upper right', fancybox = True, fontsize = 15)
leg.get_frame().set_alpha(0.4)
lines = [line1, line2, line3, line4]

lined = dict()
for legline, origline in zip(leg.get_lines(), lines):
    legline.set_picker(5)  # 5 pts tolerance
    lined[legline] = origline


def onpick(event):
    legline = event.artist
    origline = lined[legline]
    vis = not origline.get_visible()
    origline.set_visible(vis)

    if vis:
        legline.set_alpha(1.0)
    else:
        legline.set_alpha(0.2)
    fig.canvas.draw()

fig.canvas.mpl_connect('pick_event', onpick)

plt.show()




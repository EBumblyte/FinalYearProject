"""
Created on Mon Feb 26 16:16:47 2024

@author: Emilija Bumblyte(20431382)

Description:
    
    
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from astropy.io import fits

#----------------------import the h alpha spectra------------------------------

#set the h alpha wavelength and wavelength limits
h_alpha = 6562.8
h_alpha_up = h_alpha + 10
h_alpha_down = h_alpha - 10

#file location and name
filename = 'E:/SpectraXPersei/ADP_2021-04-14T141901_831.fits'


hdu_list = fits.open(filename)
hdu_list.info()

spectra = hdu_list[1].data

wave_temp = spectra.WAVE[0]
flux_temp = spectra.FLUX[0]
y_err = spectra.ERR[0]

spectra_flux = []
spectra_wave = []
spectra_error = []

#obtain the wavelength and flux values for range needed from spectral data
for i in range(len(wave_temp)):
    if wave_temp[i] > (h_alpha_down) and wave_temp[i] < (h_alpha_up) :
        spectra_wave.append(wave_temp[i]-1.09)
        spectra_flux.append(((flux_temp[i]  )-0.28*1e-11)*1e11*9.7*0.30*0.64*0.77*6.64*0.62) #*0.74
        spectra_error.append(y_err[i]*1e11)

#------------------------------------------------------------------------------

#----------------------set up the values for parameters------------------------

Rsol = 6.957e8
Msol = 1.989e30
Rper = 6.5 *Rsol
Mper = 15.5*Msol
G = 6.674e-11
c = 3e8

kb = 1.38e10-23
tk = 0.6*29500
mu = 2.35
ma = 1.674e-27

#variable parameters
rho0 = 1.9e-14
rho0 = 1 * 10**(-12.8)
#rho0 = 1*10 **(-13.4)
 
n = 3
#n = 3.05
#n = 2.5
#-----------------------------------------------------------------------------

#--------------------set up the spatial arrays---------------------------------
#resolution andthe radius of hte disk
nn = 101
grid_scale = 14*Rper

# set up a meshgrid
x = np.linspace(-grid_scale, grid_scale, nn)
y = np.linspace(-grid_scale, grid_scale, nn)
z = np.linspace(-grid_scale, grid_scale, nn)

# full coordinate arrays
xx, yy, zz = np.meshgrid(x, y, z)


# determine radius of a given cell in the grid
rr = (np.sqrt(xx**2+yy**2))

#------------------------------------------------------------------------------

#--------------------------calculate density----------------------------------- 

#calculate extra terms within density equation
cs = np.sqrt((kb*tk)/(mu*ma))
H0 = (cs*Rper)/673013
H0 = cs* Rper / np.sqrt((G*Mper)/(Rper)) 

# get density in the cell 
dd = rho0 *(np.sqrt(xx**2+yy**2)/Rper)**-n * np.exp(-zz**2 / (2*H0*(np.sqrt(xx**2+yy**2)/Rper)**3/2))


#angles = np.pi - (np.pi/2)*(1 + np.sign(xx))*(1-np.sign(yy**2)) - \
#         (np.pi/4)*(2+np.sign(xx))*np.sign(yy) - \
#         np.sign(xx*yy)*np.arctan((abs(xx)-abs(yy))/(abs(xx)+abs(yy)))
         
angles = np.arctan2(yy,xx)
for i in range(nn):
    for j in range(nn):
        for k in range(nn):
            if angles[i][j][k] <= 0:
                angles[i][j][k] += 2*np.pi
        

#remove matter from region in the star
null_area = 1.1* Rper
dd[rr < null_area] = 0


#------------------------------------------------------------------------------

#---------------------generate extra expelled matter---------------------------


#adding ring structures

x0, y0 = 0*Rper, 0*Rper
r = np.sqrt((xx - x0)**2 + (yy - y0)**2 )
r_full = np.zeros_like(r)

#set limits for rings
inner_rad = [10*Rper, 5*Rper, 1.5*Rper, 2.1*Rper]
out_rad = [14*Rper, 8.5*Rper, 3.5*Rper, 2.4*Rper]
inner_rad = [2.9*Rper]
out_rad = [4.1*Rper]

#inverting limits 
inner_rad_inv = [1.2*Rper, 1.5*Rper, 5*Rper, 10*Rper]
out_rad_inv = [1.4*Rper, 3.5*Rper, 8.5*Rper, 14*Rper]
inner_rad_inv = [ 4*Rper]
out_rad_inv = [ 4.8*Rper]

for i in range(len(inner_rad)):
    r = np.sqrt((xx - x0)**2 + (yy - y0)**2 )
    r[r>out_rad[i]] = 0
    r[r < inner_rad[i]] = 0
    
    r_full += r
    

#adding extra matter in angle to the rings

theta = [0.42*np.pi,0.57*np.pi]
theta2 = [0*np.pi,0.15*np.pi]
theta3=[1.65*np.pi,2.0*np.pi]
theta4 = [1.84*np.pi, 1.95*np.pi]

theta5 = [0.485*np.pi, 0.515*np.pi]

for i in range(nn):
    for j in range(nn):
        for k in range(nn):
            if angles[i][j][k] > theta[0] and angles[i][j][k] < theta[1]:
                r_full[i][j][k] *= 2.5
                
            if angles[i][j][k] >= theta2[0] and angles[i][j][k] < theta2[1]:
                r_full[i][j][k] *= 11
                
            if angles[i][j][k] > theta3[0] and angles[i][j][k] < theta3[1]:
                r_full[i][j][k] *= 18.0

            if angles[i][j][k] > theta4[0] and angles[i][j][k] < theta4[1]:
                r_full[i][j][k] *= 2.0
                
            if angles[i][j][k] > theta5[0] and angles[i][j][k] < theta5[1]:
                r_full[i][j][k] *= 8


for i in range(len(x)):
    for j in range(len(y)):
        for k in range(len(z)):
       
            if r_full[i][j][k] > 0:
                dd[i][j][k] = dd[i][j][k] * (r_full[i][j][k] /1e10 *1 )


#adding ring structures

x0, y0 = 0*Rper, 0*Rper
r = np.sqrt((xx - x0)**2 + (yy - y0)**2 )
r_full = np.zeros_like(r)

#set limits for rings
inner_rad = [10*Rper, 5*Rper, 1.5*Rper, 2.1*Rper]
out_rad = [14*Rper, 8.5*Rper, 3.5*Rper, 2.4*Rper]
inner_rad = [1.3*Rper]
out_rad = [2.3*Rper]


for i in range(len(inner_rad)):
    r = np.sqrt((xx - x0)**2 + (yy - y0)**2 )
    r[r>out_rad[i]] = 0
    r[r < inner_rad[i]] = 0
    
    r_full += r
    
theta = [0.42*np.pi,0.5*np.pi]

for i in range(nn):
    for j in range(nn):
        for k in range(nn):
            if angles[i][j][k] > theta[0] and angles[i][j][k] < theta[1]:
                r_full[i][j][k] *= 10
    
for i in range(len(x)):
    for j in range(len(y)):
        for k in range(len(z)):
       
            if r_full[i][j][k] > 0:
                dd[i][j][k] = dd[i][j][k] * (r_full[i][j][k] /1e10 *1 )

"""

#adding spherical extra matter

rad = (0, grid_scale, nn)
circle_rad = (1*Rper, 5*Rper, nn)


#radius of the sphere for expelled matter
radius = 2.5*Rper
#set position of the extra matter
x0, y0, z0 = 0*Rper, 4.5*Rper, 0*Rper

#set the parameters for the spheroid
r = np.sqrt((xx - x0)**2 + (yy - y0)**2 + ((zz - z0)**2)/10)
r[r > radius] = 0


for i in range(len(x)):
    for j in range(len(y)):
        for k in range(len(z)):
            if r[i][j][k] > 0:
                dd[i][j][k] = dd[i][j][k] *6.5
                
       
                

#radius of the sphere for expelled matter
radius = 5*Rper
#set position of the extra matter
x0, y0, z0 = 8.5*Rper, -3.2*Rper, 0*Rper

#set the parameters for the spheroid
r = np.sqrt((xx - x0)**2 + (yy - y0)**2 + ((zz - z0)**2)/10)
r[r > radius] = 0


for i in range(len(x)):
    for j in range(len(y)):
        for k in range(len(z)):
            if r[i][j][k] > 0:
                dd[i][j][k] = dd[i][j][k] *31
                
                
                
#radius of the sphere for expelled matter
radius =1.9*Rper
#set position of the extra matter
x0, y0, z0 = 5.4*Rper, 2.8*Rper, 0*Rper

#set the parameters for the spheroid
r = np.sqrt((xx - x0)**2 + (yy - y0)**2 + ((zz - z0)**2)/10)
r[r > radius] = 0


for i in range(len(x)):
    for j in range(len(y)):
        for k in range(len(z)):
            if r[i][j][k] > 0:
                dd[i][j][k] = dd[i][j][k] *20
              
                
#radius of the sphere for expelled matter
radius = 1*Rper
#set position of the extra matter
x0, y0, z0 = 8.5*Rper, -3.2*Rper, 0*Rper

#set the parameters for the spheroid
r = np.sqrt((xx - x0)**2 + (yy - y0)**2 + ((zz - z0)**2)/10)
r[r > radius] = 0


for i in range(len(x)):
    for j in range(len(y)):
        for k in range(len(z)):
            if r[i][j][k] > 0:
                dd[i][j][k] = dd[i][j][k] *6

"""


#------------------------------------------------------------------------------

#-----------------------------generate arms------------------------------------
"""
#adding arm 
theta = (0, 2*np.pi, nn)
#rad_place = (5*Rper,  23*Rper, nn)
#sigma = (0.1*Rper, 1*Rper, nn)

rad = (0, grid_scale, nn)
circle_rad = (1*Rper, 5*Rper, nn)
extra_dd = np.zeros_like(angles)

Oout = np.pi
rin  = 2*Rper
rout = 16*Rper


fr =  np.log(1 + rr**2) + np.cos(np.arctan2(xx, yy) + np.log(1 + rr**2))
frr =   1e12* (((fr)**(-1)))**7  *100

dd2 = 100*dd *  (frr)

null_area = 1.1* Rper
dd2[rr < null_area] = 0
"""

"""

for i in range(len(theta)):

    for j in range(len(x)):
        for k in range( len(y)):
            for l in range(len(z)):
                
                while count == 0:
                    if theta[i] == angles[j][k][l]:
                        radius = circle_rad[i]
                        x0 = xx[j]
                        z0 = zz[l]
                        y0 = yy[k]
                        #set the parameters for the spheroid
                        circle = np.sqrt((xx - x0)**2 + ((yy - y0)**2)/25 + ((zz - z0)**2))
                        circle[circle > radius] = 0
                        dd += circle
                        count = 1
"""
#------------------------------------------------------------------------------

#--------------find velocities abnd shift and put them in bins-----------------
#h alpha in si units
h_alpha_si = 6562.8e-10


vel = (np.sqrt((G*Mper/(np.sqrt(xx**2+yy**2+zz**2)))) * np.cos(np.pi/2 + angles)) * np.sin(np.pi/6)

shift = ((vel * h_alpha_si)/c) + h_alpha_si

#define the shift range to examine and number of bins 
shift_range = np.linspace(6552.8e-10, 6572.8e-10, 150)
#135
#create array to have corresponding flux values for each wavelength range
distribution = np.full_like(shift_range, 0)


total = len(shift_range-1)
for l in range(len(shift_range)-1):
    
    #progress bar for the for loop below
    progress = l/total
    arrow = '=' * int(30 * progress)
    spaces = ' ' * (30 - len(arrow))
    print(f'\r[{arrow}{spaces}] {int(progress * 100)}%', end='',flush=True)
    
    
    #for loop to sort the velocity 
    for i in range(nn):
        for j in range(nn):
            for k in range(nn):
                if shift[i][j][k] >= shift_range[l] and shift[i][j][k] < shift_range[l+1]:
                   distribution[l] = distribution[l] + (dd[i][j][k] )
                   
distribution *= 1e11
            
new_binning = np.linspace(6552.8, 6572.8, len(shift_range))
new_bin_flux = np.full_like(new_binning, 0)
new_bin_error = np.full_like(new_binning,0)
for i in range(len(new_binning)-1):
    for k in range(len(spectra_wave)):
        if spectra_wave[k] >= new_binning[i] and spectra_wave[k] <= new_binning[i+1]:
            new_bin_flux[i] += spectra_flux[k]
            new_bin_error += spectra_error[i]

new_bin_flux /= 9.7 

chi_squared = 0
for i in range(len(new_bin_flux)):
    if new_bin_flux[i] > 0:
        chi_squared += ((distribution[i]-new_bin_flux[i])**2/(new_bin_flux[i]))
        
        
def calc_reduced_chi_square(fit, x, y, yerr, N, n_free):
    '''
    fit (array) values for the fit
    x,y,yerr (arrays) data
    N total number of points
    n_free number of parameters we are fitting
    '''
    return 1.0/(N-n_free)*sum(((fit - y)/yerr)**2)

chisquared_red = calc_reduced_chi_square(distribution, new_binning, new_bin_flux, new_bin_error, len(distribution), 6)
    
#-------------------------velocity values for contour plot---------------------        

xx_vel ,yy_vel = np.meshgrid(x,y)

angles_eq = np.pi - (np.pi/2)*(1 + np.sign(xx_vel))*(1-np.sign(yy_vel**2)) - \
    (np.pi/4)*(2+np.sign(xx_vel))*np.sign(yy_vel) - \
        np.sign(xx_vel*yy_vel)*np.arctan((abs(xx_vel)-abs(yy_vel))/(abs(xx_vel)+abs(yy_vel)))
        
vel_eq = ((np.sqrt((G*Mper/(np.sqrt(xx_vel**2+yy_vel**2)))) * np.cos(np.pi/2 + angles_eq)) * np.sin(np.pi/6))/1000

#------------------------------------------------------------------------------


#----------------------------plotting------------------------------------------

# top down view by summing along z axis (max normalised to 1)
fig, axs = plt.subplots(1, 1, figsize=(6, 6))
ii = axs.imshow(np.nansum(dd, axis=2)/np.amax(np.nansum(dd, axis=2)), norm=LogNorm(vmin=0.0001, vmax=1), origin='lower')
plt.gca().set_aspect('equal')
fig.colorbar(ii, ax=axs, orientation='vertical', fraction=.045)
axs.set_xlabel('x-axis', fontsize = 15)
axs.set_ylabel('y-axis', fontsize = 15)
plt.tight_layout()
plt.show()

"""
# top down view by summing along x axis (max normalised to 1)
#note: when looking at the axes plotted, horizontal axis appears to be z and vertical apears to be y
fig, axs = plt.subplots(1, 1, figsize=(8, 6))
ii = axs.imshow(np.nansum(dd, axis=1)/np.amax(np.nansum(dd, axis=1)), norm=LogNorm(vmin=0.0001, vmax=1), origin='lower')
fig.colorbar(ii, ax=axs, orientation='vertical', fraction=.05)
axs.set_xlabel('x-axis/y-axis', fontsize = 15)
axs.set_ylabel('z-axis', fontsize = 15)
plt.tight_layout()
plt.show()


fig, axs = plt.subplots(1, 1, figsize=(8, 6))
ii = axs.imshow(dd[64,:,:])
fig.colorbar(ii, ax=axs, orientation='vertical', fraction=.05)
axs.set_xlabel('x-axis/y-axis', fontsize = 15)
axs.set_ylabel('z-axis', fontsize = 15)
plt.tight_layout()
plt.show()
"""
"""
fig, axs = plt.subplots(1, 1, figsize=(8, 6))
ii = axs.imshow(frr[:,:,60])
fig.colorbar(ii, ax=axs, orientation='vertical', fraction=.05)
axs.set_xlabel('x-axis/y-axis', fontsize = 15)
axs.set_ylabel('z-axis', fontsize = 15)
plt.tight_layout()
plt.show()
"""
# integrated density distribtion
fig, axs = plt.subplots(1, 1, figsize=(8, 6))
axs.plot(x, (np.nansum(np.nansum(dd, axis=2), axis=1)/np.amax(np.nansum(np.nansum(dd, axis=2), axis=1))))
plt.tight_layout()
plt.show()



# top down view of the velocity at equator-plot and contours 
fig, axs = plt.subplots(1, 1, figsize=(6, 6))
axs1 = axs.contourf(xx_vel,yy_vel,vel_eq, 40)
axs2 = axs.contour(axs1, colors='white', linewidth = 0.8, alpha = 0.5)
plt.clabel(axs2, inline=True, fontsize=8)
plt.gca().set_aspect('equal')
cbar =fig.colorbar(axs1)
axs.set_xlabel('x-axis', fontsize = 15)
axs.set_ylabel('y-axis', fontsize = 15)
cbar.ax.set_ylabel('velocity ($km \ s^{-1}$)')
plt.tight_layout()
plt.show()

"""
# top down view of the velocity at equator contours 
fig, axs = plt.subplots(1, 1, figsize=(6, 6))
ax = axs.imshow(vel_eq)
CS = axs.contour(vel_eq, origin='lower', linewidths=1,colors = 'white')
axs.clabel(CS,  inline=True, fontsize = 11)
plt.gca().set_aspect('equal')
cbar = fig.colorbar(ii, ax=axs, orientation='vertical', fraction=.045)
axs.set_xlabel('x-axis', fontsize = 15)
axs.set_ylabel('y-axis', fontsize = 15)
cbar.ax.set_ylabel('velocity ($km \ s^{-1}$)')
plt.tight_layout()
plt.show()
"""

#plot of modeled spectra and observed spectra
fig, axs = plt.subplots(1, 1, figsize=(8, 6))
axs.plot(shift_range* 1e10, distribution, label='Modelled Emission Line')
axs.plot(spectra_wave, spectra_flux , label = 'Emission Line Observed')
axs.set_xlabel('Wavelength ($\AA$)', fontsize = 15)
axs.set_ylabel('Flux', fontsize = 15)
plt.grid(linestyle='--', alpha=0.6)
plt.legend(fontsize=15)
plt.tight_layout()
plt.show()


#plot of modeled spectra and observed spectra
fig, axs = plt.subplots(1, 1, figsize=(8, 6))
axs.plot(shift_range* 1e10, distribution)
axs.plot(new_binning, new_bin_flux )
axs.set_xlabel('Wavelength ($\AA$)', fontsize = 15)
axs.set_ylabel('Flux', fontsize = 15)
plt.tight_layout()
plt.show()

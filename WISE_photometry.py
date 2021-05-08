


import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.visualization import ImageNormalize, ZScaleInterval, PercentileInterval
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy import units as u
from astropy.coordinates import SkyCoord

from photutils import SkyCircularAperture
from photutils import CircularAperture
from photutils import CircularAnnulus
from photutils import aperture_photometry

def wise1(image_data, radpix, wcs):
    #aperture = CircularAperture((image_data.shape[0]/2,image_data.shape[0]/2), radpix[19])
    position = [image_data.shape[0]/2,image_data.shape[0]/2]
    apertures = [CircularAperture(position, r=r) for r in radpix]
    annulus_aperture = CircularAnnulus(position, r_in=radpix[19], r_out=radpix[19]+radpix[16])
    phot_table = aperture_photometry(image_data, apertures, wcs=wcs)
    bkg_table = aperture_photometry(image_data, annulus_aperture, wcs=wcs)   
    bkg_mean = bkg_table['aperture_sum'] / annulus_aperture.area()
    #sky_aperture = to_sky(apertures,wcs)
    phot_array = np.zeros(20)
    bkg_sum_array = np.zeros(20)
    for i in range(0,20):
        phot_array[i] = phot_table['aperture_sum_' + str(i)][0]
        bkg_sum_array[i] = bkg_mean * apertures[i].area()
        phot_table['aperture_sum_' + str(i)].info.format = '%.8g'  # for consistent table output
    
    final_sum =  phot_array - bkg_sum_array   
    
    wisemags = np.array([final_sum, 0, 0, 0])
    coeffs_interp = get_coeffs(*wisemags)
    f1, f2, f3, f4 = Fsnu0/coeffs_interp * 10**(-wisemags/2.5)
    print ('Corrected WISE fluxes: %g, %g, %g, %g'%(f1, f2, f3, f4)
    phot_Jy_wise1 = f1
    
    bkg_table ['aperture_sum'].info.format = '%.8g'  # for consistent table output
    print('Backgorud outer radius =',radpix[19]+radpix[16],'pixels')
    print(bkg_table)
    print(phot_table)
    return(phot_Jy_wise1) 

def wise2(image_data, radpix, wcs):
    #aperture = CircularAperture((image_data.shape[0]/2,image_data.shape[0]/2), radpix[19])
    position = [image_data.shape[0]/2,image_data.shape[0]/2]
    apertures = [CircularAperture(position, r=r) for r in radpix]
    annulus_aperture = CircularAnnulus(position, r_in=radpix[19], r_out=radpix[19]+radpix[16])
    phot_table = aperture_photometry(image_data, apertures, wcs=wcs)
    bkg_table = aperture_photometry(image_data, annulus_aperture, wcs=wcs)   
    bkg_mean = bkg_table['aperture_sum'] / annulus_aperture.area()
    #sky_aperture = to_sky(apertures,wcs)
    phot_array = np.zeros(20)
    bkg_sum_array = np.zeros(20)
    for i in range(0,20):
        phot_array[i] = phot_table['aperture_sum_' + str(i)][0]
        bkg_sum_array[i] = bkg_mean * apertures[i].area()
        phot_table['aperture_sum_' + str(i)].info.format = '%.8g'  # for consistent table output
    
    final_sum =  phot_array - bkg_sum_array   
    
    wisemags = np.array([0, final_sum, 0, 0])
    coeffs_interp = get_coeffs(*wisemags)
    f1, f2, f3, f4 = Fsnu0/coeffs_interp * 10**(-wisemags/2.5)
    print ('Corrected WISE fluxes: %g, %g, %g, %g'%(f1, f2, f3, f4)
    phot_Jy_wise1 = f2
    
    bkg_table ['aperture_sum'].info.format = '%.8g'  # for consistent table output
    print('Backgorud outer radius =',radpix[19]+radpix[16],'pixels')
    print(bkg_table)
    print(phot_table)
    return(phot_Jy_wise2) 

def wise3(image_data, radpix, wcs):
    #aperture = CircularAperture((image_data.shape[0]/2,image_data.shape[0]/2), radpix[19])
    position = [image_data.shape[0]/2,image_data.shape[0]/2]
    apertures = [CircularAperture(position, r=r) for r in radpix]
    annulus_aperture = CircularAnnulus(position, r_in=radpix[19], r_out=radpix[19]+radpix[16])
    phot_table = aperture_photometry(image_data, apertures, wcs=wcs)
    bkg_table = aperture_photometry(image_data, annulus_aperture, wcs=wcs)   
    bkg_mean = bkg_table['aperture_sum'] / annulus_aperture.area()
    #sky_aperture = to_sky(apertures,wcs)
    phot_array = np.zeros(20)
    bkg_sum_array = np.zeros(20)
    for i in range(0,20):
        phot_array[i] = phot_table['aperture_sum_' + str(i)][0]
        bkg_sum_array[i] = bkg_mean * apertures[i].area()
        phot_table['aperture_sum_' + str(i)].info.format = '%.8g'  # for consistent table output
    
    final_sum =  phot_array - bkg_sum_array   
    
    wisemags = np.array([0, 0, final_sum, 0])
    coeffs_interp = get_coeffs(*wisemags)
    f1, f2, f3, f4 = Fsnu0/coeffs_interp * 10**(-wisemags/2.5)
    print ('Corrected WISE fluxes: %g, %g, %g, %g'%(f1, f2, f3, f4)
    phot_Jy_wise1 = f3
    
    bkg_table ['aperture_sum'].info.format = '%.8g'  # for consistent table output
    print('Backgorud outer radius =',radpix[19]+radpix[16],'pixels')
    print(bkg_table)
    print(phot_table)
    return(phot_Jy_wise3) 

def wise4(image_data, radpix, wcs):
    #aperture = CircularAperture((image_data.shape[0]/2,image_data.shape[0]/2), radpix[19])
    position = [image_data.shape[0]/2,image_data.shape[0]/2]
    apertures = [CircularAperture(position, r=r) for r in radpix]
    annulus_aperture = CircularAnnulus(position, r_in=radpix[19], r_out=radpix[19]+radpix[16])
    phot_table = aperture_photometry(image_data, apertures, wcs=wcs)
    bkg_table = aperture_photometry(image_data, annulus_aperture, wcs=wcs)   
    bkg_mean = bkg_table['aperture_sum'] / annulus_aperture.area()
    #sky_aperture = to_sky(apertures,wcs)
    phot_array = np.zeros(20)
    bkg_sum_array = np.zeros(20)
    for i in range(0,20):
        phot_array[i] = phot_table['aperture_sum_' + str(i)][0]
        bkg_sum_array[i] = bkg_mean * apertures[i].area()
        phot_table['aperture_sum_' + str(i)].info.format = '%.8g'  # for consistent table output
    
    final_sum =  phot_array - bkg_sum_array   
    
    wisemags = np.array([0, 0, 0, final_sum])
    coeffs_interp = get_coeffs(*wisemags)
    f1, f2, f3, f4 = Fsnu0/coeffs_interp * 10**(-wisemags/2.5)
    print ('Corrected WISE fluxes: %g, %g, %g, %g'%(f1, f2, f3, f4)
    phot_Jy_wise1 = f4
    
    bkg_table ['aperture_sum'].info.format = '%.8g'  # for consistent table output
    print('Backgorud outer radius =',radpix[19]+radpix[16],'pixels')
    print(bkg_table)
    print(phot_table)
    return(phot_Jy_wise4) 

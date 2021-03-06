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
from astropy.stats import sigma_clipped_stats

def photometry_MSX_A(image_data, radpix, wcs):
    #aperture = CircularAperture((image_data.shape[0]/2,image_data.shape[0]/2), radpix[19])
    position = [image_data.shape[0]/2,image_data.shape[0]/2]
    apertures = [CircularAperture(position, r=r) for r in radpix]
    annulus_aperture = CircularAnnulus(position, r_in=radpix[19], r_out=radpix[19]+radpix[16])
    phot_table = aperture_photometry(image_data, apertures, wcs=wcs)

    annulus_masks = annulus_aperture.to_mask(method='center')
    annulus_data = annulus_masks.multiply(image_data)
    mask = annulus_masks.data
    annulus_data_1d = annulus_data[mask > 0]
    _, median_sigclip, _ = sigma_clipped_stats(annulus_data_1d)
    bkg_mean = median_sigclip

    #sky_aperture = to_sky(apertures,wcs)
    phot_array = np.zeros(20)
    bkg_sum_array = np.zeros(20)
    for i in range(0,20):
        phot_array[i] = phot_table['aperture_sum_' + str(i)][0]
        bkg_sum_array[i] = bkg_mean * apertures[i].area

    final_sum =  phot_array - bkg_sum_array

    scale = np.mean(proj_plane_pixel_scales(wcs))
    phot_Jy_MSX_A = final_sum* 7.133e12 * (scale*np.pi/180.0)**2 # MSX_A W/m^2-sr to Jy *7.133e12*scale**2/((180/np.pi)**2*3600**2)

    print('Backgorud outer radius =',radpix[19]+radpix[16],'pixels')
    print(phot_table)
    return(phot_Jy_MSX_A)

def photometry_MSX_C(image_data, radpix, wcs):
    #aperture = CircularAperture((image_data.shape[0]/2,image_data.shape[0]/2), radpix[19])
    position = [image_data.shape[0]/2,image_data.shape[0]/2]
    apertures = [CircularAperture(position, r=r) for r in radpix]
    annulus_aperture = CircularAnnulus(position, r_in=radpix[19], r_out=radpix[19]+radpix[16])
    phot_table = aperture_photometry(image_data, apertures, wcs=wcs)

    annulus_masks = annulus_aperture.to_mask(method='center')
    annulus_data = annulus_masks.multiply(image_data)
    mask = annulus_masks.data
    annulus_data_1d = annulus_data[mask > 0]
    _, median_sigclip, _ = sigma_clipped_stats(annulus_data_1d)
    bkg_mean = median_sigclip

    #sky_aperture = to_sky(apertures,wcs)
    phot_array = np.zeros(20)
    bkg_sum_array = np.zeros(20)
    for i in range(0,20):
        phot_array[i] = phot_table['aperture_sum_' + str(i)][0]
        bkg_sum_array[i] = bkg_mean * apertures[i].area

    final_sum =  phot_array - bkg_sum_array

    scale = np.mean(proj_plane_pixel_scales(wcs))
    phot_Jy_MSX_C = final_sum* 2.863e13 * (scale*np.pi/180.0)**2  # MSX_C W/m^2-sr to Jy *2.863e13*scale**2/((180/np.pi)**2*3600**2)


    print('Backgorud outer radius =',radpix[19]+radpix[16],'pixels')
    print(phot_table)
    return(phot_Jy_MSX_C)

def photometry_MSX_D(image_data, radpix, wcs):
    #aperture = CircularAperture((image_data.shape[0]/2,image_data.shape[0]/2), radpix[19])
    position = [image_data.shape[0]/2,image_data.shape[0]/2]
    apertures = [CircularAperture(position, r=r) for r in radpix]
    annulus_aperture = CircularAnnulus(position, r_in=radpix[19], r_out=radpix[19]+radpix[16])
    phot_table = aperture_photometry(image_data, apertures, wcs=wcs)
    annulus_masks = annulus_aperture.to_mask(method='center')
    annulus_data = annulus_masks.multiply(image_data)
    mask = annulus_masks.data
    annulus_data_1d = annulus_data[mask > 0]
    _, median_sigclip, _ = sigma_clipped_stats(annulus_data_1d)
    bkg_mean = median_sigclip

    #sky_aperture = to_sky(apertures,wcs)
    phot_array = np.zeros(20)
    bkg_sum_array = np.zeros(20)
    for i in range(0,20):
        phot_array[i] = phot_table['aperture_sum_' + str(i)][0]
        bkg_sum_array[i] = bkg_mean * apertures[i].area

    final_sum =  phot_array - bkg_sum_array

    scale = np.mean(proj_plane_pixel_scales(wcs))
    phot_Jy_MSX_D = final_sum* 3.216e13 * (scale*np.pi/180.0)**2 # MSX_D W/m^2-sr to Jy *3.216e13*scale**2/((180/np.pi)**2*3600**2)

    print('Backgorud outer radius =',radpix[19]+radpix[16],'pixels')
    print(phot_table)
    return(phot_Jy_MSX_D)


def photometry_MSX_E(image_data, radpix, wcs):
    #aperture = CircularAperture((image_data.shape[0]/2,image_data.shape[0]/2), radpix[19])
    position = [image_data.shape[0]/2,image_data.shape[0]/2]
    apertures = [CircularAperture(position, r=r) for r in radpix]
    annulus_aperture = CircularAnnulus(position, r_in=radpix[19], r_out=radpix[19]+radpix[16])
    phot_table = aperture_photometry(image_data, apertures, wcs=wcs)
    annulus_masks = annulus_aperture.to_mask(method='center')
    annulus_data = annulus_masks.multiply(image_data)
    mask = annulus_masks.data
    annulus_data_1d = annulus_data[mask > 0]
    _, median_sigclip, _ = sigma_clipped_stats(annulus_data_1d)
    bkg_mean = median_sigclip

    #sky_aperture = to_sky(apertures,wcs)
    phot_array = np.zeros(20)
    bkg_sum_array = np.zeros(20)
    for i in range(0,20):
        phot_array[i] = phot_table['aperture_sum_' + str(i)][0]
        bkg_sum_array[i] = bkg_mean * apertures[i].area

    final_sum =  phot_array - bkg_sum_array

    scale = np.mean(proj_plane_pixel_scales(wcs))
    phot_Jy_MSX_E = final_sum * 2.476e13  * (scale*np.pi/180.0)**2 # MSX_E W/m^2-sr to Jy *2.476e13*scale**2/((180/np.pi)**2*3600**2)

    print('Backgorud outer radius =',radpix[19]+radpix[16],'pixels')
    print(phot_table)
    return(phot_Jy_MSX_E)

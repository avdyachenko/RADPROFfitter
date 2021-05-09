 
def photometry_IRAS100(image_data, radpix, wcs):
    #aperture = CircularAperture((image_data.shape[0]/2,image_data.shape[0]/2), radpix[19])
    position = [image_data.shape[0]/2,image_data.shape[0]/2]
    apertures = [CircularAperture(position, r=r) for r in radpix]
    annulus_aperture = CircularAnnulus(position, r_in=radpix[19], r_out=radpix[19]+radpix[16])
    phot_table = aperture_photometry(image_data, apertures, wcs=wcs)
    bkg_table = aperture_photometry(image_data, annulus_aperture, wcs=wcs)   
    bkg_mean = bkg_table['aperture_sum'] / annulus_aperture.area
    #sky_aperture = to_sky(apertures,wcs)
    phot_array = np.zeros(20)
    bkg_sum_array = np.zeros(20)
    for i in range(0,20):
        phot_array[i] = phot_table['aperture_sum_' + str(i)][0]
        bkg_sum_array[i] = bkg_mean * apertures[i].area
        phot_table['aperture_sum_' + str(i)].info.format = '%.8g'  # for consistent table output
    
    final_sum =  phot_array - bkg_sum_array  
    
    scale = np.mean(proj_plane_pixel_scales(wcs))    
    phot_Jy_IRAS100= final_sum* (scale*np.pi/180.0)**2 # Jy/sr to Jy
    
    
    bkg_table ['aperture_sum'].info.format = '%.8g'  # for consistent table output
    print('Backgorud outer radius =',radpix[19]+radpix[16],'pixels')
    print(bkg_table)
    print(phot_table)
    return(phot_Jy_IRAS100)

def photometry_IRAS60(image_data, radpix, wcs):
    #aperture = CircularAperture((image_data.shape[0]/2,image_data.shape[0]/2), radpix[19])
    position = [image_data.shape[0]/2,image_data.shape[0]/2]
    apertures = [CircularAperture(position, r=r) for r in radpix]
    annulus_aperture = CircularAnnulus(position, r_in=radpix[19], r_out=radpix[19]+radpix[16])
    phot_table = aperture_photometry(image_data, apertures, wcs=wcs)
    bkg_table = aperture_photometry(image_data, annulus_aperture, wcs=wcs)   
    bkg_mean = bkg_table['aperture_sum'] / annulus_aperture.area
    #sky_aperture = to_sky(apertures,wcs)
    phot_array = np.zeros(20)
    bkg_sum_array = np.zeros(20)
    for i in range(0,20):
        phot_array[i] = phot_table['aperture_sum_' + str(i)][0]
        bkg_sum_array[i] = bkg_mean * apertures[i].area
        phot_table['aperture_sum_' + str(i)].info.format = '%.8g'  # for consistent table output
    
    final_sum =  phot_array - bkg_sum_array   
    
    scale = np.mean(proj_plane_pixel_scales(wcs))  
    phot_Jy_IRAS60 = final_sum* (scale*np.pi/180.0)**2 # Jy/sr to Jy 

    
    bkg_table ['aperture_sum'].info.format = '%.8g'  # for consistent table output
    print('Backgorud outer radius =',radpix[19]+radpix[16],'pixels')
    print(bkg_table)
    print(phot_table)
    return(phot_Jy_IRAS60)

def photometry_IRAS25(image_data, radpix, wcs):
    #aperture = CircularAperture((image_data.shape[0]/2,image_data.shape[0]/2), radpix[19])
    position = [image_data.shape[0]/2,image_data.shape[0]/2]
    apertures = [CircularAperture(position, r=r) for r in radpix]
    annulus_aperture = CircularAnnulus(position, r_in=radpix[19], r_out=radpix[19]+radpix[16])
    phot_table = aperture_photometry(image_data, apertures, wcs=wcs)
    bkg_table = aperture_photometry(image_data, annulus_aperture, wcs=wcs)   
    bkg_mean = bkg_table['aperture_sum'] / annulus_aperture.area
    #sky_aperture = to_sky(apertures,wcs)
    phot_array = np.zeros(20)
    bkg_sum_array = np.zeros(20)
    for i in range(0,20):
        phot_array[i] = phot_table['aperture_sum_' + str(i)][0]
        bkg_sum_array[i] = bkg_mean * apertures[i].area
        phot_table['aperture_sum_' + str(i)].info.format = '%.8g'  # for consistent table output
    
    final_sum =  phot_array - bkg_sum_array   
    
    scale = np.mean(proj_plane_pixel_scales(wcs))  
    phot_Jy_IRAS25 = final_sum* (scale*np.pi/180.0)**2 # Jy/sr to Jy 
    
    bkg_table ['aperture_sum'].info.format = '%.8g'  # for consistent table output
    print('Backgorud outer radius =',radpix[19]+radpix[16],'pixels')
    print(bkg_table)
    print(phot_table)
    return(phot_Jy_IRAS25)


def photometry_IRAS12(image_data, radpix, wcs):
    #aperture = CircularAperture((image_data.shape[0]/2,image_data.shape[0]/2), radpix[19])
    position = [image_data.shape[0]/2,image_data.shape[0]/2]
    apertures = [CircularAperture(position, r=r) for r in radpix]
    annulus_aperture = CircularAnnulus(position, r_in=radpix[19], r_out=radpix[19]+radpix[16])
    phot_table = aperture_photometry(image_data, apertures, wcs=wcs)
    bkg_table = aperture_photometry(image_data, annulus_aperture, wcs=wcs)   
    bkg_mean = bkg_table['aperture_sum'] / annulus_aperture.area
    #sky_aperture = to_sky(apertures,wcs)
    phot_array = np.zeros(20)
    bkg_sum_array = np.zeros(20)
    for i in range(0,20):
        phot_array[i] = phot_table['aperture_sum_' + str(i)][0]
        bkg_sum_array[i] = bkg_mean * apertures[i].area
        phot_table['aperture_sum_' + str(i)].info.format = '%.8g'  # for consistent table output
    
    final_sum =  phot_array - bkg_sum_array   
    
    scale = np.mean(proj_plane_pixel_scales(wcs))  
    phot_IRAS12 = final_sum * (scale*np.pi/180.0)**2 # Jy/sr to Jy 
    
    bkg_table ['aperture_sum'].info.format = '%.8g'  # for consistent table output
    print('Backgorud outer radius =',radpix[19]+radpix[16],'pixels')
    print(bkg_table)
    print(phot_table)
    return(phot_Jy_IRAS12)

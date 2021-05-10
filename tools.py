def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx], idx

def rob_fits(user_fitsfile):
    hdul = fits.open(user_fitsfile)
    wavelength = hdul['SPECTRAL_INFO'].data
    wave = wavelength['WAVELENGTH']
    aperture = hdul['APERTURES'].data
    ap = aperture['APERTURE']
    fnu = hdul['VALUES'].data
    mod_names = hdul['MODEL_NAMES'].data
    names = np.array(mod_names.astype('<U11'))
    return(ap, wave, fnu*(1000/dist)**2/1000, names)

def user_fits(rob_fitsfile,ap,flag):
    hdul = fits.open(rob_fitsfile)
    hdul.info()
    image_data = hdul[0].data[::-1,: ] #инвертирует ось y
    wcs = WCS(hdul[0].header)
    waveleng = hdul[0].header['WAVELENG']
    scale = np.mean(proj_plane_pixel_scales(wcs))*3600
    #scale = np.abs(CDELT1)*3600  #arcsec/pix
    FOV = image_data.shape[0]*scale
    hdul.close()
    print("Scale =",scale, " arcsec/pix")
    print("FOV=", FOV/60/60, "deg")
    print('Image data shape:', image_data.shape)
    # AU           | 1.49598e+11 m
    # pc           | 3.08568e+16 m
    L = dist * FOV    #мой способ!!
    aupixel = L/image_data.shape[0]   # AU in pixel
    radius = ap/aupixel/2

    radpix = ap/2/dist/scale
    print('MAX RADIUS =', radpix[19], 'pix')

    #norm = ImageNormalize(image_data, interval=ZScaleInterval())         # Scale = ZScale
    norm = ImageNormalize(image_data, interval=PercentileInterval(99.5))  # Scale = 99.5%

    ax = plt.subplot(projection=wcs)
    for i in range(0,len(radpix)):
        ax.add_patch(plt.Circle((image_data.shape[0]/2,image_data.shape[0]/2), radpix[i] ,lw = 0.1, color='r', fill=False))
        ax.add_patch(plt.Circle((image_data.shape[0]/2,image_data.shape[0]/2), radpix[19]+ radpix[16] ,lw = 0.5, color='green', fill=False))
    plt.imshow(image_data, cmap='gray',  norm=norm)
    plt.xlabel('Galactic Longitude')
    plt.ylabel('Galactic Latitude')
    plt.colorbar()
    plt.savefig(str(obj_name +'_' + models +'_' + flag +'.png') , dpi=300)
    print('Min:', np.min(image_data))
    print('Max:', np.max(image_data))
    print('Mean:', np.mean(image_data))
    print('Stdev:', np.std(image_data))
    return(image_data, radpix, wcs, scale, waveleng,  aupixel )


def chi2_func(phot_Jy, fnu, wave_index, models,names,flag):
    chi2 = np.zeros(fnu.shape[0])
    for i in range(0,fnu.shape[0]):
        #chi2[i] = ((phot_Jy_MSX_A - fnu[i,:,wave_index])**2 *radpix**2).sum()
        chi2[i] = ((phot_Jy - fnu[i,:,wave_index])**2).sum()

    #idx = [chi2 == min(chi2)]
    #index = np.where(idx)[1]
    index = np.where(chi2 == np.nanmin(chi2))[0]
    #print('chi2 min =',np.nanmin(chi2))
    #print(index)
    plt.figure()
    plt.plot(ap , phot_Jy, 'r.') # MSX_A W/m^2-sr to Jy *7.133e12*scale**2/((180/np.pi)**2*3600**2)
    plt.plot(ap , fnu[int(index) ,:,wave_index], 'b.')
    plt.xlabel('AU')
    plt.ylabel('Jy')
    plt.title('%s  %s \n $\chi^2$ = %f %s \n MSX_%s'%(obj_name, models, chi2[index], names[index],flag))
    #plt.savefig(str('fit_'+ obj_name +'_' + models +'.png') , dpi=300)
    plt.show()
    return(chi2, wave_index)

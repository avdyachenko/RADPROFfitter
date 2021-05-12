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

def chi2_weights(phot_Jy, fnu, wave_index, models,names,flag, weights=np.ones_like(phot_Jy)):
    chi2 = np.zeros(fnu.shape[0])
    phot_Jy = phot_Jy/weights
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

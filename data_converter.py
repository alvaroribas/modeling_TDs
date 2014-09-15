######## Script to convert IRS spectra into pseudophotometric
######## datapoints for modeling the TDs

import asciitable
import numpy as np
import matplotlib.pyplot as plt
import pyfits
from scipy import interpolate

def remove_duplicates_func(seq):
    """ This function takes a list and returns
    the same without duplicate elements."""
    seen = set()
    seen_add = seen.add
    return [ x for x in seq if x not in seen and not seen_add(x)]


#####

avs_dictionary = {}
avs_dictionary['CS_Cha'] = 0.25
avs_dictionary['SZ_Cha'] = 1.90
avs_dictionary['T25'] = 0.78
avs_dictionary['T35'] = 3.5    
avs_dictionary['T56'] = 0.23
avs_dictionary['ISO_52'] = 1.3

obj = 'T25'
av = avs_dictionary[obj]

# path informations
path_main = '../objects/'
path_object = path_main + obj + '/' + obj +'_data/'

# Read the information of the filters
filters_info=asciitable.read(path_main+'filters_info.txt')
total_filters_names=filters_info['Filter']
total_filters_lmb=filters_info['Lmb'] #Armstrongs
total_filters_av_almb=filters_info['Av/Alambda']
total_filters_zp=filters_info['ZP']


# read the phot info and get it in right units
phot_file = asciitable.read(path_object + obj +'_phot.txt')
filters = phot_file['Filter']
lmb = phot_file['Lmb']
fluxes = phot_file['Value']
errors = phot_file['Error']
units = phot_file['Units']
detections = phot_file['Detection']
zp = phot_file['ZP']

filters_av_almb = list()
filters_zp = list()
for element in phot_file:
    index = np.where(total_filters_names == element['Filter'])[0][0]
    filters_av_almb.append(total_filters_av_almb[index])
    filters_zp.append(total_filters_zp[index])
filters_av_almb = np.array(filters_av_almb)
filters_zp = np.array(filters_zp)

# convert to jy
indexes_conversion = np.where(phot_file['Units'] == 'mag')[0]
error_fractions = errors[indexes_conversion] 
fluxes[indexes_conversion] = zp[indexes_conversion] * 10 ** (-0.4 * fluxes[indexes_conversion])
errors[indexes_conversion] = error_fractions * fluxes[indexes_conversion]

# derreden
fluxes = fluxes / 10**(-0.4 * av * filters_av_almb)
errors = errors / 10**(-0.4 * av * filters_av_almb)


# convert to flmb, fluxes in Jy, wavelength in microns
fluxes = (3e-8 * fluxes*1e3 / (lmb*1e4)**2)
errors = (3e-8 * errors*1e3 / (lmb*1e4)**2)

# Convert now to lmbflmb, fluxes in erg/cm2/s/A, wavelength in microns
lmb_flmb = fluxes * lmb * 1e4 #set lmb to A
lmb_flmb_e =  errors * lmb * 1e4 #set lmb to A

indexes_upper = np.where(detections == 0)[0]
lmb_flmb_e[indexes_upper] = 0.


file_to_write = open(path_object + obj +'_processed.txt','w')
file_to_write.write('#Filter\tLmb[um]\tLmb_flmb[erg/cm2/s]\tLmb_flmb_err[erg/cm2/s]\tDetection\n')
for index in xrange(len(fluxes)):
    file_to_write.write(filters[index]+'\t')
    file_to_write.write('{:.3e}\t'.format(lmb[index]))
    file_to_write.write('{:.3e}\t'.format(lmb_flmb[index]))
    file_to_write.write('{:.3e}\t'.format(lmb_flmb_e[index]))
    file_to_write.write('{:}\n'.format(detections[index]))
file_to_write.close()

file_lmb = open(path_object + obj + '.lambda','w')
lmb_unique = remove_duplicates_func(lmb)
lmb_unique = np.sort(lmb_unique)
for element in lmb_unique:
    file_lmb.write('{:.4e}\n'.format(element))
file_lmb.close()


##########################
## IRS Spectrum

# Derredening data

# Mathis1990 extinction law for spitzer (Rv=5)
mathis_lmb=[2.2,3.4,5.,7.,9.,9.7,10.,12.,15.,18.,20.,25.,35.]
mathis_alambda_aj=[0.382,0.182,0.095,0.07,0.157,0.2,0.192,0.098,0.053,0.083,0.075,0.0048,0.013]
mathis_interpol=interpolate.interp1d(mathis_lmb,mathis_alambda_aj,kind='linear')

#McClure2009 extinction law (lmb in microns)
mcclure=pyfits.open(path_main + 'McClure2009.fits')
mcclure_lmb=mcclure[1].data['lambda']
mcclure_alambda_ak1=mcclure[1].data['Al/AK1']
mcclure_alambda_ak2=mcclure[1].data['Al/AK2']
indexes=[(mcclure_lmb < 36) & (mcclure_lmb > 4)]
mcclure_lmb=mcclure_lmb[indexes]
mcclure_alambda_ak1=mcclure_alambda_ak1[indexes]
mcclure_alambda_ak2=mcclure_alambda_ak2[indexes]
mcclure_interpol1=interpolate.interp1d(mcclure_lmb,mcclure_alambda_ak1,kind='linear')
mcclure_interpol2=interpolate.interp1d(mcclure_lmb,mcclure_alambda_ak2,kind='linear')


irs=pyfits.open(path_object + obj +'_IRS.fits')
spectrum=irs[0].data
irs_lmb=spectrum[:,0] ## In microns
irs_fnu=spectrum[:,1]
irs_fnu_err = np.sqrt(spectrum[:,2]**2 + spectrum[:,3]**2 + spectrum[:,4]**2)
# get the errors in relative error
irs_fnu_rel_err = irs_fnu_err / irs_fnu


# cut the order1 between 7-14 and 20.5 - 35 microns 
order1=[(spectrum[:,8] == 1) & (((irs_lmb > 7.6) &(irs_lmb < 14.)) | ((irs_lmb > 20.5) & (irs_lmb < 35.)))]
# cut the order2 up to 20.5 microns
order2=[(spectrum[:,8] == 2) & (irs_lmb < 20.5)]
# get to corresponding values and sort them
lmb1=irs_lmb[order1]
lmb2=irs_lmb[order2]
irs_fnu1=irs_fnu[order1]
irs_fnu2=irs_fnu[order2]
irs_lmb=np.concatenate((lmb1,lmb2),axis=0)
irs_fnu=np.concatenate((irs_fnu1,irs_fnu2),axis=0)
irs_fnu=irs_fnu[np.argsort(irs_lmb)]
irs_lmb=irs_lmb[np.argsort(irs_lmb)]
#
aj=av*0.31
ak=av*0.13
print 'Aj:'+str(aj)
print 'Ak:'+str(ak)

# derreden
if (aj < 0.8):
    coeffs=mathis_interpol(irs_lmb)
    almbs=coeffs*aj
    irs_fnu_der=irs_fnu*10.**(-0.4*(-almbs))
if (aj > 0.8 and ak < 1):
    coeffs=mcclure_interpol1(irs_lmb)
    almbs=coeffs*ak
    irs_fnu_der=irs_fnu*10.**(-0.4*(-almbs))
else:
    coeffs=mcclure_interpol2(irs_lmb)
    almbs=coeffs*ak
    irs_fnu_der=irs_fnu*10.**(-0.4*(-almbs))
print 'Min:'+str(almbs.min())+'/Max:'+str(almbs.max())
# convert from Jy to erg/cm2/s/Hz
irs_fnu_der=irs_fnu_der*1e-23
# convert to erg/cm2/s/A
irs_flmb=irs_fnu_der*3.e8*1.e2/irs_lmb**2

# sort everything and get lmb_flmb
indexes_real = np.isfinite(irs_flmb)
irs_lmb = irs_lmb[indexes_real]
irs_flmb = irs_flmb[indexes_real]
irs_fnu_rel_err = irs_fnu_rel_err[indexes_real]
irs_lmbflmb = irs_lmb *1e4 * irs_flmb # irs_lmb in microns

# bring back the errors, only when the errors are real too
indexes_real = np.isfinite(irs_fnu_rel_err)
irs_lmb = irs_lmb[indexes_real]
irs_lmbflmb = irs_lmbflmb[indexes_real]
irs_fnu_rel_err = irs_fnu_rel_err[indexes_real]
irs_lmbflmb_err = irs_lmbflmb * irs_fnu_rel_err

# now bin everything
n_bins = 10
len_bins = np.int(np.floor(len(irs_lmb) / n_bins))

lmb_binned = list()
lmbflmb_binned = list()
lmbflmb_err_binned = list()

for n_bin in xrange(n_bins):
    # compute indexes for binning
    indexes_bin = np.arange(n_bin * len_bins, (n_bin+1) * len_bins)
    # in the last case, take the remaining datapoints in the last bin
    if n_bin == n_bins-1:
        indexes_bin = np.arange(n_bin * len_bins,len(irs_lmb))
    
    lmb_value = np.mean(irs_lmb[indexes_bin])
    lmbflmb_value = np.mean(irs_lmbflmb[indexes_bin])
    lmbflmb_err_value = np.std(irs_lmbflmb[indexes_bin]) / np.sqrt(len(indexes_bin))
    # append the results
    lmb_binned.append(lmb_value)
    lmbflmb_binned.append(lmbflmb_value)
    lmbflmb_err_binned.append(lmbflmb_err_value)

# finally, convert it to np arrays
lmb_binned = np.array(lmb_binned)
lmbflmb_binned = np.array(lmbflmb_binned)
lmbflmb_err_binned = np.array(lmbflmb_err_binned)

# append the results to the files
file_to_write = open(path_object + obj +'_processed.txt','a')
file_lmb = open(path_object + obj + '.lambda','a')

for index in xrange(len(lmb_binned)):
    file_lmb.write('{:.4e}\n'.format(lmb_binned[index]))
    file_to_write.write('IRS_binned\t')
    file_to_write.write('{:.3e}\t'.format(lmb_binned[index]))
    file_to_write.write('{:.3e}\t'.format(lmbflmb_binned[index]))
    file_to_write.write('{:.3e}\t'.format(lmbflmb_err_binned[index]))
    file_to_write.write('1\n')


file_to_write.close()
file_lmb.close()


# OPTIONAL: plot to check
plot_to_check = True
if plot_to_check:
    plt.errorbar(irs_lmb,irs_lmbflmb,yerr=irs_lmbflmb_err,fmt='o',mec=None, ms=1, mfc='blue')
    plt.errorbar(lmb_binned,lmbflmb_binned,yerr=lmbflmb_err_binned,fmt='o',mfc='red',mec=None,ms=8,color='red')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(4,40)
    plt.show()





        

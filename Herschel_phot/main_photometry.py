##############################################################
# Álvaro Ribas, aribas@cab.inta-csic.es
# 30/09/2014
# HIPE script: usual modules such as numpy are not available
#

"""Main script to produce photometry
 Note this is a quite specialized script and is not intended for 
 general use, it will perform photometry on Cha I maps. 
 If you are planning to use it, you should copy the functions 
 in "useful_tools.py" and "photometry_functions.py" at the beginning
 of the script.
"""

# paths
path_main = '/pcdisk/stark/aribas/Desktop/modeling_TDs/remaps_Cha/'
path_data = '/pcdisk/stark/aribas/Desktop/modeling_TDs/phot_Cha/'


# read the table with coordinates in an arcane way
file_coords = open(path_data+'coords.txt','r')
data = file_coords.readlines()
file_coords.close()


names = []
ra_sources = []
dec_sources = []
# remove first line of the file, which is the header
data = data[1:]
# extract the name and coordinates, and append them to the variables
for line in data:
    s_line = line.rsplit()
    names.append(s_line[0])
    ra_sources.append(s_line[1])
    dec_sources.append(s_line[2])

# total number of sources
n_objects = len(names)


# define some useful lists and dictionaries
instruments = ['PACS','SPIRE']                   
maps_pacs = ['scanamorphos','jscanam','unimap']
maps_spire = ['scanamorphos','destriper','unimap']
maps = {'PACS':maps_pacs, 'SPIRE':maps_spire}

                
lmbs_names = {'70':'BLUE',\
              '100':'GREEN',\
                '160_slow':'RED',\
                '160_fast':'RED',\
                '250':'PSW',\
                '350':'PMW',\
                '500':'PLW'}
psfs = {'70':6., '100':8., '160_slow':12., '160_fast':12.,'250':18., '350':25., '500':36.}
recommended_apertures = {'70':[12., 20., 30.],\
                            '100':[12., 20., 30.],\
                            '160_slow':[22., 30., 40.],\
                            '160_fast':[22., 30., 40.],\
                            '250':[22., 60., 90.],\
                            '350':[30., 60., 90.],\
                            '500':[42., 60., 90.]}
spire_beams = {'250':465.4, '350':822.5, '500':1768.6}


# get the calibration trees
obs_fast = getObservation(1342213179,instrument='PACS',useHsa=True)
obs_slow = getObservation(1342224783,instrument='PACS',useHsa=True)
cal_pacs_fast = getCalTree(obs = obs_fast)
cal_pacs_slow = getCalTree(obs = obs_slow)
cal_spire = spireCal(pool='spire_cal_12_3')
del(obs_fast,obs_slow)

# create the aperture correctors for SPIRE (interpolation within HIPE)
values = cal_spire.refs["Phot"].product.refs["RadialCorrBeam"].product["normArea"]
psw_corrector = LinearInterpolator(values['radius'].data,values['PSW'].data)
pmw_corrector = LinearInterpolator(values['radius'].data,values['PMW'].data)
plw_corrector = LinearInterpolator(values['radius'].data,values['PLW'].data)
spire_ap_correctors = {'250':psw_corrector,'350':pmw_corrector,'500':plw_corrector}

# and the color correctors. We will asume a blackbody slope (alpha=2, according to SPIRE)
spire_col_corrector = {'250':0.9454,'350':0.9481,'500':0.9432}



# We do not want to load all maps in memory at the same
# time. Therefore, we will save a photometry file per map, and then
# create a complete one using python scripts (where the world is kinder) 
for instrument in instruments:
    # Get the map names and wavelengths
    maps_instrument = maps[instrument]
    if instrument == 'PACS':
        lmbs = ['70','100','160_slow','160_fast']
    else:
        lmbs = ['250','350','500']

    # run throught the lmb
    for lmb in lmbs:
        band = lmbs_names[lmb]

        # compute the available apertures, from the psf value to the
        # recommended aperture one.
        aperture_values = range(psfs[lmb],recommended_apertures[lmb][0]+1)

        # Define the calibration tree to be used.
        # For SPIRE, also define the aperture correction
        for aperture in aperture_values:
            if (lmb == '100') or (lmb == '160_slow'):
                cal_to_use = cal_pacs_slow
            elif (lmb == '70') or (lmb == '160_fast'):
                cal_to_use = cal_pacs_fast
            else:
                ap_corrector = spire_ap_correctors[lmb]        # This is a function
                color_correction = spire_col_corrector[lmb]      # This is a value
                

        # for that particular wavelength, there are several maps
        for particular_map in maps_instrument:
            path_map = path_main+instrument+'/'+particular_map+'/ChaI_'+lmb+'.fits'                

            # if the map is a scanamorphos one, use importScanAmorphos task
            # otherwise load the image as a HIPE one
            if particular_map == 'scanamorphos':
                image = importScanAmorphos(file = path_map)
            else:
                image = fitsReader(file = path_map)

            # set image units correctly just in case
            if instrument == 'PACS':
                # for any weird reason, this makes HIPE slower. Use it only for unimap
                if particular_map == 'unimap':
                    image.setUnit('Jy/pixel')
                image = convertImageUnit(image,newUnit = 'Jy/pixel')
            else:
                image.setUnit('Jy/beam')
                image = convertImageUnit(image,newUnit='Jy/pixel',beamArea=spire_beams[lmb])

            # Open the file in which we will save the photometry, and save the header
            file_photometry = open(path_data+'phot_'+lmb+'_'+particular_map+'.txt','w')
            file_photometry.write('Name\tRA\tDec\t')
            for aperture in aperture_values:
                file_photometry.write('F_'+lmb+'_'+particular_map+'_'+str(int(aperture))+'\t')
            # also include one column for the stdev of the sky
            file_photometry.write('Sky_stdev_'+lmb+'_'+particular_map)
            file_photometry.write('\n')
           
            # Everything is set up, so we can start to extract the photometry
            for index_source in xrange(n_objects):
                name = names[index_source]
                ra_source = float(ra_sources[index_source])
                dec_source = float(dec_sources[index_source])
                # The coordinates have to be given in string to the photometry functions
                ra_source_string = convertRA(ra_source)
                dec_source_string = convertDec(dec_source)
                
                # write that in the file
                file_photometry.write(name+'\t'+str(ra_source)+'\t'+str(dec_source)+'\t')
                
                # Run through all the apertures to get the photometry
                for aperture_value in aperture_values:
                    radii_phot = [aperture_value,recommended_apertures[lmb][1],recommended_apertures[lmb][2]]
                    try:
                        if instrument == 'PACS':
                            value = pacsPhotometry(ra_source_string, dec_source_string,\
                                                    radii_phot, image, band, calTree = cal_to_use)
                        if instrument == 'SPIRE':
                            # color and aperture corrections have to be computed for SPIRE
                            ap_correction = ap_corrector(aperture_value)
                            value = spirePhotometry(ra_source_string, dec_source_string,\
                                                     radii_phot, image, 
                                                    ap_correction = ap_correction, color_correction = color_correction)
                        # discard photometry if the coordinates are outside the image
                        is_good = True
                        if value == 0.0:
                            is_good = False
                            value = 'NaN'
                    except:
                        # This means the object is outside the map, and it is trying to compute photometry on NaNs. That is why it crashed.
                        is_good = False
                        value = 'NaN'

                    file_photometry.write(str(value)+'\t')


                sky_std = 'NaN'
                # Background stddev calculation, only if the were no problems during
                if is_good:
                    # we will try generating 1000 positions around, and compute the stddev
                    sky_values = []
                    # aperture value is now the recommended one
                    r_min = 2 *aperture_value
                    r_max = 4 *aperture_value
                    try:
                        is_out = False
                        for index_aperture in xrange(1000):
                            position = coords_randomizer(ra_source, dec_source, r_min, r_max)
                            ra_new = convertRA(position[0])
                            dec_new = convertDec(position[1])
                            if instrument == 'PACS':
                                value = pacsPhotometry(ra_new, dec_new,\
                                                    radii_phot, image, band, calTree = cal_to_use, is_for_bg = True)
                            if instrument == 'SPIRE':
                                # ap correction already computed before
                                value = spirePhotometry(ra_new, dec_new,\
                                                     radii_phot, image, 
                                                    ap_correction = ap_correction, color_correction = color_correction, is_for_bg = True)
                            # If *ANY* value is out, we will leave
                            # that source. It probably means wrong
                            # background estimation from the annulus
                            if value == 0:
                                is_out = True
                            sky_values.append(value)
                        # and now check if everything is ok and assing the stdev
                        if not(is_out):
                            sky_std = STDDEV(sky_values)
                                
                    except:
                        # The source is too close to a border
                        sky_std = 'NaN'
                        
                # append the value
                file_photometry.write(str(sky_std)+'\n')
            # finally close the file
            file_photometry.close()

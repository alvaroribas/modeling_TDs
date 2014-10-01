##############################################################
# Álvaro Ribas, aribas@cab.inta-csic.es
# 30/09/2014
# HIPE functions: usual modules such as numpy are not available

def pacsPhotometry(ra, dec, radii, image_pacs, band_pacs, calTree):
    """ Function to compute annular aperture photometry on PACS data.
    It returns photometry already aperture-corrected photometry. 
    - coords_source: strings in hh:mm:ss.ss and (-) dd:mm:ss.s format
    - radii: aperture, and inner/outer radii for the background calculation.
             in arcsecs.
    - image_pacs: the image on which to compute the photometry
    - band_pacs: corresponding band, "Blue", "Green" or "Red"
    - calTree: calibration tree of the observations
    """
    # coordinates, must be given as strings in hh:mm:ss.ss 
    # and (-) dd:mm:ss.s format
    ra_photometry = ra
    dec_photometry = dec
    # aperture values
    raper = radii[0]
    rskyin = radii[1]
    rskyout = radii[2]
    # obtain photometry
    phot_value =  annularSkyAperturePhotometry(image=image_pacs,centroid=False,\
                                               fractional=True,algorithm=4,\
                                               centerRA=ra_photometry ,centerDec=dec_photometry,\
                                               radiusArcsec=raper,\
                                               innerArcsec=rskyin,outerArcsec=rskyout)
    # obtain aperture correction factor
    corrected_phot_value = photApertureCorrectionPointSource(apphot=phot_value, band=band_pacs,\
                                                             calTree = calTree, \
                                                             responsivityVersion=6)
    # We will NOT make any color correction for PACS, 
    # where the SED is flat.
    measured_flux = corrected_phot_value["Results table"]["Total flux"].data[2]
    #
    return measured_flux



def spirePhotometry(ra, dec, radii, image_spire, ap_correction = 1., color_correction = 1.):
    """ Function to compute annular aperture photometry on PACS data.
    It returns photometry already aperture-corrected photometry. 
    - coords_source: strings in hh:mm:ss.ss and (-) dd:mm:ss.s format
    - radii: aperture, and inner/outer radii for the background calculation.
             in arcsecs.
    - image_spire: the image on which to compute the photometry
    - ap_correction: aperture correction, to be computed outside this function. Set to 1 by default.
    - color_correction: color correction, to be computed outside this function. Set to 1 by default.
    """
    # coordinates, must be given as strings in hh:mm:ss.ss 
    # and (-) dd:mm:ss.s format
    ra_photometry = ra
    dec_photometry = dec
    # aperture values
    raper = radii[0]
    rskyin = radii[1]
    rskyout = radii[2]
    # obtain photometry
    phot_value =  annularSkyAperturePhotometry(image=image_spire,centroid=False,\
                                               fractional=True,algorithm=4,\
                                               centerRA=ra_photometry ,centerDec=dec_photometry,\
                                               radiusArcsec=raper,\
                                               innerArcsec=rskyin,outerArcsec=rskyout)

    value = phot_value["Results table"]["Total flux"].data[2]
    # apply aperture correction (divide)
    value = value / ap_correction
    # apply color correction (multiply)
    value = value * color_correction
    return value

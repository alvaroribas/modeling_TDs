def pacsPhotometry(coords_source, radii, image_pacs, band_pacs, calTree):
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
    ra_photometry = coords_source[0]
    dec_photometry = coords_source[1]
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
    corrected_phot_value = photApertureCorrectionPointSource(apphot=phot_value, band=pacs_band,\
                                                             calTree = calTree, \
                                                             responsivityVersion=6)
    # We will NOT make any color correction for PACS, 
    # where the SED is flat.
    measured_flux = corrected_phot_value["Results table"]["Total flux"].data[2]
    #
    return measured_flux



def spirePhotometry():(coords_source, radii, image_pacs, band_pacs, calTree):
    """ Function to compute annular aperture photometry on PACS data.
    It returns photometry already aperture-corrected photometry. 
    - coords_source: strings in hh:mm:ss.ss and (-) dd:mm:ss.s format
    - radii: aperture, and inner/outer radii for the background calculation.
             in arcsecs.
    - image_spire: the image on which to compute the photometry
    - band_spire: corresponding band, "PSW", "PMW" or "PLW"
    - band_spire_lmb: corresponding band, "250", "350", "500"
    - calTree: calibration tree of the observations
    """
    # coordinates, must be given as strings in hh:mm:ss.ss 
    # and (-) dd:mm:ss.s format
    ra_photometry = coords_source[0]
    dec_photometry = coords_source[1]
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
    # obtain aperture correction factor.
    # For SPIRE, it has to be obtained from the calTree



    corrected_phot_value = photApertureCorrectionPointSource(apphot=phot_value, band=pacs_band,\
                                                             calTree = calTree, \
                                                             responsivityVersion=6)
    # We will NOT make any color correction for PACS, 
    # where the SED is flat.
    measured_flux = corrected_phot_value["Results table"]["Total flux"].data[2]
    #
    return measured_flux
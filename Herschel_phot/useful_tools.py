##############################################################
# Álvaro Ribas, aribas@cab.inta-csic.es
# 30/09/2014
# HIPE script: usual modules such as numpy are not available


def convertRA(ra):
    """ Convert a RA value in format dd.dddd to dd:mm:ss.ss 
	Needed for the AnnularSkyAperturePhotometry tool"""
    hours = ra * 24. / 360
    minutes = (hours - int(hours)) * 60.
    seconds = (minutes - int(minutes)) * 60.
    
    hours = int(hours)
    minutes = int(minutes)
    aux_hours = ''
    aux_minutes = ''
    aux_seconds = ''
    if hours < 10:
        aux_hours = '0'
    if minutes < 10:
        aux_minutes = '0'
    if seconds < 10:
        aux_seconds = '0'

    result = aux_hours + str(hours) +':'+\
             aux_minutes + str(minutes) +':'+\
             aux_seconds + str(seconds)
    return result


def convertDec(dec):
    """ Convert a Dec value in format dd.dddd to +-dd:mm:ss.ss
		Needed for the AnnularSkyAperturePhotometry tool"""

    sign = 1
    aux_sign = ''
    if dec < 0:
        sign = -1
        aux_sign = '-' 

    hours = sign * dec
    minutes = (hours - int(hours)) * 60.
    seconds = (minutes - int(minutes)) * 60.
    
    hours = int(hours)
    minutes = int(minutes)
    aux_hours = ''
    aux_minutes = ''
    aux_seconds = ''
    if hours < 10:
        aux_hours = '0'
    if minutes < 10:
        aux_minutes = '0'
    if seconds < 10:
        aux_seconds = '0'

    result = aux_sign + aux_hours + str(hours) +':'+\
             aux_minutes + str(minutes) +':'+\
             aux_seconds + str(seconds)
    return result


def coords_randomizer(ra_center, dec_center, r_min, r_max):
    """Function to randomize a central coordinate with a central coordinate and
    given inner and outer radii).  It will do it with a uniform
    probability.  Note that this is NOT the correct formula, it will
    only work for small r_min y r_max.
    - r_min and r_max, in arcsec
    """
    # convert r_max r_min to deg
    r_min /= 3600.
    r_max /= 3600.
    random1 = random.random()
    random2 = random.random()
    radius_displacement = r_min + (r_max - r_min)*random1
    orientation = 2 * 3.1415 * random2

    new_ra = ra_center + radius_displacement * math.cos(orientation)
    new_dec = dec_center + radius_displacement * math.sin(orientation)
    new_coords = [new_ra, new_dec]
    return new_coords

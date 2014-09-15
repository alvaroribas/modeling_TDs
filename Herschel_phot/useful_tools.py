##############################################################
# Álvaro Ribas, aribas@cab.inta-csic.es
# 15/09/2014
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

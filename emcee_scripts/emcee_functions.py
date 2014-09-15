import os
import numpy as np
import matplotlib
import logging
import glob
import astropy, astropy.io.ascii
_log = logging.getLogger('mcfostpy')


# Define the probability function as likelihood * prior.
def lnprior(params):
    """ This function creates all priors for the MCMC phase 
    - prior for **TOTAL** mass (convert to dust mass and BEWARE of log scale!):         0.001Msun < Mdisk < 0.1 Mcorresponding star
    - prior for rin:                                                                    5 AU < Rin < 100 AU
    - prior for rot:                                                                    rin + 1 < Rout < 1000 AU
    - prior for surface density profile:                                                -3 < surf_exp < 1  
    - prior for flaring exponent:                                                       -3 < flaring_exp < -1  
    - prior for scale height (at 5AU):                                                  0.01 AU < scale_height < 5 AU
    - prior for dust max size (in um and log scale!):                                   -2 < a_max < 2 AU

    ALL OF THEM WILL BE FLAT PRIORS OR 0 VALUE
    """
    
    # create prior for disk mass
    dust_mass = 10.**params['disk_dust_mass']
    disk_total_mass = dust_mass*(1+params['gas_to_dust_ratio'])
    prior_mass = (disk_total_mass > -3) & (disk_total_mass < params['stellar_mass'])

    # create prior for rin
    prior_rin = (params['rin'] > 5.) & (params['rin'] < 100.)

    # create prior for rout
    prior_rout = (params['rout'] > params['rin'] + 1.) & (params['rout'] < 1000.)

    # create prior for surface density profile
    prior_surf_density = (params['surf_density_profile'] > -3.) & (params['surf_density_profile'] < 1.)

    # create prior for flaring angle exponent
    prior_flaring_angle = (params['flaring_angle_exp'] > -3.) & (params['flaring_angle_exp'] < -1.)

    # create prior for scale height
    prior_scale_height = (params['scale_height'] > 0.01) & (params['scale_height'] < 5.)    

    # create prior for dust max size
    prior_a_max = (np.log10(params['a_max']) > -2) & (np.log10(params['a_max']) < 2)    
    
  
    if prior_mass and prior_rin and prior_rout and prior_surf_density and prior_flaring_angle and prior_scale_height and prior_a_max:
        return 0.0
    return -np.inf



def lnlike(theta, x, y, yerr):
    m, b, lnf = theta
    model = m * x + b
    inv_sigma2 = 1.0/(yerr**2 + model**2*np.exp(2*lnf))
    return -0.5*(np.sum((y-model)**2*inv_sigma2 - np.log(inv_sigma2)))

def lnprob(theta, x, y, yerr):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x, y, yerr)

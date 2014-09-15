from my_create_params import *
import numpy



def generate_mcfost_file(path,params):
    mcfost_file = Paramfile(path)


    ## Global parameters: we may want to leave them constant for all the disks

    # Photons
    mcfost_file.nbr_photons_eq_th = params['nbr_photons_eq_th']
    mcfost_file.nbr_photons_lambda = params['nbr_photons_lambda']
    mcfost_file.nbr_photons_image = params['nbr_photons_image']
    # #Grid geometry and size
    mcfost_file.grid_geometry = params['grid_geometry']       # 1 = cylindrical, 2 = spherical, 3 = Voronoi tesselation (this is in beta, please ask Christophe)
    mcfost_file.grid_n_rad = params['grid_n_rad']             # n_rad (log distribution)
    mcfost_file.grid_nz = params['grid_nz']                   # nz (or n_theta)
    mcfost_file.grid_n_az = params['grid_n_az']               # n_az
    mcfost_file.grid_n_rad_in = params['grid_n_rad_in']       # n_rad_in
    #Maps
    mcfost_file.map_grid_nx = params['map_grid_nx']
    mcfost_file.map_grid_ny = params['map_grid_ny']
    mcfost_file.map_grid_size = params['map_grid_size'].      # AU
    mcfost_file.mc = params['mc']
    mcfost_file.mc_n_bin_incl = params['mc_n_bin_incl']
    mcfost_file.mc_n_bin_az = params['mc_n_bin_az']
    mcfost_file.rt_imin = params['rt_imin']
    mcfost_file.rt_imax = params['rt_imax']
    mcfost_file.rt_n_inlc = params['rt_n_inlc']
    mcfost_file.rt_is_centered = params['rt_is_centered']
    mcfost_file.distance = params['distance']                 # distance (pc)
    mcfost_file.pa = params['pa']                             # disk PA


    ## Particular parameters

    # Wavelength
    mcfost_file.n_lambda = params['n_lambda']                                         # micron
    mcfost_file.lambda_min = params['lambda_min']                                     # micron
    mcfost_file.lambda_max = params['lambda_max']                                     # micron
    mcfost_file.use_default_wavelength_grid = params['use_default_wavelength_grid']
    mcfost_file.wavelength_file = params['wavelength_file']

    # Zones
    mcfost_file.n_zones = params['n_zones']
    # Zone 1, always needed
    mcfost_file.zone1_type = params['zone1_type']                                         # zone type : 1 = disk, 2 = tappered-edge disk, 3 = envelope, 4 = debris disk, 5 = wall
    mcfost_file.zone1_dust_mass = params['zone1_dust_mass']                               # dust mass
    mcfost_file.zone1_gas_to_dust_mass_ratio = params['zone1_gas_to_dust_mass_ratio']      
    mcfost_file.zone1_scale_height = params['zone1_scale_height']                         # scale height
    mcfost_file.zone1_ref_radius = params['zone1_ref_radius']                             # reference radius (AU), unused for envelope
    mcfost_file.zone1_vert_profile_exp = params['zone1_vert_profile_exp']                 # vertical profile exponent (only for debris disk)
    mcfost_file.zone1_rin = params['zone1_rin']                                           # inner radius 
    mcfost_file.zone1_edge = params['zone1_edge']
    mcfost_file.zone1_rout = params['zone1_rout']                                         # inner radius
    mcfost_file.zone1_flaring_exp = params['zone1_flaring_exp']                           # flaring exponent, unused for envelope
    mcfost_file.zone1_surf_density_exp = params['zone1_surf_density_exp']                 # surface density exponent (or -gamma for tappered-edge disk), usually < 0
    mcfost_file.zone1_minusgamma_exp = params['zone1_minusgamma_exp']                     # -gamma_exp (or alpha_in & alpha_out for debris disk)
    # Zone 2 if needed:
    if mcfost_file.n_zones > 1:
        mcfost_file.zone2_type = params['zone2_type']                                         # zone type : 1 = disk, 2 = tappered-edge disk, 3 = envelope, 4 = debris disk, 5 = wall
        mcfost_file.zone2_dust_mass = params['zone2_dust_mass']                               # dust mass
        mcfost_file.zone2_gas_to_dust_mass_ratio = params['zone2_gas_to_dust_mass_ratio']      
        mcfost_file.zone2_scale_height = params['zone2_scale_height']                         # scale height
        mcfost_file.zone2_ref_radius = params['zone2_ref_radius']                             # reference radius (AU), unused for envelope
        mcfost_file.zone2_vert_profile_exp = params['zone2_vert_profile_exp']                 # vertical profile exponent (only for debris disk)
        mcfost_file.zone2_rin = params['zone2_rin']                                           # inner radius 
        mcfost_file.zone2_edge = params['zone2_edge']
        mcfost_file.zone2_rout = params['zone2_rout']                                         # inner radius
        mcfost_file.zone2_flaring_exp = params['zone2_flaring_exp']                           # flaring exponent, unused for envelope
        mcfost_file.zone2_surf_density_exp = params['zone2_surf_density_exp']                 # surface density exponent (or -gamma for tappered-edge disk), usually < 0
        mcfost_file.zone2_minusgamma_exp = params['zone2_minusgamma_exp']                     # -gamma_exp (or alpha_in & alpha_out for debris disk)
    # Zone 3 if needed:
    if mcfost_file.n_zones > 2:
        mcfost_file.zone3_type = params['zone3_type']                                         # zone type : 1 = disk, 2 = tappered-edge disk, 3 = envelope, 4 = debris disk, 5 = wall
        mcfost_file.zone3_dust_mass = params['zone3_dust_mass']                               # dust mass
        mcfost_file.zone3_gas_to_dust_mass_ratio = params['zone3_gas_to_dust_mass_ratio']      
        mcfost_file.zone3_scale_height = params['zone3_scale_height']                         # scale height
        mcfost_file.zone3_ref_radius = params['zone3_ref_radius']                             # reference radius (AU), unused for envelope
        mcfost_file.zone3_vert_profile_exp = params['zone3_vert_profile_exp']                 # vertical profile exponent (only for debris disk)
        mcfost_file.zone3_rin = params['zone3_rin']                                           # inner radius 
        mcfost_file.zone3_edge = params['zone3_edge']
        mcfost_file.zone3_rout = params['zone3_rout']                                         # inner radius
        mcfost_file.zone3_flaring_exp = params['zone3_flaring_exp']                           # flaring exponent, unused for envelope
        mcfost_file.zone3_surf_density_exp = params['zone3_surf_density_exp']                 # surface density exponent (or -gamma for tappered-edge disk), usually < 0
        mcfost_file.zone3_minusgamma_exp = params['zone3_minusgamma_exp']                     # -gamma_exp (or alpha_in & alpha_out for debris disk)

    # Grains
    # Grain 1, always needed
    mcfost_file.grain1_n_species = params['grain1_n_species']    # Number of species
    mcfost_file.grain1_amin = params['grain1_amin']              # amin (um)
    mcfost_file.grain1_amax = params['grain1_amax']              # amax (um)
    mcfost_file.grain1_aexp = params['grain1_aexp']              # aexp
    mcfost_file.grain1_n_grains = params['grain1_n_grains']      # n_grains (log distrib)
    # Grain 2 if needed:
    if mcfost_file.n_zones > 1:
        mcfost_file.grain2_n_species = params['grain2_n_species']    # Number of species
        mcfost_file.grain2_amin = params['grain2_amin']              # amin (um)
        mcfost_file.grain2_amax = params['grain2_amax']              # amax (um)
        mcfost_file.grain2_aexp = params['grain2_aexp']              # aexp
        mcfost_file.grain2_n_grains = params['grain2_n_grains']      # n_grains (log distrib)
    # Grain 3 if needed:
    if mcfost_file.n_zones > 2:
        mcfost_file.grain3_n_species = params['grain3_n_species']    # Number of species
        mcfost_file.grain3_amin = params['grain3_amin']              # amin (um)
        mcfost_file.grain3_amax = params['grain3_amax']              # amax (um)
        mcfost_file.grain3_aexp = params['grain3_aexp']              # aexp
        mcfost_file.grain3_n_grains = params['grain3_n_grains']      # n_grains (log distrib)

    # Star
    mcfost_file.n_stars = params['n_stars']                          # Number of stars
    mcfost_file.star1_teff = params['star1_teff'] 
    mcfost_file.star1_radius = params['star1_radius'] 
    mcfost_file.star1_mass = params['star1_mass'] 
    mcfost_file.star1_x = params['star1_x']                          # AU
    mcfost_file.star1_y = params['star1_y']                          # AU
    mcfost_file.star1_z = params['star1_z']                          # AU
    mcfost_file.star1_blackbody = params['star1_blackbody'] 
    mcfost_file.star1_spectrum = params['star1_spectrum'] 
    mcfost_file.star1_fUV = params['star1_fUV']                                    
    mcfost_file.star1_slopefUV = params['star1_slopefUV'] 

    mcfost_file.writeFile()

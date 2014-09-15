import os
import numpy as np
import matplotlib
import logging
import glob
import astropy, astropy.io.ascii





class Paramfile:
    """ Object class interface to MCFOST parameter files 

    It does help to change and write MCFOST parameter files.

    """

    _minimum_version = 2.15  # minimum MCFOST version for this code to run 
    def __init__(self,  path):
        """ Set initial parameter values.
        
        THESE ARE JUST DEFINITIONS, NO NEED TO CHANGE THIS PART UNLESS A VERY
        SPECIFIC PARAMETER MUST BE CHANGED
        """

        self.path = path       # path to save the file to
        self.version = '2.19'  # MCFOST version.

        
        # Number of photon packages
        self.nbr_photons_eq_th = np.nan   # T computation
        self.nbr_photons_lambda = np.nan  # SED computation
        self.nbr_photons_image = np.nan   # image computation

        
        # Wavelength
        self.n_lambda = np.nan                    # micron
        self.lambda_min = np.nan                  # micron
        self.lambda_max = np.nan                  # micron
        self.use_default_wavelength_grid = 'T'
        self.compute_teff = 'T'
        self.compute_sed = 'T'
        self.use_default_wavelength_grid = 'T'      # compute temperature, compute sed, use default wavelength grid
        self.wavelength_file = 'chose_it.lambda'        # wavelength file (if previous parameter is F)
        self.separate_contributions = 'F'  # Separate contributions
        self.compute_stokes = 'F'

        
        # Grid geometry and size
        self.grid_geometry = np.nan                              # 1 = cylindrical, 2 = spherical, 3 = Voronoi tesselation (this is in beta, please ask Christophe)
        self.grid_n_rad = np.nan   # n_rad (log distribution)
        self.grid_nz = np.nan      # nz (or n_theta)
        self.grid_n_az = np.nan         # n_az
        self.grid_n_rad_in = np.nan     # n_rad_in

        
        #Maps
        self.map_grid_nx = np.nan
        self.map_grid_ny = np.nan
        self.map_grid_size = np.nan  # AU
        self.mc = np.nan             # N bins inclination for the MC phase
        self.mc_n_bin_incl = 1
        self.mc_n_bin_az = 1         # This parameter is not used by the user
        self.rt_imin = 45.
        self.rt_imax = 45.
        self.rt_n_incl = 1.
        self.rt_is_centered = 'F'
        self.distance = np.nan                               # distance (pc)
        self.pa = 0.                                         # disk PA


        # Scattering method
        self.scattering_method = 0  # 0=auto, 1=grain prop, 2=cell prop
        self.scattering_theory = 1  # 1=Mie, 2=hg (2 implies the loss of polarizarion)

        
        # Symetries
        self.image_symmetry = 'T'	  # image symmetry
        self.central_symmetry = 'T'   # central symmetry
        self.axial_symmetry = 'T'     # axial symmetry (important only if N_phi > 1)

        
        # Disk physics
        self.dust_settling = 0
        self.dust_exp_strat = 0.50
        self.dust_a_strat = 1.
        self.dust_radial_migration = 'F'               # dust radial migration
        self.dust_sublimate_dust = 'F'                 # sublimate dust
        self.dust_hydrostatic_eq = 'F'                 # hydostatic equilibrium
        self.dust_viscous_heating = 'F'     		   # viscous heating
        self.dust_alpha_viscosity = 0.2

        
        # Number of zones : 1 zone = 1 density structure + corresponding grain properties
        self.n_zones = 1 # number of zones

            
        # ZONES
        #### We will define here 3 zones just in case, and use them only if n_zones is > 1
        # Zone1
        self.zone1_type = 1                           # zone type : 1 = disk, 2 = tappered-edge disk, 3 = envelope, 4 = debris disk, 5 = wall
        self.zone1_dust_mass = np.nan                 # dust mass
        self.zone1_gas_to_dust_mass_ratio = 100.      
        self.zone1_scale_height = np.nan              # scale height
        self.zone1_ref_radius = np.nan                # reference radius (AU), unused for envelope
        self.zone1_vert_profile_exp = 2               # vertical profile exponent (only for debris disk)
        self.zone1_rin = np.nan
        self.zone1_edge = 0.0
        self.zone1_rout = np.nan
        self.zone1_rc = 80.                           # Rc (AU) Rc is only used for tappered-edge disks (Rout set to 8*Rc if Rout==0)
        self.zone1_flaring_exp = np.nan               # flaring exponent, unused for envelope
        self.zone1_surf_density_exp = np.nan          # surface density exponent (or -gamma for tappered-edge disk), usually < 0
        self.zone1_minusgamma_exp = 0.0               # -gamma_exp (or alpha_in & alpha_out for debris disk)
        # Zone2
        self.zone2_type = 1                           # zone type : 1 = disk, 2 = tappered-edge disk, 3 = envelope, 4 = debris disk, 5 = wall
        self.zone2_dust_mass = np.nan                 # dust mass
        self.zone2_gas_to_dust_mass_ratio = 100.      
        self.zone2_scale_height = np.nan              # scale height
        self.zone2_ref_radius = np.nan                # reference radius (AU), unused for envelope
        self.zone2_vert_profile_exp = 2               # vertical profile exponent (only for debris disk)
        self.zone2_rin = np.nan
        self.zone2_edge = 0.0
        self.zone2_rout = np.nan
        self.zone2_rc = 80.                           # Rc (AU) Rc is only used for tappered-edge disks (Rout set to 8*Rc if Rout==0)
        self.zone2_flaring_exp = np.nan               # flaring exponent, unused for envelope
        self.zone2_surf_density_exp = np.nan          # surface density exponent (or -gamma for tappered-edge disk), usually < 0
        self.zone2_minusgamma_exp = 0.0               # -gamma_exp (or alpha_in & alpha_out for debris disk)
        # Zone3
        self.zone3_type = 1                           # zone type : 1 = disk, 2 = tappered-edge disk, 3 = envelope, 4 = debris disk, 5 = wall
        self.zone3_dust_mass = np.nan                 # dust mass
        self.zone3_gas_to_dust_mass_ratio = 100.      
        self.zone3_scale_height = np.nan              # scale height
        self.zone3_ref_radius = np.nan                # reference radius (AU), unused for envelope
        self.zone3_vert_profile_exp = 2               # vertical profile exponent (only for debris disk)
        self.zone3_rin = np.nan
        self.zone3_edge = 0.0
        self.zone3_rout = np.nan
        self.zone3_rc = 80.                           # Rc (AU) Rc is only used for tappered-edge disks (Rout set to 8*Rc if Rout==0)
        self.zone3_flaring_exp = np.nan               # flaring exponent, unused for envelope
        self.zone3_surf_density_exp = np.nan          # surface density exponent (or -gamma for tappered-edge disk), usually < 0
        self.zone3_minusgamma_exp = 0.0               # -gamma_exp (or alpha_in & alpha_out for debris disk)


        # Cavity : everything is empty above the surface
        self.is_cavity = 'F'	  	     # cavity 
        self.cavity_heigh = 15.          # height
        self.cavity_ref_radius = 50.	 # reference radius (AU)
        self.cavity_flaring = 1.5 	     # flaring exponent


        # GRAINS
        ### Grain properties suffer from the same "problem" that disk zones, so will put three.
        # Grain 1
        self.grain1_n_species = 1                         # Number of species
        self.grain1_type = 'Mie'                          # Grain type (Mie or DHS)
        self.grain1_n_components = 1
        self.grain1_mixing_rule = 2                       # Mixing rule (1 = EMT or 2 = coating)
        self.grain1_porosity = 0.0
        self.grain1_mass_fract = 1.0
        self.grain1_vmax = 0.9                            # Vmax (for DHS)
        self.grain1_dust_file = 'Draine_Si_sUV.dat'       # Optical indices file
        self.grain1_dust_volume = 1.0                     # Volume fraction                     
        self.grain1_heating_method = 1                    # Heating method : 1 = RE + LTE, 2 = RE + NLTE, 3 = NRE
        self.grain1_amin = 0.03                           # amin (um)
        self.grain1_amax = 1000.0                         # amax (um)
        self.grain1_aexp = 3.5                            # aexp
        self.grain1_n_grains = 50                         # n_grains (log distrib)
        # Zone 2
        self.grain2_n_species = 1                         # Number of species
        self.grain2_type = 'Mie'                          # Grain type (Mie or DHS)
        self.grain2_n_components = 1
        self.grain2_mixing_rule = 2                       # Mixing rule (1 = EMT or 2 = coating)
        self.grain2_porosity = 0.0
        self.grain2_mass_fract = 1.0
        self.grain2_vmax = 0.9                            # Vmax (for DHS)
        self.grain2_dust_file = 'Draine_Si_sUV.dat'       # Optical indices file
        self.grain2_dust_volume = 1.0                     # Volume fraction                     
        self.grain2_heating_method = 1                    # Heating method : 1 = RE + LTE, 2 = RE + NLTE, 3 = NRE
        self.grain2_amin = 0.03                           # amin (um)
        self.grain2_amax = 1000.0                         # amax (um)
        self.grain2_aexp = 3.5                            # aexp
        self.grain2_n_grains = 50                         # n_grains (log distrib)
        # Zone 3
        self.grain3_n_species = 1                         # Number of species
        self.grain3_type = 'Mie'                          # Grain type (Mie or DHS)
        self.grain3_n_components = 1
        self.grain3_mixing_rule = 2                       # Mixing rule (1 = EMT or 2 = coating)
        self.grain3_porosity = 0.0
        self.grain3_mass_fract = 1.0
        self.grain3_vmax = 0.9                            # Vmax (for DHS)
        self.grain3_dust_file = 'Draine_Si_sUV.dat'       # Optical indices file
        self.grain3_dust_volume = 1.0                     # Volume fraction                     
        self.grain3_heating_method = 1                    # Heating method : 1 = RE + LTE, 2 = RE + NLTE, 3 = NRE
        self.grain3_amin = 0.03                           # amin (um)
        self.grain3_amax = 1000.0                         # amax (um)
        self.grain3_aexp = 3.5                            # aexp
        self.grain3_n_grains = 50                         # n_grains (log distrib)
        
        ### I WILL NOT DO ANY MOLECULAR STUFF. INSTEAD, I WILL COPY THE PARAGRAPHS EXACTLY AS IN THE ORIGINAL FILE. NO USE ANYWAY
        
        #Star properties
        self.n_stars = 1                                      # Number of stars
        # star 1
        self.star1_teff = np.nan
        self.star1_radius = np.nan
        self.star1_mass = np.nan
        self.star1_x = 0.                                     # AU
        self.star1_y = 0.                                     # AU
        self.star1_z = 0.                                     # AU
        self.star1_blackbody = 'F'
        self.star1_spectrum = np.nan
        self.star1_fUV = 0.                                   
        self.star1_slopefUV = 2.2
        # star 2
        self.star2_teff = np.nan
        self.star2_radius = np.nan
        self.star2_mass = np.nan
        self.star2_x = 0.                                     # AU
        self.star2_y = 0.                                     # AU
        self.star2_z = 0.                                     # AU
        self.star2_blackbody = 'F'
        self.star2_spectrum = np.nan
        self.star2_fUV = 0.                                   
        self.star2_slopefUV = 2.2
        # star 3
        self.star3_teff = np.nan
        self.star3_radius = np.nan
        self.star3_mass = np.nan
        self.star3_x = 0.                                     # AU
        self.star3_y = 0.                                     # AU
        self.star3_z = 0.                                     # AU
        self.star3_blackbody = 'F'
        self.star3_spectrum = np.nan
        self.star3_fUV = 0.                                   
        self.star3_slopefUV = 2.2



        
    def writeFile(self):
        """ Write an MCFOST file.

        It assumes I am clever enough to get all the variable definitions right... """

        mcfost_file = open(self.path,'w')

        # version
        mcfost_file.write(str(self.version) +'          mcfost version \n')
        mcfost_file.write('\n')

        # Number of photon packages
        mcfost_file.write('#Number of photon packages\n')
        mcfost_file.write('  {:.3e}'.format(float(self.nbr_photons_eq_th))+'           nbr_photons_eq_th  : T computation\n')
        mcfost_file.write('  {:.3e}'.format(float(self.nbr_photons_lambda))+'           nbr_photons_lambda : SED computation\n')
        mcfost_file.write('  {:.3e}'.format(float(self.nbr_photons_image))+'           nbr_photons_image  : images computation\n')
        mcfost_file.write('\n')

        # Wavelength
        mcfost_file.write('#Wavelength\n')
        values = '  {:} '.format(int(self.n_lambda)) + '{:2e} '.format(float(self.lambda_min)) + '{:.2e} '.format(float(self.lambda_max))
        mcfost_file.write(values + '           n_lambda, lambda_min, lambda_max [mum]\n')
        values = '  ' + self.use_default_wavelength_grid + ' ' + self.compute_teff + ' ' + self.compute_sed
        mcfost_file.write(values + '           compute temperature?, compute sed?, use default wavelength grid ?\n')
        mcfost_file.write('  ' + self.wavelength_file + '           wavelength file (if previous parameter is F)\n')
        values = '  ' + self.separate_contributions + ' ' +self.compute_stokes
        mcfost_file.write('  ' + values + '           separation of different contributions?, stokes parameters?\n') 
        mcfost_file.write('\n')

        # Grid geometry and size
        mcfost_file.write('#Grid geometry and size\n')
        mcfost_file.write('  {:} '.format(int(self.grid_geometry)) + \
                          '  1 = cylindrical, 2 = spherical, 3 = Voronoi tesselation (this is in beta, please ask Christophe)\n')
        values = '  {:} '.format(int(self.grid_n_rad)) + '{:} '.format(int(self.grid_nz)) +\
            '{:} '.format(int(self.grid_n_az)) + '{:} '.format(int(self.grid_n_rad_in))
        mcfost_file.write(values + '    n_rad (log distribution), nz (or n_theta), n_az, n_rad_in\n')
        mcfost_file.write('\n')

        # Maps
        mcfost_file.write('#Maps\n')
        values = '  {:} '.format(int(self.map_grid_nx)) + '{:} '.format(int(self.map_grid_nx)) + '{:.3} '.format(float(self.map_grid_size))
        mcfost_file.write(values + '           grid (nx,ny), size [AU]\n')
        values = '  {:} '.format(int(self.mc)) + '{:} '.format(int(self.mc_n_bin_incl)) + '{:} '.format(int(self.mc_n_bin_az))
        mcfost_file.write(values + '              MC : N_bin_incl, N_bin_az\n')
        values = '  {:.3} '.format(float(self.rt_imin)) + '{:.3} '.format(float(self.rt_imax)) +\
            '{:} '.format(int(self.rt_n_incl)) + ' ' + self.rt_is_centered
        mcfost_file.write(values + '          RT: imin, imax, n_incl, centered ?\n')
        mcfost_file.write('  {:.3} '.format(float(self.distance)) + '  		  distance (pc)\n')
        mcfost_file.write('  {:.3} '.format(float(self.pa)) + ' 			  disk PA\n')
        mcfost_file.write('\n')

        # Scattering method
        mcfost_file.write('#Scattering method\n')
        mcfost_file.write('  {:} '.format(int(self.scattering_method)) + '	                  0=auto, 1=grain prop, 2=cell prop\n')
        mcfost_file.write('  {:} '.format(int(self.scattering_theory)) + '	                  1=Mie, 2=hg (2 implies the loss of polarizarion)\n')
        mcfost_file.write('\n')

        # Symmetries
        mcfost_file.write('#Symmetries\n')
        mcfost_file.write('  ' + self.image_symmetry + ' 	                  image symmetry\n')
        mcfost_file.write('  ' + self.central_symmetry + '	                  central symmetry\n')
        mcfost_file.write('  ' + self.axial_symmetry +  '	                  axial symmetry (important only if N_phi > 1)\n')
        mcfost_file.write('\n')

        # Disk physics
        mcfost_file.write('#Disk physics\n')
        values = '  {:} '.format(int(self.dust_settling)) + '{:.3} '.format(float(self.dust_exp_strat)) + '{:.3} '.format(float(self.dust_a_strat))
        mcfost_file.write(values + '	  dust_settling (0=no settling, 1=parametric, 2=Dubrulle, 3=Fromang), exp_strat, a_strat (for parametric settling)\n')
        mcfost_file.write('  ' + self.dust_radial_migration + '                       dust radial migration\n')
        mcfost_file.write('  ' + self.dust_sublimate_dust + '                         sublimate\n')
        mcfost_file.write('  ' + self.dust_hydrostatic_eq + '                        hydostatic equilibrium\n')
        mcfost_file.write('  ' + self.dust_viscous_heating + ' '+'{:1e}'.format(float(self.dust_alpha_viscosity)) + '		  viscous heating, alpha_viscosity\n')
        mcfost_file.write('\n')

        # Number of zones
        mcfost_file.write('#Number of zones : 1 zone = 1 density structure + corresponding grain properties\n')
        mcfost_file.write('  {:} '.format(int(self.n_zones))+'\n')
        mcfost_file.write('\n')

        # Density structure
        mcfost_file.write('#Density structure\n')

        ## Zone 1, which exisits for sure
        mcfost_file.write(' {:} '.format(int(self.zone1_type)) + '                       zone type : 1 = disk, 2 = tappered-edge disk, 3 = envelope, 4 = debris disk, 5 = wall\n')
        values = '  {:.3e} '.format(float(self.zone1_dust_mass)) + '{:.3} '.format(float(self.zone1_gas_to_dust_mass_ratio))
        mcfost_file.write(values + '	dust mass,  gas-to-dust mass ratio\n')
        values = '  {:.3e} '.format(float(self.zone1_scale_height)) + '{:.3} '.format(float(self.zone1_ref_radius)) +\
           '{:.3} '.format(float(self.zone1_vert_profile_exp))
        mcfost_file.write(values + '           scale height, reference radius (AU), unused for envelope, vertical profile exponent (only for debris disk)\n')
        values = '  {:.3e} '.format(float(self.zone1_rin)) + '{:.3} '.format(float(self.zone1_edge)) +\
           '{:.3} '.format(float(self.zone1_rout))+ '{:.3} '.format(float(self.zone1_rc))
        mcfost_file.write(values + '  Rin, edge, Rout, Rc (AU) Rc is only used for tappered-edge & debris disks (Rout set to 8*Rc if Rout==0)\n')
        mcfost_file.write('  {:.3} '.format(float(self.zone1_flaring_exp)) + '                      flaring exponent, unused for envelope\n')
        values = '  {:.3} '.format(float(self.zone1_surf_density_exp)) + '{:.3} '.format(float(self.zone1_minusgamma_exp))
        mcfost_file.write(values + '   surface density exponent (or -gamma for tappered-edge disk or volume density for envelope),'+\
                                   ' usually < 0, -gamma_exp (or alpha_in & alpha_out for debris disk)\n')
        ## Zone 2 if needed
        if self.n_zones > 1:
            mcfost_file.write(' {:} '.format(int(self.zone2_type)) +  '                       zone type : 1 = disk, 2 = tappered-edge disk,'+\
                                                                               ' 3 = envelope, 4 = debris disk, 5 = wall\n')
            values = '  {:.3e} '.format(float(self.zone2_dust_mass)) + '{:.3} '.format(float(self.zone2_gas_to_dust_mass_ratio))
            mcfost_file.write(values + '	dust mass,  gas-to-dust mass ratio\n')
            values = '  {:.3e} '.format(float(self.zone2_scale_height)) + '{:.3} '.format(float(self.zone2_ref_radius)) +\
            '{:.3} '.format(float(self.zone2_vert_profile_exp))
            mcfost_file.write(values + '           scale height, reference radius (AU), unused for envelope, vertical profile exponent (only for debris disk)\n')
            values = '  {:.3e} '.format(float(self.zone2_rin)) + '{:.3} '.format(float(self.zone2_edge)) +\
                '{:.3} '.format(float(self.zone2_rout))+ '{:.3} '.format(float(self.zone2_rc))
            mcfost_file.write(values + '  Rin, edge, Rout, Rc (AU) Rc is only used for tappered-edge & debris disks (Rout set to 8*Rc if Rout==0)\n')
            mcfost_file.write('  {:.3} '.format(float(self.zone2_flaring_exp)) + '                      flaring exponent, unused for envelope\n')
            values = '  {:.3} '.format(float(self.zone2_surf_density_exp)) + '{:.3} '.format(float(self.zone2_minusgamma_exp))
            mcfost_file.write(values + '   surface density exponent (or -gamma for tappered-edge disk or volume density for envelope),'+\
                                       ' usually < 0, -gamma_exp (or alpha_in & alpha_out for debris disk)\n')
        ## Zone 3 if needed
        if self.n_zones > 2:
            mcfost_file.write('  {:} '.format(int(self.zone3_type))  + '                       zone type : 1 = disk, 2 = tappered-edge disk,'+\
                                                                               ' 3 = envelope, 4 = debris disk, 5 = wall\n')
            values = '  {:.3e} '.format(float(self.zone_3dust_mass)) + '{:.3} '.format(float(self.zone_3gas_to_dust_mass_ratio))
            mcfost_file.write(values + '	dust mass,  gas-to-dust mass ratio\n')
            values = '  {:.3e} '.format(float(self.zone_3scale_height)) + '{:.3} '.format(float(self.zone_3ref_radius)) +\
            '{:.3} '.format(float(self.zone_3vert_profile_exp))
            mcfost_file.write(values + '           scale height, reference radius (AU), unused for envelope, vertical profile exponent (only for debris disk)\n')
            values = '  {:.3e} '.format(float(self.zone_3rin)) + '{:.3} '.format(float(self.zone_3edge)) +\
                '{:.3} '.format(float(self.zone_3rout))+ '{:.3} '.format(float(self.zone_3rc))
            mcfost_file.write(values + '  Rin, edge, Rout, Rc (AU) Rc is only used for tappered-edge & debris disks (Rout set to 8*Rc if Rout==0)\n')
            mcfost_file.write('  {:.3} '.format(float(self.zone_3flaring_exp)) + '                      flaring exponent, unused for envelope\n')
            values = '  {:.3} '.format(float(self.zone_3surf_density_exp)) + '{:.3} '.format(float(self.zone_3minusgamma_exp))
            mcfost_file.write(values + '   surface density exponent (or -gamma for tappered-edge disk or volume density for envelope),'+\
                                       ' usually < 0, -gamma_exp (or alpha_in & alpha_out for debris disk)\n')            
        mcfost_file.write('\n')

        # Cavity
        mcfost_file.write('#Cavity : everything is empty above the surface\n')
        mcfost_file.write('  ' + self.is_cavity + '	  	     	  cavity ?\n')
        values = '  {:.3} '.format(float(self.cavity_heigh)) + '{:.3} '.format(float(self.cavity_ref_radius))
        mcfost_file.write(values + '		  height, reference radius (AU)\n')
        mcfost_file.write('  {:.3} '.format(float(self.cavity_flaring))+ ' 			  flaring exponent\n')
        mcfost_file.write('\n')

        # Grains
        mcfost_file.write('#Grain properties\n')
        ## Grain 1, which exist for sure
        mcfost_file.write('  {:} '.format(int(self.grain1_n_species))+ '   Number of species\n')
        values = '  ' + self.grain1_type + ' {:} '.format(int(self.grain1_n_components)) + '{:} '.format(int(self.grain1_mixing_rule)) +\
            '{:.3} '.format(float(self.grain1_porosity)) + '{:.3} '.format(float(self.grain1_mass_fract))+ '{:.3} '.format(float(self.grain1_vmax))
        mcfost_file.write(values + ' Grain type (Mie or DHS), N_components, mixing rule (1 = EMT or 2 = coating),'+\
                                   '  porosity, mass fraction, Vmax (for DHS)\n')
        values = '  ' + self.grain1_dust_file + ' {:.3} '.format(float(self.grain1_dust_volume))
        mcfost_file.write(values + '  Optical indices file, volume fraction\n')
        mcfost_file.write('  {:} '.format(int(self.grain1_heating_method)) + '	                  Heating method : 1 = RE + LTE, 2 = RE + NLTE, 3 = NRE\n')
        values = '  {:.3} '.format(float(self.grain1_amin )) + '{:.3} '.format(float(self.grain1_amax)) +\
          '{:.3} '.format(float(self.grain1_aexp)) + '{:} '.format(int(self.grain1_n_grains)) 
        mcfost_file.write(values + ' 	  amin, amax [mum], aexp, n_grains (log distribution)\n')
        ## Grain 2 if needed
        if self.n_zones > 1:
            mcfost_file.write('  {:} '.format(int(self.grain2_n_species))+ '   Number of species\n')
            values = '  ' + self.grain2_type + ' {:} '.format(int(self.grain2_n_components)) + '{:} '.format(int(self.grain2_mixing_rule)) +\
                '{:.3} '.format(float(self.grain2_porosity)) + '{:.3} '.format(float(self.grain2_mass_fract))+ '{:.3} '.format(float(self.grain2_vmax))
            mcfost_file.write(values + ' Grain type (Mie or DHS), N_components, mixing rule (1 = EMT or 2 = coating),'+\
                                       '  porosity, mass fraction, Vmax (for DHS)\n')
            values = '  ' + self.grain2_dust_file + ' {:.3} '.format(float(self.grain2_dust_volume))
            mcfost_file.write(values + '  Optical indices file, volume fraction\n')
            mcfost_file.write('  {:} '.format(int(self.grain2_heating_method)) + '	                  Heating method : 1 = RE + LTE, 2 = RE + NLTE, 3 = NRE\n')
            values = '  {:.3} '.format(float(self.grain2_amin )) + '{:.3} '.format(float(self.grain2_amax)) +\
                '{:.3} '.format(float(self.grain2_aexp)) + '{:} '.format(int(self.grain2_n_grains)) 
            mcfost_file.write(values + ' 	  amin, amax [mum], aexp, n_grains (log distribution)\n')
        ## Grain 3 if needed
        if self.n_zones > 1:
            mcfost_file.write('  {:} '.format(int(self.grain3_n_species))+ '   Number of species\n')
            values = '  ' + self.grain3_type + ' {:} '.format(int(self.grain3_n_components)) + '{:} '.format(int(self.grain3_mixing_rule)) +\
                '{:.3} '.format(float(self.grain3_porosity)) + '{:.3} '.format(float(self.grain3_mass_fract))+ '{:.3} '.format(float(self.grain3_vmax))
            mcfost_file.write(values + ' Grain type (Mie or DHS), N_components, mixing rule (1 = EMT or 2 = coating),'+\
                                       '  porosity, mass fraction, Vmax (for DHS)\n')
            values = '  ' + self.grain3_dust_file + ' {:.3} '.format(float(self.grain3_dust_volume))
            mcfost_file.write(values + '  Optical indices file, volume fraction\n')
            mcfost_file.write('  {:} '.format(int(self.grain3_heating_method)) + '	                  Heating method : 1 = RE + LTE, 2 = RE + NLTE, 3 = NRE\n')
            values = '  {:.3} '.format(float(self.grain3_amin )) + '{:.3} '.format(float(self.grain3_amax)) +\
                '{:.3} '.format(float(self.grain3_aexp)) + '{:} '.format(int(self.grain3_n_grains)) 
            mcfost_file.write(values + ' 	  amin, amax [mum], aexp, n_grains (log distribution)\n')
        mcfost_file.write('\n')

        # Molecular RT settings. This was fast! :)
        mcfost_file.write('#Molecular RT settings\n'+\
            '  T T T 15.	          lpop, laccurate_pop, LTE, profile width (km.s^-1)\n'+\
            '  0.2 			  v_turb (delta)\n'+\
            '  1			  nmol\n'+\
            '  co@xpol.dat 6           molecular data filename, level_max\n'+\
            '  1.0 20     	  	  vmax (km.s^-1), n_speed\n'+\
            '  T 1.e-6 abundance.fits.gz   cst molecule abundance ?, abundance, abundance file\n'+\
            '  T  3                       ray tracing ?,  number of lines in ray-tracing\n'+\
            '  1 2 3	 		  transition numbers\n')
        mcfost_file.write('\n')

        # Star properties
        mcfost_file.write('#Star properties\n')
        # star 1, always present
        mcfost_file.write('  {:} '.format(int(self.n_stars)) +' Number of stars\n')
        values = '  {:.3} '.format(float(self.star1_teff)) + '{:.3} '.format(float(self.star1_radius)) + '{:.3} '.format(float(self.star1_mass)) +\
                   ' {:.3} '.format(float(self.star1_x)) + '{:.3} '.format(float(self.star1_y)) + '{:.3} '.format(float(self.star1_z)) + ' '+ self.star1_blackbody
        mcfost_file.write(values + '      Temp, radius (solar radius),M (solar mass),x,y,z (AU), is a blackbody?\n')
        mcfost_file.write('  ' + self.star1_spectrum +'\n')
        values = '  {:.3} '.format(float(self.star1_fUV)) + '{:.3} '.format(float(self.star1_slopefUV))
        mcfost_file.write(values + '  fUV, slope_fUV\n')

                 
        mcfost_file.close()


        
def generate_mcfost_file(path,params):

    """ This function takes a path and a dictionary with the parameters and creates the Paramfile class    """
    
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


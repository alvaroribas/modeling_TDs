# 
#  This file is part of Herschel Common Science System (HCSS).
#  Copyright 2001-2014 Herschel Science Ground Segment Consortium
# 
#  HCSS is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as
#  published by the Free Software Foundation, either version 3 of
#  the License, or (at your option) any later version.
# 
#  HCSS is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU Lesser General Public License for more details.
# 
#  You should have received a copy of the GNU Lesser General
#  Public License along with HCSS.
#  If not, see <http://www.gnu.org/licenses/>.
# 
###########################################################################
###                        SPIRE Map Merging Script                     ###
###########################################################################
# Title: Photometer_MapMerge.py
#
# Purpose:
# This script will merge two or more observations performed in SPIRE Large Map or 
# Small Map or Parallel mode to produce a single map. 
#
# Usage:
# 1) Load multiple observations (require obs ID and Pool name for each observation),
# 2) Collect the level 1 products (scan lines) together from each observation, 
# 3) Perform baseline subtraction on level 1 timelines.
# 4) Perform Map Making on scan lines from all observations for each SPIRE array.
#
# Required Inputs:
# a) Array of Observation IDs
# b) Corresponding array of Pool names (if not using the Herschel Science Archive)
# c) output directory for writing the maps FITS files
#
# Assumptions: A data pool has already been created on disk. 
#               The data has already been processed to Level 1
# 
# Last Updated: 27 Feb 2014
#
#
#   Written for HIPE v.12.0
#    
###########################################################################
#
###########################################################################
###                         Required inputs                             ###
###########################################################################
# (A) List of OBSIDs in the form of an integer or hexadecimal (0x) number:
#     (Example data taken from OT1 User Workshop, ESAC, March 2011) 
# (B) List of the corresponding data Pools in your Local Store:
# (C) Specify the output directory for writing the maps FITS files:
#
obsids  = [1342213178, 1342213179]

outDir     = "/pcdisk/stark/aribas/Desktop/modeling_TDs/remaps_Cha/SPIRE/destriper/"

###########################################################################
###                            user Options                             ###
###########################################################################
# 1) Apply relative bolometer gains for extended emission maps: Default = True
# 2) Use Destriper (else use simple median baseline removal): Default = True
# 3) Get Data from Herschel Science Archive (else get data from Local Pool): Default = False
#
# applyExtendedEmissionGains = True  : Map optimized for extended emission
# applyExtendedEmissionGains = False : Map optimized for Point Sources
applyExtendedEmissionGains = False     
useDestriper               = True
useHsa                     = True

# Applying relative bolometer gains requires chanRelGain Calibration Product
# Assume SPIRE Calibration Tree is present in your Local Store
if applyExtendedEmissionGains:
	cal          = spireCal(pool="spire_cal_12_3")
	chanRelGains = cal.phot.chanRelGain


###########################################################################
    
scanLines = 0

# Collect all scan lines from observation list.
for i in range(len(obsids)):
	print "Processing ObsID =", obsids[i],"("+hex(obsids[i])+")"
	#
	if useHsa:
		# Load observation from Herschel Science Archive
		obs = getObservation(obsids[i], useHsa=True, instrument="SPIRE")
	else:
		# Load observation from Local Pool
		obs = getObservation(obsids[i], poolName=pools[i],instrument="SPIRE")
	#
	# Fix for pre HIPE 10 version Level 2 map names
	prefix=""
	if obs.level2.refs['psrcPSW']!=None: prefix='psrc'
	#
	# Display PLW map of the input observations
    # To display a different band, change 2 instances of "PLW" in the line below
	Display(obs.refs["level2"].product.refs[prefix+"PLW"].product, title= "%i PLW Map"%obs.obsid)
	#
	if scanLines == 0:
		# Copy level 1 product from 1st observation.
		scanLines = obs.level1
                print obs.level1.count, "Scan lines from  observation ", obsids[i]
	else:
		for ref in obs.level1.refs:
			# Attach level 1 product from next observation.
			scanLines.addRef(ref)
                print obs.level1.count, "Scan lines from  observation ", obsids[i]

###########################################################################

# Apply relative gains for bolometers for better extended maps
if applyExtendedEmissionGains:
    #
    for i in range(scanLines.getCount()):
        psp = scanLines.getProduct(i)
        psp.meta['type'].string = "PSP"
        scanLines.getRefs().set(i, ProductRef(\
            applyRelativeGains(psp, gains = chanRelGains)))
    print "Finished applying relative gains"
    print


###########################################################################

# Using new Level 1 context. Run baseline removal  as an input to the map making
arrays = ["PSW","PMW","PLW"]   # Edit to select fewer bands
pixelSize = [6,10,14]          # Map pixel size in arcsec for PSW, PMW, PLW respectively
maps  = {}                     #  Empty Dictionary for maps  
diagnostics = {}               #  Empty Dictionary for Destriper diagnostic product 



# ####   Run the Destriper
if useDestriper:
  # Using all scan lines. Run destriper to remove offsets and make combined maps
  for iArray in range(len(arrays)):
      scanLines,map,diag,p4,p5 = destriper(level1=scanLines,\
          pixelSize=pixelSize[iArray], array=arrays[iArray], l2DeglitchRepeat=0)
      diagnostics[arrays[iArray]] = diag    # Keep diagnostic products
  print "Finished the Destriper run"
# 
# ####   or Run the Median Baseline Subtraction 
else:
# Using Level 1 context. Run median offset removal and make combined maps
  for iArray in range(len(arrays)):
      scanLines=baselineRemovalMedian(scanLines)
  print "Finished Median Baseline Removal"

###########################################################################


#Create, Display and save the final maps
for iArray in range(len(arrays)):
	# Create Map
	maps[arrays[iArray]] = naiveScanMapper(scanLines, array=arrays[iArray],resolution=pixelSize[iArray])
	# Display Map
	Display(maps[arrays[iArray]], title='Combined %s Map'%(arrays[iArray]))
	#Save map as a FITS file
	simpleFitsWriter(maps[arrays[iArray]], outDir+"%s_combined.fits"%(arrays[iArray]))
	simpleFitsWriter(diagnostics[arrays[iArray]], outDir+"%s_diagnostic.fits"%(arrays[iArray]))
print "Map saved as FITS files to %s"%(outDir)



print
print "********* End of Baseline Removal / Destriper Examples Script *********"
print
##########################  END  OF SCRIPT  #################################




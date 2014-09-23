### Main script to produce photometry
# This is a HIPE script, most python libraries cannot be included (e.g. numpy)

# read the table somehow
ra_sources = []
dec_sources = []

# total number of sources
n_objects = len(ra_sources)

# define some useful lists
instruments = ['PACS','SPIRE']
maps_pacs = ['scanamorphos','jscanam','unimap']
maps_spire = ['scanamorphos','destriper','unimap']
psfs = {'70':6., '100':8., '160':12.,\
        '250':18., '350':25., '500':36.}
recommended_apertures = {'70':[12., 20., 30.], '100':[12., 20., 30.], '160':[22., 30., 40.],\
                              '250':[22., 60., 90.], '350':[30., 60., 90.], '500':[42., 60., 90.]}

for index_object in 
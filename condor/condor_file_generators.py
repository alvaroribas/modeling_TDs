# Functions to produce the condor files.
# We will generate a whole condor file per generation.

def condor_header_generator(obj,path_file,path_generation,generation_number,kind,path_mcfost_executable):
    """ This function creates the header of a condor file for a MCMC generation """
    ### BEWARE!!!! ###
    ## executable and initialdir are still to be defined!!
    condor_file = open(path_file,'w')
    condor_file.write("# Condor submission for modeling TDs with MCFOST\n")
    condor_file.write("# A. Ribas\n")
    condor_file.write("# Last modification: "+str(date.today().day)+"/"+str(date.today().month)+"/"+str(date.today().year)+"\n")
    condor_file.write("# Update initialdir to match your desired path\n")
    condor_file.write("# Object: %s \n" %obj)
    condor_file.write("# Generation: %i \n" %generation_number)
    condor_file.write("# Kind: %s \n" %kind)
    condor_file.write('executable              = '+path_mcfost_executable+'\n')
    condor_file.write("getenv                  = True\n")
    condor_file.write("universe                = vanilla\n")
    #condor_file.write("requirements            = (machine != \"kool.cab.inta-csic.es\")\n")
    condor_file.write("transfer_executable     = False\n")
    condor_file.write("should_transfer_files   = False\n")
    condor_file.write("nice_user               = True\n")
    condor_file.write("notification 	          = Never\n" )
    condor_file.write("log                     = " +obj+".log\n")
    condor_file.write( 'initialdir              = '+path+"models/\n")
    condor_file.close()
    
def condor_include_object(obj,path_file,path_generation,generation_number,kind):
    """ This function adds the corresponding lines to add an object to the condor file for MCFOST"""
    ### BEWARE!!!!
    ## all these still to be decided, not sure if needed
    condor_file = open(path_file,'a')
    condor_file.write( 'initialdir              = '+path+"\n")
#   condor_file.write("arguments               = \n")
    condor_file.write("output                  = model%i_generation%i.out\n", %(model_number,generation_number))
    condor_file.write("error                   = model%i_generation%i.error\n", %(model_number,generation_number))
    condor_file.write("queue\n")
    condor_file.close()

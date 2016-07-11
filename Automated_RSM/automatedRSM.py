import math
import os
import numpy as np
import imp
import glob
lhsu = imp.load_source("lhsu", "./LHS/lhsu.py")

'''Push button generation of a RSM for the CD of a satellite'''


###################### CLEAN-UP TPM DIAGNOSTIC FILES #######################
def tpm_cleanup(tpmdir):
    diagfiles = glob.glob(tpmdir+"/tpmout*")
    if not diagfiles == []:
        cmd = "rm "+tpmdir+"/tpmout*"
        os.system(cmd)
 

#################### CHANGE "ADDPATH" FOR MATLAB CODES #####################
def matlab_addpath(mcmcdir, gpmsadir):
    #Define addpath string
    addpath = "addpath('"+gpmsadir+"');\n"

    #Open master.m
    masterpath = os.path.join(mcmcdir, "master.m")
    fmaster = open(masterpath, "r")

    #Read master.m
    lines = []
    for line in fmaster:
        lines.append(line)
    fmaster.close()

    #Change addpath
    for i in range(len(lines)):
        if "addpath" in lines[i]:
            lines[i] = addpath

    #Rewrite file
    fmaster = open(masterpath, "w")
    for line in lines:
        fmaster.write(line)
    fmaster.close()  


#################### CHANGE "ADDPATH" FOR MATLAB CODES #####################
def matlab_makeRSM(NENS, mcmcdir):
    #Define makeRSM modified string
    NENSstr = "M = "+str(NENS)+";\n"

    #Open makeRSM.m
    makeRSMpath = os.path.join(mcmcdir, "makeRSM.m")
    fRSM = open(makeRSMpath, "r")

    #Read makeRSM.m
    lines = []
    for line in fRSM:
        lines.append(line)
    fRSM.close()

    #Change addpath
    for i in range(len(lines)):
        if "% Gather some stuff together" in lines[i]:
            lines[i+1] = NENSstr

    #Rewrite file
    fRSM = open(makeRSMpath, "w")
    for line in lines:
        fRSM.write(line)
    fRSM.close()          


########################### READ INPUT PARAMETERS ############################
def modify_input(NENS, tpmdir, species, D2R):
    #Open the input file
    inputpath = os.path.join(tpmdir, "tpm.inp")
    fin = open(inputpath, "r")

    #Read all the lines and store them
    lines = []
    for line in fin:
        lines.append(line)

    #Close the file
    fin.close()

    #Define search strings
    searchnames = ["Number of Ensemble Points", "Gas Surface Interaction Model", "Species Mole Fractions"]

    #Define Number of ensemble members and Gas-Surface Interaction Model
    for iline, line in enumerate(lines):
        data = line.strip().split("#")
        if searchnames[0] in line:
            data[1] = "# %d\n" % NENS
            lines[iline] = data[0] + data[1]
        if searchnames[1] in line:
            GSI_MODEL = int(data[1])
        if searchnames[2] in line:
            data[1] = "#  %lf %lf %lf %lf %lf %lf\n" % (species[0], species[1], species[2], 
                                                        species[3], species[4], species[5])
            lines[iline] = data[0] + data[1]

    #Change mole fraction array
    inputpath = os.path.join(tpmdir, "tpm.inp")
    fin = open(inputpath, "w")
    for line in lines:
        fin.write(line)
    fin.close()         

    #Define variable names
    variablenames = ["Magnitude of Bulk Velocity", "Yaw Angle Orientation", "Pitch Angle Orientation",
                     "Satellite Surface Temperature", "Atmospheric Translational Temperature"]

    #Append GSI specific variables
    if GSI_MODEL == 0:
        variablenames.append("Specular Fraction")
        NVAR = 6
    if GSI_MODEL == 1:
        variablenames.append("Energy Accommodation Coefficient")
        NVAR = 6
    if GSI_MODEL == 2:
        variablenames.append("Normal Energy Accommodation Coefficient")
        variablenames.append("Tangential Momentum Accommodation Coefficient")
        NVAR = 7

    #Initialize min and max arrays
    xmin = np.zeros(NVAR)
    xmax = np.zeros(NVAR)

    #Read the file and assign min and max arrays
    for iline, line in enumerate(lines):
        data = line.strip().split("#")
        if len(data) > 1:
            data = data[1].split()
            for i, name in enumerate(variablenames):
                if name in line:
                    xmin[i] = float(data[0])
                    xmax[i] = float(data[1])
                    break

    #Convert degrees to radians
    xmin[1] *= D2R
    xmax[1] *= D2R
    xmin[2] *= D2R
    xmax[2] *= D2R

    return xmin, xmax, GSI_MODEL

###################### GENERATE LHS ENSEMBLE  ##########################
def lhs_ensemble(tpmdir, xmin, xmax, NENS, GSI_MODEL):

    #CREATE LATIN HYPERCUBE SAMPLE
    LHS = lhsu.lhsu(xmin, xmax, NENS)

    #OPEN ENSEMBLE FILE
    ensemblepath = os.path.join(tpmdir, "tpm.ensemble")
    f = open(ensemblepath, "w")

    #DEFINE HEADER
    if GSI_MODEL == 0:
        header = "Umag [m/s]     theta [radians]     phi [radians]       Ts [K]     Ta [K]     epsilon"
    if GSI_MODEL == 1:
        header = "Umag [m/s]     theta [radians]     phi [radians]       Ts [K]     Ta [K]     alpha"
    if GSI_MODEL == 2:
        header = "Umag [m/s]     theta [radians]     phi [radians]       Ts [K]     Ta [K]     alphan     sigmat"

    #WRITE ENSEMBLES
    np.savetxt('./tpm.ensemble', LHS, header=header)


######################## RUN TPMC ENSEMBLE #############################
def run_tpmc(NPROCS, speciesnames, ispec, rtype, outdir):

    #RUN THE CODE
    print "Starting simulation for "+speciesnames[ispec]+" "+rtype+" set\n"
    cmd = "mpiexec -n "+str(NPROCS)+" ./tpm"
    os.system(cmd)

    #COPY THE OUTPUT TO A NEW FILE
    print "Copying output data\n"
    outname = "CYGNSS_"+speciesnames[ispec]+"_"+rtype+"_set.dat"
    outpath = os.path.join(outdir, outname)
    cmd = 'cp Cdout.dat "'+outpath+'"'
    os.system(cmd)



############### LOOP OVER TEST/TRAINING AND MOLE FRACTIONS ##############
def tpmc_loop(NPROCS, NENS, mcmcdir, tpmdir):

    ####################### CONSTANTS #######################
    D2R = (math.pi/180.0)   #DEGREES TO RADIANS CONVERSION
    NSPECIES = 6            #NUMBER OF SPECIES
    speciesnames = ["O", "O2", "N", "N2", "He", "H"]
    species = np.zeros(NSPECIES)

    os.chdir(tpmdir)

    for idx in range(2):
        for ispec in range(NSPECIES):
            #Define training and test sets
            if idx == 0:
                rtype = "training"
                outdir = os.path.join(mcmcdir, "data/Training Set/")
                
            if idx ==1 :
                rtype = "test"
                outdir = os.path.join(mcmcdir, "data/Test Set/")
            
            #Create directory if it doesn't exist
            if not os.path.exists(outdir):
                    os.makedirs(outdir)

            #Define Species Array
            for i in range(NSPECIES):
                species[i] = 0.0
                if i == ispec:
                    species[i] = 1.0
      
            #Read input
            xmin, xmax, GSI_MODEL = modify_input(NENS, tpmdir, species, D2R)

            #Create LHS ensemble
            lhs_ensemble(tpmdir, xmin, xmax, NENS, GSI_MODEL)

            #Run TPMC code
            run_tpmc(NPROCS, speciesnames, ispec, rtype, outdir)

########################## RUN THE MCMC FITTING #########################
def run_mcmc(mcmcdir):

    #Clean-up parent directory
    cmd = "rm ./Cdout*.dat"
    os.system(cmd)

    #Check for required directories; create them if they don't exist
    poutdir = os.path.join(mcmcdir, "pout")
    if not os.path.exists(poutdir):
        os.makedirs(poutdir)
    rhodir = os.path.join(mcmcdir, "rho")
    if not os.path.exists(rhodir):
        os.makedirs(rhodir)
    testpreddir = os.path.join(mcmcdir, "testpred")
    if not os.path.exists(testpreddir):
        os.makedirs(testpreddir)

    #Spawn master.m
    os.chdir(mcmcdir)
    cmd = "matlab -nodisplay -nosplash -r master"
    os.system(cmd)


################# MOVE AND RENAME RSM OUTPUT DATE FILE ###################
def move_output(RSMNAME, basedir, mcmcdir):

    oldRSMpath = os.path.join(mcmcdir, "RSM.dat")
    newRSMpath = os.path.join(basedir, RSMNAME+"_RSM.dat")

    cmd = "mv "+oldRSMpath+" "+newRSMpath
    os.system(cmd)


##########################################################################
##################### RUN AUTOMATED RSM GENERATION #######################
##########################################################################

if __name__ == '__main__':

    basedir = os.getcwd()
    tpmdir  = os.path.join(basedir, "TPM")
    mcmcdir = os.path.join(basedir, "MCMC")
    gpmsadir = os.path.join(mcmcdir, "gpmsa/matlab")
    
    ####################### INPUTS ##########################
    NENS = 1000             #NUMBER OF ENSEMBLE MEMBERS
    NPROCS = 20             #NUMBER OF PROCESSORS FOR SIMULATION
    RSMNAME = "CYGNSS"      #NAME OF OUTPUT RSM FILE

    ##################### START CODE ########################

    tpm_cleanup(tpmdir)                        #Cleanup TPMC diagnostic files (If it crashed last time)
    matlab_addpath(mcmcdir, gpmsadir)          #Change matlab "addpath" directory
    matlab_makeRSM(NENS, mcmcdir)              #Change number of ensembles in RSM write-out
    tpmc_loop(NPROCS, NENS, mcmcdir, tpmdir)   #Run TPMC test and training sets
    run_mcmc(mcmcdir)                          #Run MCMC fitting to data
    move_output(RSMNAME, basedir, mcmcdir)     #Move final output to base directory
    tpm_cleanup(tpmdir)                        #Cleanup TPMC diagnostic files (Again)

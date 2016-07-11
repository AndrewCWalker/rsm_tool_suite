import math
import os
import numpy as np
import subprocess
import glob


'''Computes satellite projected area as a function of orientation'''


############# COMPUTE NUMBER OF SIMULATION POINTS ###############
def compute_numpoints(phi_min, phi_max, theta_min, theta_max, dphi, dtheta):

    #COMPUTE TOTAL ANGLE RANGE
    deltaphi = phi_max - phi_min
    deltatheta = theta_max - theta_min

    #COMPUTE NUMBER OF ANGLE POINTS IN EACH DIRECTION
    Nphi = int(deltaphi/dphi)+1
    Ntheta = int(deltatheta/dtheta)+1

    return Nphi, Ntheta


####################### WRITE ENSEMBLE FILE #####################
def write_ensemble(areadir, phi_min, phi_max, theta_min, theta_max, Nphi, Ntheta, NUMPOINTS, MESHNAME):

    ##### CONSTANTS #####
    R2D = math.pi/180.0

    #CREATE ANGLE ARRAYS
    phi = np.zeros(NUMPOINTS)
    theta = np.zeros(NUMPOINTS)

    #CONVERT TO RADIANS
    phi_min *= R2D
    phi_max *= R2D
    theta_min *= R2D
    theta_max *= R2D

    #ASSIGN ANGLE ARRAYS
    for i in range(Nphi):
        for j in range(Ntheta):
            k = i*Ntheta + j 
            phi[k] = phi_min + (phi_max - phi_min)/(Nphi-1)*i
            theta[k] = theta_min + (theta_max - theta_min)/(Ntheta-1)*j

    #OPEN FILE
    ensemblepath = os.path.join(areadir, "area.ensemble")
    f = open(ensemblepath, "w")

    #WRITE HEADER
    meshpath = os.path.join(areadir, "Mesh_Files", MESHNAME)
    f.write("Number of Ensemble Members          # %d\n" % NUMPOINTS)
    f.write("Mesh Filename                       # %s\n" % meshpath)
    f.write("Umag [m/s]      theta [radians]     phi [radians]\n")

    #WRITE ANGLES
    for k in range(NUMPOINTS):
        f.write("7500.0          %1.10e          %1.10e\n" % (theta[k], phi[k]))

    #CLOSE FILE
    f.close()


########################### RUN AREA CODE #######################
def run_area_code(areadir, NPROCS):
    os.chdir(areadir)
    cmd = "mpiexec -n {0} ./area".format(NPROCS)
    os.system(cmd)


########################### RENAME OUTPUT #######################
def rename_output(basedir, SATNAME):
    oldpath = os.path.join(basedir, "Aout.dat")
    newfilename = SATNAME+"_area.dat"
    newpath = os.path.join(basedir, newfilename)
    cmd = "mv "+oldpath+" "+newpath 
    os.system(cmd)


########################### CLEAUP OUTPUT #######################
def cleaup_output(areadir):
    outputfiles = glob.glob(areadir+"/Aout*.dat")
    if not outputfiles == []:
        cmd = "rm "+areadir+"/Aout*.dat"
        os.system(cmd)


########################### MAIN FUNCTION #######################

if __name__ == '__main__':

    basedir = os.getcwd()
    areadir = os.path.join(basedir, "Area_Code")
    
    ########## INPUT VARIABLES ##########

    MESHNAME = "CYGNSS_final_ascii.stl"   #MESH FILENAME
    SATNAME = "CYGNSS"                    #SATELLITE NAME
    NPROCS = 30                           #NUMBER OF PROCESSORS

    dphi = 0.1                            #PITCH ANGLE RESOLUTION (IN DEGREES)
    dtheta = 0.1                          #YAW ANGLE RESOLUTION (IN DEGREES)

    phi_min = 0.0                         #MINIMUM PITCH ANGLE (IN DEGREES)
    phi_max = 90.0                        #MAXIMUM PITCH ANGLE (IN DEGREES)
    theta_min = 0.0                       #MINIMUM YAW ANGLE (IN DEGREES)
    theta_max = 3.0                       #MAXIMUM YAW ANGLE (IN DEGREES)

    ######## END INPUT VARIABLES ########

    #COMPUTE NUMBER OF DATA POINTS
    Nphi, Ntheta = compute_numpoints(phi_min, phi_max, theta_min, theta_max, dphi, dtheta)
    NUMPOINTS = Nphi*Ntheta  #NOTE: MAKE SURE NUM_POINTS IN TPM.C MATCHES THIS VALUE
    
    #WRITE ENSEMBLE FILE
    write_ensemble(areadir, phi_min, phi_max, theta_min, theta_max, Nphi, Ntheta, NUMPOINTS, MESHNAME)

    #RUN AREA CODE
    run_area_code(areadir, NPROCS)

    #RENAME OUTPUT FILE
    rename_output(basedir, SATNAME)

    #CLEANUP OUTPUT FILES
    cleaup_output(areadir)

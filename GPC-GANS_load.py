## This is a python 3 implementation of GPC-JE to the GANS problem.

import sys
import json
import os
import itertools
import numpy

from CallPSGFunctions import *


# The purpose of this program is to define the parameter space, create configuration files, and submit GPC jobs to process those config files through PSG.

project = '121619/' # unique for each run. This is a drectory! Include a '/'!

# Here we begin assigning parameter values.
Rlist = [ 1.0 ] # Object radius in Earth radii
glist = [ 9.8 ] # Surface gravity in m/s2
P0list = [ -1.0,-0.5,]# 0.0,0.5, 1.0,1.5, 2.0 ] # Surface pressure in log(bars)
Aslist = [ 0.05, 0.2,]# 0.35, 0.5 ]  # Surface grey albedo. Dimensionless

cloudstate = ['clouds',]#'haze','clear'] # Name the cloud state. Include 'clouds' for clouds, 'clear' for clear sky, 'haze' for faux haze
rcldlist = [ 10 ] # Radius of cloud particulates in microns
lcldlist = [ -3.0,]#-2.5,-2.0,-1.5,-1.0 ] # Linear cloud density in log(kg/m2) 
Ptoplist = [ -2.0,]#-1.5,-1.0,-0.5,0.0, 0.5, 1.0 ] # Pressure level of the top of the cloud in bars
fadelist = [ 0.1] # The rate at which the clouds fade outside of the "core" layers. Dimensionless.

gasNames = ['H2O']#,'CO2','CH4','O2','SO2'] # names of gasses
gasTypes = ['HIT[1]']#,'HIT[2]','HIT[6]','HIT[7]','HIT[9]'] # database to use for gas line features. one entry for each entry in gasNames
gasMixRat = [ -7.0]#,-6.0, -5.0, -4.0, -3.0, -2.0, -1.0, -0.01 ] # log(mixing ratio)

H2mix = [ 0.0]#,0.02, 0.05, 0.15, 0.5 ] # The quantity H2/(H2+N2). Dimensionless.

T = 300

gasAbundanceList = []

for i,j,k in itertools.product(range(len(gasNames)),range(len(gasMixRat)),range(len(H2mix))): # constructing appropriate mixing ratio assemblies
    gasAbundanceList.append({'gas':gasNames[i],'gasType':gasTypes[i],'gasMix':10**gasMixRat[j],'H2Type':'HIT[45]','N2Type':'HIT[22]','H2Mix':(1-10**gasMixRat[j])*H2mix[k],'N2Mix':(1-10**gasMixRat[j])*(1-H2mix[k])})

configList = []

print('Processing...')
print()

iterator = 0 # This is just a counter, to track how many loops we've been through.

import os
# Creating all the directories required for computation.
if not os.path.exists('Cfgs/'+project):
    os.mkdir('Cfgs/'+project)


# Beginning creation of config files and GPC jobs.
# Handle cloudy atmospheres.
if ('clouds' in cloudstate):

    clouds = 'clouds'
    print('Starting clouds...')
    print()
    
    # loop over possible atmosphere configurations
    for R,g,P0,As,atm,rcld,lcld,Ptop,fade in itertools.product(Rlist,glist,P0list,Aslist,gasAbundanceList,rcldlist,lcldlist,Ptoplist,fadelist):
        if Ptop >= P0: # catch impossible cloud locations
            continue

        # construct the config file name
        filename = project # start with project name
        filename += atm['gas']+str(round(numpy.log10(atm['gasMix']),3))+'H2'+str(round(atm['H2Mix']/(1-atm['gasMix']),3))+'_R'+str(R)+'g'+str(g)+'P'+str(P0)+'A'+str(As) # add most parameter values
        filename += 'r'+str(rcld)+'l'+str(lcld)+'ct'+str(Ptop)+'f'+str(fade)+'.cfg' # add cloud parameter values and filetype name
        
        # move pressure parameters out of log space
        P0 = round(10**P0,3) 
        Ptop = round(10**Ptop,3)

        try: # running the configuration file maker, and adding the output to our future joblist
            sys.stdout.write(' '+str(iterator)+': '+MakeCFG(filename,[R,g,P0,As,atm,rcld,lcld,Ptop,clouds,fade,T])+'        ')
            sys.stdout.write('\r')
            sys.stdout.flush()
            configList.append(filename)
        except BadClouds: # BadClouds should only come up from MakeCFG, if the cloud top pressure is too high (i.e. would produce underground clouds)
            continue
        iterator = iterator+1
    print()
    print('Clouds finished!') # Status tracker



# Handle clear atmospheres.
if 'clear' in cloudstate:
    clouds = 'clear'
    print('Starting clearsky...')
    print()

    # loop over possible atmosphere configuraitons.
    for R,g,P0,As,atm in itertools.product(Rlist,glist,P0list,Aslist,gasAbundanceList):

        # construct the config file name
        filename = project # start with project name
        filename += atm['gas']+str(round(numpy.log10(atm['gasMix']),3))+'H2'+str(round(atm['H2Mix']/(1-atm['gasMix']),3))+'_R'+str(R)+'g'+str(g)+'P'+str(P0)+'A'+str(As) # add the parameter values

        filename += '.cfg' # tack on the file type name
        
        P0 = round(10**P0,3) # move pressure out of log space
        # create the configuration file
        sys.stdout.write(' '+str(iterator)+': '+MakeCFG(filename,[R,g,P0,As,atm,1,1,1,clouds,1,T])+'        ')
        sys.stdout.write('\r')
        sys.stdout.flush()
        configList.append(filename) # add it to the future joblist
        iterator = iterator+1
    print()
    print('Clearsky finished!') # Status tracker


# Handle "hazy" atmospheres.
if 'haze' in cloudstate:
    clouds = 'haze'
    print('Starting haze...')
    print()

    # loop over possible atmosphere configurations
    for R,g,P0,As,atm,rcld,lcld,fade in itertools.product(Rlist,glist,P0list,Aslist,gasAbundanceList,rcldlist,lcldlist,fadelist):

        # construct the config file name
        filename = project # start with project name
        filename += atm['gas']+str(round(numpy.log10(atm['gasMix']),3))+'H2'+str(round(atm['H2Mix']/(1-atm['gasMix']),3))+'_R'+str(R)+'g'+str(g)+'P'+str(P0)+'A'+str(As) # add most parameter values
        filename += 'r'+str(rcld)+'l'+str(lcld)+'f'+str(fade) # add aerosol parameter values

        filename += '.cfg' # add filetype name
        
        # move pressure out of log space    
        P0 = round(10**P0,3)
        # create the configuration file
        sys.stdout.write(' '+str(iterator)+': '+MakeCFG(filename,[R,g,P0,As,atm,rcld,lcld,1,clouds,fade,T])+'        ')
        sys.stdout.write('\r')
        sys.stdout.flush()
        configList.append(filename) # add it to the future joblist
        iterator = iterator+1

    print()
    print('Haze finished!') # Status tracker

# Congrats, we've finished making our configuration files!
print()
print(iterator,'Config Files Prepared!') # How many did we make?
quit()
print('Building the job list. Please stand by, this can take a while...')
print()

import python_client
url = os.getenv('JE_HOST','http://localhost:8007/') # where to find PSG, should point to a locally installed Docker
jec = python_client.GpcJobEngineClient(url) # for GPC

if not os.path.exists('Data/'+project):
    os.mkdir('Data/'+project)
    if 'clouds' in cloudstate:
        os.mkdir('Data/'+project+'Clouds/')
    if 'clear' in cloudstate:
        os.mkdir('Data/'+project+'NoClouds/')
    if 'haze' in cloudstate:
        os.mkdir('Data/'+project+'Haze/')

i = 0

for configFile in configList: # Running over our job/payload list
    # This little code block just makes a spinny animation while it works.
    if i == 0:
        j = '  -  ' 
    if i == 1:
        j = '  /  '
    if i == 2:
        j = '  |  '
    if i == 3:
        j = '  \\  '
    i = (i+1)%4
    sys.stdout.write('\r')
    sys.stdout.write(j)
    sys.stdout.flush()
    # End animation code block.
    payload = {}
    payload['name'] = configFile # Define the payload.
    jec.addJob(json.dumps(payload), configFile, jobFamily='GPC-PSG-GANS') # Create the GPC job. Be careful mucking with this line, suggest you talk to GPC people before changing anything here.

print('All jobs submitted.')
# Alldone, goodbye

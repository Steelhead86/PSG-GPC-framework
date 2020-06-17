
import shutil
import copy
import math
import os
import astropy.units as u
import scipy.constants as cons
import numpy as np
from tempfile import mkstemp


class BadClouds(Exception):
    '''
    This exception exists soley for diagnostic purposes.
    '''
    pass


def NormalizeAtmosphereAbund(abundances):
    '''
    This very basic function exists to verify that the sum of mixing ratios for each element in the atmosphere is equal to 100%.
    If it is not, the function will correct it while mantaining the ratios.
    '''
    total = sum([ float(i) for i in abundances ])
    if total != 1:
        normalization = 0.01*total
        newAbunds = [ str(float(i)/normalization) for i in abundances ]
    else:
        newAbunds = abundances
    return newAbunds


def MakeCFG(file_path, parameters):
    '''
    The core function for constructing configuration files based on a template.
    The file_path variable should be a string dictating what we want the new file to be named.
    The other variables are parameters to be fed into the configuration file.
    '''
    
    # Doing some pre-computation before moving into the read/write.
    # This is to save processor time later.
    #
    # Extracting parameter values from the bulk parameter variable.
    P0  = parameters[2]
    atm = parameters[4]

    # Declare number of layers in our initial atmosphere. This may dynamically change with clouds.
    nlayersinit = 51

    # High density atmospheres require more layers. Therefore, we increase the declared nlayers based on surface pressure.
    if P0 >= 10:
        nlayersinit += 20
    if P0 >= 30:
        nlayersinit += 35
    if P0 >= 100:
        nlayersinit += 50
    
    # Define the molecular mass of our absorbing chemicals. The units are grams per mole.
    molMass = {'H2O':18.01528, 'CO2':44.009, 'SO2':64.066, 'O2':31.9988, 'CH4':16.043, 'N2':28.0134, 'H2':2.01588 }

    # Use this to calculate the mean molecular weight of the atmosphere in grams per mole.
    mmw = atm['gasMix']*molMass[atm['gas']] + atm['H2Mix']*molMass['H2'] + atm['N2Mix']*molMass['N2']



    # We now begin actual file construction.
    #
    # Define the name of our default .cfg. This should be compatible with the current PSG version.
    cfgdefault = 'template.cfg'
    # Duplicate the default configuration file.
    directory = './Cfgs/'
    shutil.copy(cfgdefault,directory+file_path) 
    # Create temporary working file.
    fh, abs_path = mkstemp()

    # We now open both the temporary and the template configuration file, and run through the template line-by-line. Each line is fed through the
    # ScanLine() subroutine to determine what needs to be done to it.
    with os.fdopen(fh,'w') as new_file: # Opening the working file
        with open(directory+file_path) as old_file: # Opening the duplicate template.
            # Go through each line of the configuration file.
            for line in old_file:
                ScanLine(new_file,line,parameters,nlayersinit,mmw) # Run this line through 'Scanline()' to change it as required
                                                                                       # by the parameters.

    # We now overwrite the duplicate configuration file with our modified version, and finish the function.
    shutil.move(abs_path, directory+file_path)
    return file_path


def ScanLine(file,line,parameters,nlayersinit,mmw):
    '''
    This function reads an output file, an input line (string), and a set of parameters. It checks the line to see if it needs to be modified,
    makes any appropriate modifications, and prints it to the output file.
    Note that the current incarnation is hard coded to use an atmosphere of N2, H2, and 1 other gaseous component, as well as up to one aerosol type.
    '''

    # Extract all parameter values from the 'parameters' value list.
    (Rp,g,P0,As,atm,rcld,lcld,Ptop,clouds,fade,T) = (*parameters,)

    # What follows are a number of 'if' statements, to look for and replace specific lines in the text. They are encased in a "loop" to allow use
    # of the 'break' command.
    # Future consideration: This could be made more compact with a dict-type object.
    for _ in [1]:
        # First category: Atmospheres.
        if '<ATMOSPHERE-LAYERS-MOLECULES>' in line: # This line indicates the start of the Atmosphere layers definition.
            # The Atmosphere-Layers line construction is involved, so we handle it in a separate function.
            line = AtmoLineCreation(nlayersinit,P0, atm, rcld, lcld, Ptop, clouds, fade, mmw, T)
            # We replaced a line, so we set this to True. All the other flags will get skipped.
            break 

        if '<ATMOSPHERE-NGAS>' in line: # ATMOSPHERE-NGAS tells the program how many gasses are present.
            line = '<ATMOSPHERE-NGAS>'+'3\n' # We have 3.
            break
        
        # These atmosphere lines relate to aerosols.
        #
        # The global abundance of atmosphere aerosols. With scaler set in AABUN, this is a scaler multiplier on the abundances we wrote for clouds above.
        if '<ATMOSPHERE-AABUN>' in line and clouds != 'clear':
            line = '<ATMOSPHERE-AABUN>1\n' 
            break

        # The units of each aerosol unit written in AABUN. Should be set to scalar ('scl')
        if '<ATMOSPHERE-AUNIT>' in line and clouds != 'clear':
            line = '<ATMOSPHERE-AUNIT>scl\n'
            break

        # The radius of aerosol particles.
        if '<ATMOSPHERE-ASIZE>' in line and clouds != 'clear':
            line = '<ATMOSPHERE-ASIZE>'+str(rcld)+'\n' 
            break

        # The number of aerosol types in the atmosphere. Note that this is either 1 or 0.
        if '<ATMOSPHERE-NAERO>' in line:
            # A catch to determine if there are clouds or not.
            if clouds == 'clear':
                naero = 0
            else:
                naero = 1 # Currently, we only handle one-aerosol-type cases.
            line = '<ATMOSPHERE-NAERO>'+str(naero)+'\n' # Set the numer of aerosols
            break

        # Non-cloud atmosphere elements
        #
        # Set the mean molecular weight of the atmosphere.
        if '<ATMOSPHERE-WEIGHT>' in line:
            line = '<ATMOSPHERE-WEIGHT>'+str(mmw)+'\n'
            break

        # Name the gasses (not aerosols!) present in the atmosphere.
        if '<ATMOSPHERE-GAS>' in line:
            line = '<ATMOSPHERE-GAS>'+atm['gas']+',H2,N2'
            line = line+'\n'
            break
        
        # Name the line lists to be used for each gas named in ATMOSPHERE-GAS.
        if '<ATMOSPHERE-TYPE>' in line:
            line = '<ATMOSPHERE-TYPE>'+atm['gasType']+','+atm['H2Type']+','+atm['N2Type']
            line = line+'\n'
            break

        # As with AABUN, this is the global abundance of the gasses. With ATMOSPHERE-UNIT set to 'scl' this acts as a multiplier on the abundances
        # written above.
        if '<ATMOSPHERE-ABUN>' in line:
            line = '<ATMOSPHERE-ABUN>'
            for abun in range(3):
                line +='1,' # Should be 1 for everything; these set scalar relations to the quantities defined in layers.
            line += '\n'
            break
        
        # As with AUNIT above, this is the units of each value written in ABUN. Should be set to scalar ('scl').
        if '<ATMOSPHERE-UNIT>' in line:
            line = '<ATMOSPHERE-UNIT>'
            for unit in range(3):
                line +='scl,' # Should be 'scl' to cooperate with ATMOSPHERE-ABUN above.
            line +='\n'
            break

        # Surface pressure of the atmosphere.
        if '<ATMOSPHERE-PRESSURE>' in line:
            line = '<ATMOSPHERE-PRESSURE>'+str(P0)+'\n'
            break

        # Set the unit of pressure for the atmosphere.
        if '<ATMOSPHERE-PUNIT>' in line:
            line = '<ATMOSPHERE-PUNIT>bar'+'\n'
            break

        # Second category: Object. Primarily bulk properties of the planet, such as size, mass, orientation.
        #
        # Surface gravity of the object is applied here.
        if '<OBJECT-GRAVITY>' in line:
            line = '<OBJECT-GRAVITY>'+str(g)+'\n'
            break

        # Diameter of the object here. Currently set to be a multiple of Earth radius.
        if '<OBJECT-DIAMETER>' in line:
            line = '<OBJECT-DIAMETER>'+str(2*Rp*u.earthRad.to(u.km))+'\n' 
            break

        # Third category: Surface. For handling surface properties of the planet.
        #
        # Isotropic grey albedo for the planet surface.
        if '<SURFACE-ALBEDO>' in line:
            line = '<SURFACE-ALBEDO>'+str(As)+'\n' 
            break

        # Inverse of SURFACE-ALBEDO.
        if '<SURFACE-EMISSIVITY>' in line:
            line = '<SURFACE-EMISSIVITY>'+str(1-As)+'\n' 
            break


    file.write(line) # We have modified the line; now write it to the file.


def AtmoLineCreation(nlayersinit,P0, atm, rcld, lcld, ctop, clouds, fade, mmw, T):
    # The first Atmospheres line declares the composite molecules of our atmosphere - that is, which component is in which column.
    line = '<ATMOSPHERE-LAYERS-MOLECULES>'+atm['gas']+',H2,N2,'

    # Before continuing, we must compute the atmosphere pressure grid.
    AGrid = BuildAtmosphereGrid(P0,nlayersinit)

    #  Catching certain problems with the atmosphere shape. Possibly obsolete?
    if ( type(AGrid) != type(np.array([0,1,2])) ) or ( ( type(AGrid[0]) == int ) or ( type(AGrid[0]) == float ) ):
        print('Bad atmo grid!                        ')
        print(type(AGrid))
        print(P0, atm, rcld, lcld, ctop, clouds)
        exit()


    # Identify if we need to handle clouds or not, and act appropriately.
    if clouds == 'clouds':
        # If this is a cloudy atmosphere, build clouds,
        AGrid = BuildClouds(ctop, lcld, AGrid, mmw, fade, T)
        # and append them to the declaration line.
        line = line+'Cloud'
    elif clouds == 'haze': # If this is a hazy atmosphere, build haze,
        AGrid = BuildHaze(lcld, AGrid, mmw, fade, T)
        line = line+'Cloud' # and add cloud element to declaration line.
    else: # Not cloudy, not hazy, must be clear.
        AGrid = [ [i,0] for i in AGrid ]

    # We are now done with the first (declaration) Atmospheres line. Add an endline and continue.
    line+='\n'
    
    # The number of atmosphere layers may have changed since initialization. Here we extract the final number of layers.
    nlayers = str(len(AGrid))
    # The second Atmosphere line declares the number of atmosphere layers. We create that text, append it to the line and include an endline.
    line+= '<ATMOSPHERE-LAYERS>'+nlayers+'\n'

    # We now loop over each layer in the atmosphere, adding its components to the line.
    for i in range(len(AGrid)): # For each layer in the atmosphere:

        # Error handling. Possibly obsolete?
        try:
            # Start the layer text. Write the layer number (starting with 1, at the bottom of the atmosphere), pressure, and temperature.
            line+= '<ATMOSPHERE-LAYER-'+str(i+1)+'>'+str(AGrid[i][0])+','+str(T)+',' 
        except:
            print()
            print('Bad Atmo Grid!')
            print(AGrid)
            print(AGrid.shape)
            print(P0,atm,rcld,lcld,ctop,clouds)
            exit()

        # Add gas mixing ratios to the text.
        line+= str(atm['gasMix'])+','+str(atm['H2Mix'])+','+str(atm['N2Mix'])+',' 

        # Include aerosols if appropriate.
        if clouds in ['clouds','haze']:
            line+=str(AGrid[i][-1])

        # Add an endline and continue to the next layer.
        line+=',\n' 

    # We have finished this process. Return the completed line.
    return line


def BuildAtmosphereGrid(P0,layers):
    '''
    This function constructs the pressure layer grid of the atmosphere. It is defined by a number of 'layers', and the pressure at the surface.
    '''

    # Initialize an empty grid
    grid = np.zeros(int(layers))
    
    # Declare the top-of-atmosphere pressure (in bars), then translate that number into scale heights (relative to P0)
    ToA = 10**-5
    maxScale = np.log(P0/ToA)

    # Given the number of layers declared, calculate how many scale heights per layer. Maximum of 1/3 of a scale height per layer.
    HperLayer = min(1./3., maxScale/layers)

    # Construct the grid.
    for i in range(len(grid)):
        height = HperLayer * i # The number of scale heights at this layer.
        grid[i] = P0*np.exp(-1*height) # The atmospheric pressure at this many scale heights.
    
    # Finished.
    return grid


def BuildHaze(totalMass,structure,mmw,fade,T):
    '''
    This function adds a fading "haze"-like cloud structure to an atmosphere.
    The haze is computed such that the column mass of cloud matter in the atmosphere sums to a predefined (totalMass) value.
    Since PSG treats aerosols in a mass mixing ratio (kg/kg) we need to understand the mass structure of the atmosphere to correctly construct the
    haze. Since this cloud structure is quite simple, mass-structure computation dominates this function.
    '''


    # Begin mass-structure computation. First, a bit of bookkeeping.
    #
    # Create a grid from the atmosphere grid, including temperature at each layer.
    structure = [ [i, T] for i in structure]
    # Extract the pressure levels from the atmosphere grid.
    # Pstruct = [ i[0] for i in structure ] # NOT USED RIGHT NOW


    # Second, we create a layer list in distance-space, rather than in pressure-space.
    #
    # Compute the length of one scale height (H), in meters. For 1 bar this should be about 8.5 km.
    SclHgt = (cons.R*T)/(mmw*cons.g)
    # For any given altitude z, the pressure is P(z) = P0*exp(z/H). We invert this to get the altitude, in meters, from the pressure at each
    # atmosphere layer. The surface pressure P0 is the pressure of the atmosphere at level 0, so we snag that and then compute.
    P0 = structure[0][0]
    atmoLayerAlt = np.array([ ( -SclHgt ) * math.log(Plevel[0]/P0) for Plevel in structure ])


    # Next, we use this altitude information to compute the mass density of gas within each atmosphere layer.
    # 
    # The linear thickness (in meters) of a layer is the distance from the bottom of the layer to the bottom of the next layer. Compute this for
    # each layer.
    atmoLayerThick = np.array([ atmoLayerAlt[i+1] - atmoLayerAlt[i] for i in range(len(atmoLayerAlt)-1) ])
    # The volumetric mass density (kg/m^3) of an ideal gas is equal to P * (mu/RT); compute for each layer.
    atmoLayerDensity = np.array([ (mmw*Plevel[0])/(cons.R*T)*10**5 for Plevel in structure ])
    # We also need column mass density (in kg/m^2) for each layer. We simply multiply the volume mass density by the linear thickness of each layer.
    Aseq = atmoLayerDensity[:-1] * atmoLayerThick
    # Finished with mass-structure computation.
    

    # Beginning cloud structure computation. First, determine the relative structure of the aerosol layers.
    #
    # Starting with the bottom layer, each subsequent layer has less aerosol matter. The reduction from one layer to the next is controlled by the
    # 'fade' parameter passed into this function. With, say, fade=0.1, each subsequent layer has 10% of the aerosol of the layer below it; or, for
    # layer 'i', there will be 0.1**i times as much aerosol as there is in layer 0. This is summarized in the PowerSet:
    PowerSet = np.array([fade**i for i in range(len(structure)-1)])
    PowerSet[-1] = 0.0 # We want the last (top) layer to go to 0 for computational reasons.

    # Next, compute what the mass mixing ratio must be at the bottom layer to pruduce the correct total column mass.
    #
    # By multipling the gas column densities by the power set, and adding them up, we get an "effective" atmosphere mass.
    Asum = sum(PowerSet*Aseq)
    # Convert total mass from log space (log10(kg/m2)) to linear space (kg/m2)
    totalMass = 10**totalMass
    # Finally, simply dividing the total mass by the effective atmosphere mass gives us the mass mixing ratio for the bottom layer.
    baseDensity = totalMass/(Asum)


    # With the base layer mass mixing ratio computed, we now simply multiply that value by the earlier PowerSet to get the mass mixing ratio at
    # all layers in the atmosphere.
    massStruc = PowerSet * baseDensity


    # We are done! Below, we manually compute the total column mass of clouds, to verify that it is correct.
    CMassSum = sum(massStruc * atmoLayerDensity[:-1] * atmoLayerThick)
    
    # Now we compare against the desired mass. Tolerance currently set at 1%.
    if abs(CMassSum/totalMass - 1) > 0.01:
        print('Warning: Cloud mass discrepency detected.')
        print('Found ', CMassSum, ' instead of ', totalMass)
        exit() # We quit out if the check fails.

    # Everything is good. We now append our cloud mixing ratios onto the original atmosphere pressure structure and return it.
    for i in range(len(structure)-1):
        structure[i].append(massStruc[i])
    structure[-1].append(0.0) # Make sure the top layer is 0 for computation reasons.

    return structure


def BuildClouds(Pc,totalMass,structure,mmw,fade,T):
    '''
    This function is handles the creation of the more complex, and generally thicker, "cloud bank" aerosol distribution. It shares some processes
    in common with BuildHaze.
    The clouds are computed such that the column mass of cloud matter in the atmosphere sums to a predefined (totalMass) value.
    Since PSG treats aerosols in a mass mixing ratio (kg/kg) we need to understand the mass structure of the atmosphere to correctly construct the
    haze.
    '''

    # Begin mass-structure computation. First, a bit of bookkeeping.
    #
    # Create a grid from the atmosphere grid, including temperature at each layer.
    structure = [ [i, T] for i in structure]


    # Second, we initialize the cloud structure. These values will be modified later.
    #
    # Placeholder value. This will eventually tell us the layer number of the cloud top.
    baseIndex = 0
    # Declare the thickness, in (initial) pressure layers, of the cloud deck.
    cloudCoreLayers = 3
    # We look for the correct base index using this loop.
    for i in range(len(structure)):
        # If the pressure too high, go to the next level.
        if structure[i][0] > Pc:
            continue
        # If the pressure isn't too high, we just passed the correct level. (Very unlikely that it is exact.) In this case, the previous layer
        # is the new base index.
        else:
            baseIndex = i-1
            break # Exit the loop immediately.
    # A catch to ensure we aren't putting clouds underground if we place the cloud deck at the chosen layer.
    if baseIndex-cloudCoreLayers < 0:
        raise BadClouds

    # Initialize the "power set" that determines the basic structure of clouds.This will contain multipliers for the cloud top level cloud
    # mass mixing ratio, to define each layer in relation to that layer.
    #
    # Create an empty structure, to loop through.
    PowerSet = np.zeros(len(structure)-1)
    # Many if statements in this loop, as there are numerous 'regions' in this structure..
    for i in range(len(PowerSet)):

        # At very low altitudes, there are no clouds.
        if i < baseIndex-4-cloudCoreLayers-1:
            PowerSet[i] = 0

        # Just above the "very low altitudes", the clouds are very thin.
        if i == baseIndex-4-cloudCoreLayers: 
            PowerSet[i] = 10**-10

        # Clouds fade in as we approach the bottom of the cloud deck.
        if i >= (baseIndex-4-cloudCoreLayers) and i < baseIndex-cloudCoreLayers:
            PowerSet[i] = fade**(baseIndex-i)

        # Within the core of the cloud, full density clouds.
        if i >= baseIndex-cloudCoreLayers and i <= (baseIndex):
            PowerSet[i] = 10**0

        # Above the top of the clouds, the clouds fade away again for some layers.
        if i > (baseIndex) and PowerSet[i-1] > 10**-7: 
            PowerSet[i] = fade**(i-baseIndex)

        # Beyond a certain point, we start checking for the clouds to vanish.
        if i > baseIndex and PowerSet[i-1] <= 10**-7:
            # Once the clouds are thin enough, they vanish.
            if PowerSet[i-1] <= 10**-9:
                PowerSet[i] = 0

            # Otherwise, make them thin enough. The transition should be subtle enough to not cause issues.
            else:
                PowerSet[i] = 10**-9


    # All basic computation is complete. We need to do some extra bookkeeping now, before making adjustments for optical depth.
    #
    # The first parameters and calculations are done to ensure our layers are not too optically thick.
    # Cloud particle effective cross section. Note that this assumes a grey cloud. This is used to calculate optical depth.
    EXcld = 150 
    # The number of angles PSG is set to compute (parameter <GEOMETRY-DISK-ANGLES>).
    DAngles = 6.0
    # Calculate the 'critical' optical depth. If a layer has an optical depth greater than this value, PSG will have issues computing radiative
    # transfer through the layer at high angles of incidence (close to the limb). Influenced by DAngles.
    tauCrit = math.cos( (np.pi/2.0)*(2*DAngles-1.0)/(2*DAngles) ) 
    # Now, some more general bookkeeping.

    # Reassign the pressure level at the baseIndex to be exactly equal to the cloud top pressure.
    structure[baseIndex+1:][0] += structure[baseIndex+1][0] - Pc

    # Placeholder pressure structure for manipulation in the subsequent loop.
    Pstruct = [ i[0] for i in structure ]

    # Extract surface pressure for computational purposes.
    P0 = structure[0][0]
    # Move total mass out of log space kg/m2
    totalMass = 10**totalMass
    # Compute one scale height, in meters.
    SclHgt = (cons.R*T)/(mmw*cons.g)
    # Right now, we do not know if our clouds are reasonable. 
    ReasonableClouds = False


    # Bookeeping and initialization done. The following loop makes slight adjustments to the layering until there are no optical depth conflicts,
    # as represented by the 'Reasonable Clouds' variable.
    while ReasonableClouds == False:

        # First, we create a layer list in distance-space, rather than in pressure-space.
        #
        # For any given altitude z, the pressure is P(z) = P0*exp(z/H). We invert this to get the altitude, in meters, from the pressure at each
        # atmosphere layer. The surface pressure P0 is the pressure of the atmosphere at level 0, so we snag that and then compute.
        atmoLayerAlt = np.array([ ( -SclHgt )*math.log(Plevel[0]/P0) for Plevel in structure ])

        # Next, we use this altitude information to compute the mass density of gas within each atmosphere layer.
        # 
        # The linear thickness (in meters) of a layer is the distance from the bottom of the layer to the bottom of the next layer. Compute this for
        # each layer.
        # heightStructure = [ atmoLayerAlt[i+1] - atmoLayerAlt[i] for i in range(len(structure)-1) ] Possibly obsolete variable.
        atmoLayerThick = np.array([ atmoLayerAlt[i+1] - atmoLayerAlt[i] for i in range(len(atmoLayerAlt)-1) ] )
        # The volumetric mass density (kg/m^3) of an ideal gas is equal to P * (mu/RT); compute for each layer.
        atmoLayerDensity = np.array([ (mmw*Plevel[0])/(cons.R*T)*10**5 for Plevel in structure ])
        # We need column mass density (in kg/m^2) for each layer. We simply multiply the volume mass density by the linear thickness of each layer.
        Aseq = atmoLayerDensity[:-1] * atmoLayerThick


        # Next, compute what the mass mixing ratio must be at the bottom layer to pruduce the correct total column mass.
        #
        # By multipling the gas column densities by the power set, and adding them up, we get an "effective" atmosphere mass.
        Asum = sum(PowerSet*Aseq)
        # Finally, simply dividing the total mass by the effective atmosphere mass gives us the mass mixing ratio for the bottom layer.
        baseDensity = totalMass/(Asum)


        # Here things begin to diverge from BuildHaze(). We need to assess the optical depth situation.
        #
        # Combining baseDensity with PowerSet we have a first-pass cloud structure for the atmosphere; multiply by atmosphere density to find
        # volumetric mass density of the aerosol; multiply by the linear thickness of the atmosphere layer to find column mass density.
        baseColumn = (baseDensity * PowerSet) * (atmoLayerThick * atmoLayerDensity[:-1])
        # Compute the optical depth of each atmosphere layer by multiplying by the aerosol effective cross section.
        taucore = EXcld * baseColumn


        # These boolean flags will identify what "regions" need modification, handled later.
        coreGood = True
        topGood = True
        bottomGood = True
        # Build a framework for tracking layers which are too optically thick.
        badLayerList = []
        # Examine each layer's optical depth, flagging any that are too thick.
        for i in range(len(taucore)):
            if taucore[i] > tauCrit:
                # The clouds have been flagged as too thick here. Add this index to the list.
                badLayerList.append(i)
                # Identify which region this layer belongs to, and flag it as "bad".
                if i < baseIndex: # Bottom section.
                    bottomGood = False
                if i >= baseIndex - cloudCoreLayers and i <= baseIndex: # Core section.
                    coreGood = False
                if i > baseIndex: # Top section.
                    topGood = False

        # Check if we're done and exit the loop if we are.
        if len(badLayerList) == 0:
            ReasonableClouds = True # All good, clouds are reasonable.
            break
        

        # If we made it this far, there must be layers that need to be adjusted.
        #
        # These loops adjust cloud layers. Any layer which is too optically thick, will be split into two layers. The new layer will have a pressure
        # level calculated halfway between the flagged layer and the next-lowest-pressure layer.


        # We handle the cloud core region in the first loop. The layers are run in reverse order in this loop, for indexing reasons.
        if not coreGood:
            for i in badLayerList[::-1]:
                # The extra if statement avoids working unintended layers.
                if i >= baseIndex - cloudCoreLayers and i <= baseIndex:
                    # This layer is in the core.

                    # Calculate the new pressure.
                    newP = np.sqrt(structure[i][0] * structure[i+1][0])
                    # newP = 10**((np.log10(structure[i][0])+np.log10(structure[i+1][0]))/2.)

                    # Insert our new pressure level into the atmosphere structure.
                    structure.insert(i+1,[newP,T]) 

                    # The aerosol density of the new level will be calculated similarly to the new pressure; here, we calculate that density.
                    # Note that in the core, the density at PowerSet[i] and PowerSet[i+1] will often be the same.
                    newPwr = np.sqrt(PowerSet[i] * PowerSet[i+1])

                    # Insert the new power number into power set for next loop.
                    PowerSet = np.insert(PowerSet,i+1,newPwr)

                    # The core of the cloud now covers one addional layer; account for this.
                    # NOTE: Shouldn't the cloud top layer also change?
                    cloudCoreLayers += 1


        # We handle the below-the-cloud region in the second loop. The layers are run in reverse order in this loop, for indexing reasons.
        # To avoid indexing issues, we only address these layers if the core region is confirmed "good".
        if coreGood and not bottomGood:
            for i in badLayerList[::-1]:
                # The extra if statement avoids working unintended layers.
                if i < baseIndex - cloudCoreLayers:
                    # This layer is under the core.
                    # Calculate the new pressure.
                    newP = np.sqrt(structure[i][0] * structure[i+1][0])
                    # Insert our new pressure level into the atmosphere structure.
                    structure.insert(i+1,[newP,T])
                    # The aerosol density of the new level will be calculated similarly to the new pressure; here, we calculate that density.
                    newPwr = np.sqrt(PowerSet[i] * PowerSet[i+1])
                    # Insert the new power number into power set for next loop.
                    PowerSet = np.insert(PowerSet,i+1,newPwr)
                    # The index of the cloud top layer goes up by one when we insert a new layer below it.
                    baseIndex += 1
                    # We break because we moved some indexes of potentially-bad cloud layers. If we continued, we could have indexing issues.
                    break

        # We handle the above-the-cloud region in the third loop. We no longer need to run these in reverse order.
        # To avoid indexing issues, we only address these layers if the core and bottom regions are confirmed "good".
        if coreGood and bottomGood and not topGood:
            for i in badLayerList:
                # The extra if statement avoids working unintended layers.
                if i > baseIndex:
                    # This layer is under the core.
                    # Calculate the new pressure.
                    newP = np.sqrt(structure[i][0] * structure[i+1][0])
                    # Insert our new pressure level into the atmosphere structure.
                    structure.insert(i+1,[newP,T])
                    # The aerosol density of the new level will be calculated similarly to the new pressure; here, we calculate that density.
                    newPwr = np.sqrt(PowerSet[i] * PowerSet[i+1])
                    # Insert the new power number into power set for next loop.
                    PowerSet = np.insert(PowerSet,i+1,newPwr)
                    # We break because we moved some indexes of potentially-bad cloud layers. If we continued, we could have indexing issues.
                    break

        # PSG can only handle up to 500 pressure levels. Here we count the layers to make sure we are under this limit.
        if len(structure[0]) > 500:
            print('Warning: data overflow, too many layers')
            print('l=',totalMass,'P0=',P0,'Pc=',Pc)
            exit()

    # Our clouds are successfully defined in units of kg/kg; append them to the structure.
    massStruc = [ (PowerSet[i])*baseDensity for i in range(len(structure)-1) ]

    # We are done! Below, we manually compute the total column mass of clouds, to verify that it is correct.
    CMassSum = sum(massStruc * atmoLayerDensity[:-1] * atmoLayerThick)
    # Now we compare against the desired mass. Tolerance currently set at 1%.
    if abs(CMassSum/totalMass - 1) > 0.01:
        print('Warning: Cloud mass discrepency detected.')
        print('Found ', CMassSum, ' instead of ', totalMass)
        exit()

    # Everything is good. We now append our cloud mixing ratios onto the original atmosphere pressure structure and return it.
    for i in range(len(structure)):
        try:
            structure[i].append(massStruc[i]) # If clouds exist in the layer, put them in.
        except:
            structure[i].append(0.0) # Otherwise, set to 0.

    return structure


def CallPSG(CFGname):
    '''
    This function actually calls the Docker PSG via curl command, giving it the input file name (payload) and target file name (destination).
    CFGname is the name of the configuration file to be run through PSG.
    May need adjustment if you play with the 'Clear' 'Cloudy' 'Hazy' file type structure, or the CFG locations.
    '''

    # Determine the destination filename.
    if 'ct' in CFGname: # Is this a cloudy model?
        dataName = './Data/Clouds/'+CFGname[0:-4]+'.dat' # Send it to the cloudy folder.
    elif 'r' in CFGname and 'l' in CFGname: # It's not a cloudy model, but does it have cloud particulate at all?
        dataName = './Data/Haze/'+CFGname[0:-4]+'.dat' # Must be hazy, send it to the hazy folder.
    else: # It's not cloudy, and it's not hazy, so it must be clear.
        dataName = './Data/NoClouds/'+CFGname[0:-4]+'.dat' # Send it to the clearsky folder.

    # Name and location of the source configuration file.
    inputName = './Cfgs/'+CFGname

    # Feed the source file and destination path into a curl PSG call.
    Command = 'curl --data-urlencode file@'+inputName+' http://localhost:3000/api.php > '+dataName
    os.system(Command)

    return dataName
# -*- coding: utf-8 -*-
# encoding: utf-8

from yade import utils,ymport,export,plot
import numpy as np
import os, sys, importlib, shutil

# optional use of table in batch-mode
#utils.readParamsFromTable(PACK=None,LITH=None,ANIINCL=None,KINEB=None,noTableOk=True) # batch mode: define parameters from table
#from yade.params.table import *

# no batch-mode use : yadedaily --cores '4' loadTile.py  |& tee terminalOutput_loadTile.txt
PACK='tile_75-4-22_0.02_89212'  
LITH='isoLith19'
ANIINCL='iso'
KINEB=2e10

print('lithology (LITH): ', LITH)
print('inclination of the weakness plane (ANIINCL): ', ANIINCL)


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##### basic setup

#### packing used
packing=PACK                                                            # don't forget to adjust the total volume calculation (totVol) when the input geometry is changing !!! 
print('packing used: ', packing)


#### Scale Factor
scaleFactor=1000.                                                       # set the hight of the slope after cutting the tile in next step, in [m]


#### stability criterion
kinE_bound=KINEB                                                        # kinetic energy bound as stabilisation criterion [J] (after gravity and weakness plane introduction)
#kinE_bound=2e10 
#kinE_bound=1e16
print('kinetic energy bound:',"{:.1e}".format(kinE_bound),'[J]')


#### gravity increase step                                                        # gravity increase step [m/s²] for applying gravity, 9.81 sudden increase to final value
gStep=9.81
#gStep=GSTEP
print('gravity increase step (gStep)', gStep)


#### predefining moving average characteristics (applied on kinetic energy and used for stabilisation function)
movAveInt=1000
movAveSamplingNb=10                                                 # arbitrarily set to 10
movAveSamplingStep=int(movAveInt/movAveSamplingNb)                  # do not change that

print('movAveInt=',movAveInt)
print('movAveSamplingNb=',movAveSamplingNb)
print('movAveSamplingStep=',movAveSamplingStep)


#### interparticle properties, load lithology
lithology=LITH
path=os.getcwd()                                                        # extract path
sys.path.append(path)
filename='%s' % lithology
lith_module = importlib.import_module(filename)
module_dict = lith_module.__dict__
try:
    to_import = lith_module.__all__
except AttributeError:
    to_import = [name for name in module_dict if not name.startswith('_')]
globals().update({name: module_dict[name] for name in to_import})

particleDensity=4000                                                    # preliminary value, will be altered below.


#### Weakness plane orientation: RELEVANT ONLY IF ANI=True -> SET smoothJoint=True in Law2_ScGeom_JCFpmPhys_JointedCohesiveFrictionalPM
if ANI==True:
    GAMMA=ANIINCL               # angle of weakness plane with respect to the horizontal direction (0 if planes are horizontal, 90 if planes are vertical  - X direction)


#### assign output file name here
if ANI==True :
    output=packing+'_'+lithology+'_'+str(ANIINCL)+'_loaded' #+'_IncFac'+str(strengthIncreaseFactor)   #+'_gStep'+str(gStep)
if ANI==False :
    output=packing+'_'+lithology+'_noANI'+'_loaded' #+'_IncFac'+str(strengthIncreaseFactor)   #+'_gStep'+str(gStep)


#### creating a folder dedicated to a single run for output files
mainDir=output                                  # create main directory
mainPath=os.path.join(path,mainDir)
if os.path.exists(mainPath) == True:
    shutil.rmtree(mainPath)
    os.mkdir(mainPath)
else:
    os.mkdir(mainPath)

posDir='positions'                              # create positions directory
posPath=os.path.join(mainPath,posDir)
if os.path.exists(posPath) == True:
    shutil.rmtree(posPath)
    os.mkdir(posPath)
else:
    os.mkdir(posPath)


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##### calculate dimensions and calculate the correct density for the particles

#### import sphere packing for pre-processing
def sphereMat0(): return JCFpmMat(type=1,density=particleDensity,young=YOUNG,poisson=ALPHA,tensileStrength=TENS,cohesion=COH,frictionAngle=radians(FRICT))
O.bodies.append(ymport.text(packing+'.spheres',scale=scaleFactor,shift=Vector3(0,0,0),material=sphereMat0))


#### get dimensions of the packing
def dimensionsPack():
    global dim, xinf, xsup, X, yinf, ysup, Y, zinf, zsup, Z
    dim=utils.aabbExtrema()
    xinf=dim[0][0]
    xsup=dim[1][0]
    X=xsup-xinf
    yinf=dim[0][1]
    ysup=dim[1][1]
    Y=ysup-yinf
    zinf=dim[0][2]
    zsup=dim[1][2]
    Z=zsup-zinf

dimensionsPack()
print('xinf=',xinf,' | yinf=',yinf,' | zinf=',zinf)
print('xsup=',xsup,' | ysup=',ysup,' | zsup=',zsup)
print('X=',X,' | Y=',Y,' | Z=',Z)


#### get dimensions of spheres
def dimensionsSph():
    global R, nbSpheres, Rmax, Rmin, Rmean
    R=0
    Rmax=0
    Rmin=1e6
    nbSpheres=0
    Rmean=0
    for o in O.bodies:
        if isinstance(o.shape,Sphere):
            nbSpheres+=1
            R+=o.shape.radius
            if o.shape.radius>Rmax:
                Rmax=o.shape.radius
            if o.shape.radius<Rmin:
                Rmin=o.shape.radius
    Rmean=R/nbSpheres

dimensionsSph()
print('nbSpheres=',nbSpheres)
print('Rmax=',Rmax,' | Rmin=',Rmin,' | Rmean=',Rmean)


#### compute volume function
def volume():
    global packingVolume, volSpheres
    packingVolume=0
    volSpheres=0

    packingVolume=X*Y*Z
    for o in O.bodies:
        if isinstance(o.shape,Sphere):
            volSpheres+=(4./3.)*pi*(o.shape.radius)**3.

volume()
print('packingVolume=',packingVolume)
print('volSpheres=',volSpheres)


#### calculate compacity
comp=volSpheres/packingVolume
print('initial compacity=',comp)


#### particle density calculation
particleDensity=rho_init*(packingVolume/volSpheres)
print('particleDensity=',particleDensity)


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##### set material properties and boundary conditions

#### import things in reset simulation
O.reset()

### particles
def sphereMat(): return JCFpmMat(type=1,density=particleDensity,young=YOUNG,poisson=ALPHA,tensileStrength=TENS,cohesion=COH,frictionAngle=radians(FRICT),jointNormalStiffness=YOUNG/(pi*Rmean),jointShearStiffness=ALPHA*YOUNG/(pi*Rmean),jointTensileStrength=TENS,jointCohesion=COH,jointFrictionAngle=radians(FRICT),jointDilationAngle=radians(0))
#jointShearStiffness=Young/(pi*Rmean) enables to define similar stiffness for structures and matrix as default value
O.bodies.append(ymport.text(packing+'.spheres',scale=scaleFactor,shift=Vector3(0,0,0),material=sphereMat))

### add wall to bottom
def wallMat(): return JCFpmMat(type=0,density=particleDensity,young=0.01*YOUNG,poisson=ALPHA,tensileStrength=0,cohesion=0,frictionAngle=radians(FRICT))
O.bodies.append(wall((xinf+X/2.,yinf+Y/2.,zinf),2,sense=0,color=None,material=wallMat))

### add wall to all 4 vertical sides (to avoid particles "jumping" over the boundaries
O.bodies.append(wall((xinf,yinf+Y/2.,zinf+Z/2.),0,sense=0,color=None,material=wallMat))
O.bodies.append(wall((xsup,yinf+Y/2.,zinf+Z/2.),0,sense=0,color=None,material=wallMat))
O.bodies.append(wall((xinf+X/2.,yinf,zinf+Z/2.),1,sense=0,color=None,material=wallMat))
O.bodies.append(wall((xinf+X/2.,ysup,zinf+Z/2.),1,sense=0,color=None,material=wallMat))

#### define boundary conditions
e=3*Rmean
baseBodies=[]
for o in O.bodies:
    if isinstance(o.shape,Sphere):
        
        ### these ensure the boundaries to not move in their normal direction ("roller" like BC)
        ## front particles
        if o.state.pos[0]<(xinf+e) :
            o.state.blockedDOFs='x'
            o.shape.color=(1,0,0)
        ### back particles
        if o.state.pos[0]>(xsup-e) :
            o.state.blockedDOFs='x'
            o.shape.color=(1,0,0)
        ## left particles
        if o.state.pos[1]<(yinf+e) :
            o.state.blockedDOFs+='y'
            o.shape.color=(0,0,1)
        ## right particles
        if o.state.pos[1]>(ysup-e) :
            o.state.blockedDOFs+='y'
            o.shape.color=(0,0,1)
        ## ground particles
        if o.state.pos[2]<(zinf+e) :
            o.state.blockedDOFs='z'
            o.shape.color=(0,0,0)
            baseBodies.append(o.id)


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##### other stuff

#### Add here colouring of spheres, for better visualising weakness plane introduction

## assign stripe colours
stripeColour1=(0,0,0)
stripeColour2=(1,1,1)

#list of colours:       
                        #Black:    (0,0,0)
                        #White:    (1,1,1)
                        #Red:      (1,0,0)
                        #Green:    (0,1,0)
                        #Blue:     (0,0,1)
                        #Yellow:   (1,1,0)
                        #Cyan:     (0,1,1)
                        #Magenta:  (1,0,1)


## insertion of stripe colouring properties and loop that colours particles
if ANI==True:
    numStripe=2                                         # number of colour stripes per colour interval
    numColourInt=75                                     # number of colour intervals
    colourIntThickness=Z/22                             # colour interval thickness
    stripeThickness=colourIntThickness/numStripe        # (do not change)
    
    for o in O.bodies:                                  
        if isinstance(o.shape,Sphere):
            for i in range(numStripe):
                for j in range(numColourInt):
                    if GAMMA>=0 and GAMMA<90 :
                        if i==0 and o.state.pos[2]>(((i*stripeThickness/cos(radians(GAMMA)))+(tan(radians(GAMMA))*o.state.pos[0]*-1))+(((numColourInt/2)-1-j)*colourIntThickness/cos(radians(GAMMA)))) and o.state.pos[2]<=(((i+1)*stripeThickness/cos(radians(GAMMA)))+(tan(radians(GAMMA))*o.state.pos[0]*-1)+(((numColourInt/2)-1-j)*colourIntThickness/cos(radians(GAMMA)))):
                            o.shape.color=stripeColour1 # colour can be changed here
                        if i==1 and o.state.pos[2]>(((i*stripeThickness/cos(radians(GAMMA)))+(tan(radians(GAMMA))*o.state.pos[0]*-1))+(((numColourInt/2)-1-j)*colourIntThickness/cos(radians(GAMMA)))) and o.state.pos[2]<=(((i+1)*stripeThickness/cos(radians(GAMMA)))+(tan(radians(GAMMA))*o.state.pos[0]*-1)+(((numColourInt/2)-1-j)*colourIntThickness/cos(radians(GAMMA)))):
                            o.shape.color=stripeColour2 # colour can be changed here
                    if GAMMA==90 :
                        if i==0 and o.state.pos[0]>(i*stripeThickness)+(j-(numColourInt/2))*colourIntThickness and o.state.pos[0]<=(i+1)*stripeThickness+(j-(numColourInt/2))*colourIntThickness:
                            o.shape.color=stripeColour1 # colour can be changed here
                        if i==1 and o.state.pos[0]>(i*stripeThickness)+(j-(numColourInt/2))*colourIntThickness and o.state.pos[0]<=(i+1)*stripeThickness+(j-(numColourInt/2))*colourIntThickness:
                            o.shape.color=stripeColour2 # colour can be changed here
                    if GAMMA>=91 and GAMMA<180 :
                        if i==0 and o.state.pos[2]>(((i*stripeThickness/cos(radians(180-GAMMA)))+(tan(radians(180-GAMMA))*o.state.pos[0]*1))+(((numColourInt/2)-1-j)*colourIntThickness/cos(radians(180-GAMMA)))) and o.state.pos[2]<=(((i+1)*stripeThickness/cos(radians(180-GAMMA)))+(tan(radians(180-GAMMA))*o.state.pos[0]*1)+(((numColourInt/2)-1-j)*colourIntThickness/cos(radians(180-GAMMA)))):
                            o.shape.color=stripeColour1 # colour can be changed here
                        if i==1 and o.state.pos[2]>(((i*stripeThickness/cos(radians(180-GAMMA)))+(tan(radians(180-GAMMA))*o.state.pos[0]*1))+(((numColourInt/2)-1-j)*colourIntThickness/cos(radians(180-GAMMA)))) and o.state.pos[2]<=(((i+1)*stripeThickness/cos(radians(180-GAMMA)))+(tan(radians(180-GAMMA))*o.state.pos[0]*1)+(((numColourInt/2)-1-j)*colourIntThickness/cos(radians(180-GAMMA)))):
                            o.shape.color=stripeColour2 # colour can be changed here


## include a checkerboard colouring scheme for the isotropic case
if ANI==False:
    squareSize=200 # choose the length and height of the checkerboard squares
    
    x_lower = (xinf//squareSize)*squareSize
    x_upper = ((xsup//squareSize)*squareSize)+squareSize
    
    cb_nbOfIntervX = int(abs(x_lower/squareSize)+abs(x_upper/squareSize))
    
    z_lower = (zinf//squareSize)*squareSize
    z_upper = ((zsup//squareSize)*squareSize)+squareSize
    
    cb_nbOfIntervZ = int(abs(z_lower/squareSize)+abs(z_upper/squareSize))
    
    cb_matrix = np.zeros((cb_nbOfIntervX,cb_nbOfIntervZ),dtype=int) # creation of a matrix with zeros and ones like a checkerboard
    cb_matrix[1::2,::2] = 1
    cb_matrix[::2,1::2] = 1

    for o in O.bodies:                                  
        if isinstance(o.shape,Sphere):
            for i in range(len(cb_matrix)):
                for j in range(len(cb_matrix[i])):
                    if o.state.pos[0] > x_lower+i*squareSize and o.state.pos[0] <= x_lower+(i+1)*squareSize and o.state.pos[2] > z_lower+j*squareSize and o.state.pos[2] <= z_lower+(j+1)*squareSize and cb_matrix[i][j] == 0 :
                        o.shape.color=stripeColour1
                    if o.state.pos[0] > x_lower+i*squareSize and o.state.pos[0] <= x_lower+(i+1)*squareSize and o.state.pos[2] > z_lower+j*squareSize and o.state.pos[2] <= z_lower+(j+1)*squareSize and cb_matrix[i][j] == 1 :
                        o.shape.color=stripeColour2


## colouring particles uniformly
#if ANI==False:
    #for o in O.bodies:
        #if isinstance(o.shape,Sphere):
            #o.shape.color=(1,0,1) # colour can be changed here


#### Identify indicator on top of tile
refPoint=0
Xref=xinf+X/2.
Yref=yinf+Y/2.
Zref=zsup-2*Rmax
for o in O.bodies:
    if isinstance(o.shape,Sphere):
        if o.state.pos[0]>(Xref-Rmax) and o.state.pos[0]<(Xref+Rmax) and o.state.pos[1]>(Yref-Rmax) and o.state.pos[1]<(Yref+Rmax) and o.state.pos[2]>(Zref-Rmax) and o.state.pos[2]<(Zref+Rmax) :
            refPoint=o.id
            px0=o.state.pos[0]
            py0=o.state.pos[1]
            pz0=o.state.pos[2]

print('refPoint=',refPoint,' | Xref=',Xref,',Yref=',Yref,', Zref=',Zref)
O.bodies[refPoint].shape.color=(1,0,0)

### compute the vertical stress
sigmaZbase=0
def vertStress():
    global sigmaZbase
    sigmaZbase=utils.sumForces(baseBodies,(0,0,1))/(X*Y)


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##### function for moving average of kinetic energy

#### creating numpy array for moving average of kinetic energy
movAve_kinE=0.
dynList_kinE=np.empty(movAveSamplingNb)
dynList_kinE.fill(kinE_bound*1.1)

#### the function for the moving average of kinetic energy
def movAverage_kinEnergy():
    global dynList_kinE, movAve_kinE
    kinE=utils.kineticEnergy()
    
    dynList_kinE=np.append(dynList_kinE,kinE)
    dynList_kinE=np.delete(dynList_kinE,0)
    movAve_kinE=np.mean(dynList_kinE)


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##### define the engines and corresponding functions

#### define engines
O.engines=[
    ForceResetter(),
    InsertionSortCollider([Bo1_Sphere_Aabb(aabbEnlargeFactor=intR,label='aabb'),Bo1_Wall_Aabb()]),
    InteractionLoop(
        [Ig2_Sphere_Sphere_ScGeom(interactionDetectionFactor=intR,label='isssg'),Ig2_Wall_Sphere_ScGeom()],
        [Ip2_JCFpmMat_JCFpmMat_JCFpmPhys(cohesiveTresholdIteration=1,label='interactionPhys')],
        [Law2_ScGeom_JCFpmPhys_JointedCohesiveFrictionalPM(recordCracks=True,smoothJoint=ANI,Key=(output),label='interactionLaw')]),
        # set here smoothJoint=True for introducing anisotropy/weakness plane

    GlobalStiffnessTimeStepper(defaultDt=0.1*utils.PWaveTimeStep(),timestepSafetyCoefficient=0.3), # timestepSafetyCoefficient=0.5 used for earlier simulations
    NewtonIntegrator(damping=0.5,gravity=(0.,0.,0.),label='newton',dead=1),

    #recording functions
    PyRunner(iterPeriod=50,initRun=True,command='myRecorder()',label='recData'),
    PyRunner(iterPeriod=1,initRun=True,command='positions()',label='recPositions'),

    #VTK recorder
    VTKRecorder(iterPeriod=1,initRun=True,fileName=(output+'/'+output+'_'),Key=(output),recorders=['spheres','velocity','colors','cracks','jcfpm','bstresses'],label='recVTK'),

    #moving average pf kinetic energy function
    PyRunner(command='movAverage_kinEnergy()',label='movAveKinE')]


#### Record data
V=0
def myRecorder():
    global V
    for o in O.bodies :                             # |
        if abs(o.state.vel[2])>0.01:                # | It calculates the volume that moves in z-direction
            V+=1.333*pi*pow(o.shape.radius,3)       # |
    rId=O.bodies[refPoint]
    plot.addData(Ek=utils.kineticEnergy(),v=rId.state.vel[2],g=newton.gravity[2],px=rId.state.pos[0]-px0,py=rId.state.pos[1]-py0,pz=rId.state.pos[2]-pz0,iterations=O.iter,t=O.time,vol=V,tc=interactionLaw.nbTensCracks,sc=interactionLaw.nbShearCracks,unbF=utils.unbalancedForce())
    plot.saveDataTxt(output+'/'+output)
    V=0                             #reset Volume calculation

def positions():
    export.text(output+'/'+posDir+'/'+output+'-pos_'+str(O.iter))


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##### set up model

#### recorders off
recPositions.dead=1 #|  -> first, set all recorders off
recVTK.dead=1       #|

#### set sampling step for moving average of kinetic energy
movAveKinE.iterPeriod=movAveSamplingStep

#### create interactions between particles
recVTK.dead=0
O.step()
recVTK.dead=1


#### coordination number verification
numCohesivelinks=0
for i in O.interactions:
    if not i.isReal : continue
    if isinstance(O.bodies[i.id1].shape,Sphere) and isinstance(O.bodies[i.id2].shape,Sphere) and i.phys.isCohesive :
      numCohesivelinks+=1

print('total nb of bonds=',numCohesivelinks)
print('coordination number =',2.0*numCohesivelinks/nbSpheres)


#### calculate vertical stress
vertStress()
print('vertStress = ', sigmaZbase)


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##### stabilisation function

#### defining the stabilisation loop
printInt=1          #printing interval for stabilisation loop, printInt=1 prints out every sampling step

def stabilisation():
    count=printInt
    while 1:
        O.run(movAveSamplingStep,True)           # the moving average sampling step (movAveSamplingStep) is used to define the interval for the stability check

        if movAve_kinE>=kinE_bound:
            count-=1
            if count==0:
                print('iter= ',O.iter,', stabilisation, movAve_kinE=', round(movAve_kinE))
                count=printInt
        elif movAve_kinE<kinE_bound:
            print('iter= ',O.iter,', stability criterion reached, movAve_kinE=', round(movAve_kinE))
            print('slope is stable!')
            break


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##### Introduce weakness plane in the packing -> SET smoothJoint=TRUE in engines->Law2_ScGeom_JCFpmPhys_JointedCohesiveFrictionalPM()

#### the loop for weakness plane introduction
if ANI==True:
    print('IDENTIFYING WEAKNESS PLANE BONDS')
    nbWPB=0
    BEDNORM=Vector3(sin(radians(GAMMA)),0,cos(radians(GAMMA))) # orientation of normal to weakness plane -> CAREFUL WITH THE MODEL ORIENTATION
    print('angle of weakness plane with respect to horizontal =', GAMMA, ' | normal =', BEDNORM)
    # search for weakness plane contacts and modification of properties
    for i in O.interactions:
        if isinstance(O.bodies[i.id1].shape,Sphere) and isinstance(O.bodies[i.id2].shape,Sphere): # particles only
            pdct=(i.geom.normal).dot(BEDNORM)
            if abs(pdct) > cos(radians(DGAMMA)): # condition that ensures the contact normal belongs to the cone defined by the weakness plane's normal with more or less dGamma degrees
                ### contact reorientation -> use of smooth contact model
                # bodies
                O.bodies[i.id1].state.onJoint=True
                O.bodies[i.id2].state.onJoint=True
                O.bodies[i.id1].state.joint=1
                O.bodies[i.id2].state.joint=1
                O.bodies[i.id1].state.jointNormal1=BEDNORM
                O.bodies[i.id2].state.jointNormal1=BEDNORM
                
                ## interactions
                #if STENS==0 and SCOH==0:
                    #i.phys.isCohesive=False
                i.phys.isOnJoint=True
                i.phys.jointNormal=BEDNORM
                i.phys.jointNormal=BEDNORM*np.sign(i.phys.jointNormal.dot(i.geom.normal))
                i.phys.initD = abs((O.bodies[i.id1].state.pos - O.bodies[i.id2].state.pos).dot(i.phys.jointNormal))
                
                
                nbWPB+=1


    print('number of weakness plane bonds:',nbWPB)


#### setting the weakness plane's stiffness
if ANI==True:
    print('SETTING THE WEAKNESS PLANE STIFFNESS')
    for i in O.interactions :
        if i.phys.isOnJoint==True :
            # bodies
            O.bodies[i.id1].mat.jointNormalStiffness=SNSTIFF*i.phys.kn/i.phys.crossSection
            O.bodies[i.id2].mat.jointNormalStiffness=SNSTIFF*i.phys.kn/i.phys.crossSection
            O.bodies[i.id1].mat.jointShearStiffness=SSSTIFF*i.phys.ks/i.phys.crossSection
            O.bodies[i.id2].mat.jointShearStiffness=SSSTIFF*i.phys.ks/i.phys.crossSection
            # interactions
            i.phys.kn*=SNSTIFF
            i.phys.ks*=SSSTIFF

    ## VTK recording

    recVTK.dead=0
    O.step()
    recVTK.dead=1


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##### gravity loading

#### initialize stresses in slope (gravity)
print('STRESS INITIALIZATION  (apply gravity)')

newton.dead=0
gravityZ=0

#gStep=0.1

while 1:
    if gravityZ>=-9.81+gStep:
        gravityZ-=gStep
        newton.gravity=(0.,0.,gravityZ)
        print('iter= ',O.iter,', gravity increase, g=', gravityZ)
    elif gravityZ<=-9.81+gStep and gravityZ!=-9.81:
        gravityZ=-9.81
        newton.gravity=(0.,0.,gravityZ)
        print('iter= ',O.iter,', last gravity increase, g=', gravityZ)
    elif gravityZ==-9.81:
        print('iter= ',O.iter,', STRESS INITIALIZATION (apply gravity) completed')
        break

    ## stabilisation and a VTK recording
    stabilisation()

    vertStress()
    print('vertStress = ', sigmaZbase)

    recVTK.dead=0
    O.step()
    recVTK.dead=1


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##### save and stop simulation

#### recalculate density of packing after gravitational settling (and weakness plane introduction)
dimensionsPack()
print('final dimensions of packing:')
print('xinf=',xinf,' | yinf=',yinf,' | zinf=',zinf)
print('xsup=',xsup,' | ysup=',ysup,' | zsup=',zsup)
print('X=',X,' | Y=',Y,' | Z=',Z)
volume()
print('final packingVolume=',packingVolume)
comp=volSpheres/packingVolume
print('final compacity=',comp)
rho_settled=particleDensity*comp
print('final density of loaded tile =',rho_settled)


#### save simulation
#O.save(output+'/'+output+'.yade')
O.save(output+'.yade')


#### stop simulation
O.pause()


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##### move all output vtk-files in dedicated folders

#### create folders for different vtk-files
cracksDir='cracks'                              # create cracks directory
cracksPath=os.path.join(mainPath,cracksDir)
if os.path.exists(cracksPath) == True:
    shutil.rmtree(cracksPath)
    os.mkdir(cracksPath)
else:
    os.mkdir(cracksPath)

intrsDir='interactions'                         # create interactions directory
intrsPath=os.path.join(mainPath,intrsDir)
if os.path.exists(intrsPath) == True:
    shutil.rmtree(intrsPath)
    os.mkdir(intrsPath)
else:
    os.mkdir(intrsPath)

spheresDir='spheres'                            # create spheres directory
spheresPath=os.path.join(mainPath,spheresDir)
if os.path.exists(spheresPath) == True:
    shutil.rmtree(spheresPath)
    os.mkdir(spheresPath)
else:
    os.mkdir(spheresPath)


#### find files in main directory
os.chdir(mainPath) 
VTKfiles = os.listdir(mainPath)


#### moving files
for f in VTKfiles:
    if (f.startswith(output + '_cracks')) :
        shutil.move(os.path.abspath(f), cracksPath)
    elif (f.startswith(output + '_intrs')) :
        shutil.move(os.path.abspath(f), intrsPath)
    elif (f.startswith(output + '_spheres')) :
        shutil.move(os.path.abspath(f), spheresPath)

if ANI==True:
    shutil.move(path + '/cracks_'+output+'.txt', cracksPath)
    os.chdir(path)


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#end

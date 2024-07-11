# -*- coding: utf-8 -*-
# encoding: utf-8

from yade import utils,ymport,export,plot
import numpy as np
import os, sys, importlib, shutil

### running simulation in batch-mode
#utils.readParamsFromTable(INP=None,GEOMETRY=None,SLOPEINCL=None,recINT=None,ITERMAX=None,KINEB=None,MAV_INT=None,noTableOk=True) # define parameters from table
#from yade.params.table import *  # uncomment if using batch-mode with table

# test with no table (comment if table is used), ( yadedaily --cores '4' cutAndFail.py |& tee terminalOutput_cutAndFail.txt)
INP='tile_75-4-22_0.02_89212_isoLith19_noANI_loaded'
GEOMETRY='slopeStep'
SLOPEINCL=60
recINT=1000
ITERMAX=350000
KINEB=6e10
MAV_INT=3000


print('input file: ', INP)
print('geometry: ', GEOMETRY)
print('slope inclination: ', SLOPEINCL)


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##### first steps, defining controlling parameters and output

#### import yade-file, the loaded tile
O.load(INP+'.yade')


#### length of simulation
ITERMAX=ITERMAX
#ITERMAX=400000                                  # after topography is cut and slope is stable, for the strength reduction part)
print('maximum iteration steps: ', ITERMAX)


####  recording interval
recINT=recINT
#recINT=1000                                     # for vtk_files and position_files)
print('recording interval: ', recINT)


#### Scale Factor
scaleFactor=1000                                # redefine scale factor


#### Plateau Length plateauLengthRatio
plateauLengthRatio=2.0                          # sets the ratio of plateau length to slope height (which is 1*scaleFactor), applies only for slopeStep geometry (for now)
forelandLengthRatio=2.0                         # sets the ratio of foreland length to slope height (which is 1*scaleFactor), applies only for slopeStep geometry (for now)


#### stability criterion
kinE_bound=KINEB
#kinE_bound=6e10                                 # kinetic energy bound as stability criterion [J] (after gravity and weakness plane introduction)
print('kinetic energy bound:',"{:.1e}".format(kinE_bound),'[J]')


#### strength update interval (for strength reduction)
strengthUpdateInterval=100


#### predefining moving average characteristics (applied on kinetic energy and used for both stabilisation and strength reduction functions)
movAveInt=MAV_INT
#movAveInt=3000
movAveSamplingNb=int(movAveInt/strengthUpdateInterval)      # this defines the length of the moving average array, we have set it in a way that the moving average sampling step always 
                                                            #   matches the strength update interval, but it can be set to a fixed value, too
#movAveSamplingNb=10
movAveSamplingStep=int(movAveInt/movAveSamplingNb)          # do not change that

print('movAveInt=',movAveInt)
print('movAveSamplingNb=',movAveSamplingNb)
print('movAveSamplingStep=',movAveSamplingStep)


#### number of steps for carving out topography
carvingSteps=10
carvingSteps_overburden=2


#### extract information from input file name
extract=INP.rsplit('_', 5)
#print(extract)
Rpacking=extract[1]
print('Rpacking: ',Rpacking)     #defined particle radius at packing creation
NofParticles=extract[2]
print('NofParticles: ',NofParticles)     #number of particles of packing
lithology=extract[3]
print('lithology: ',lithology)
ANIINCL=extract[4]
print('weakness plane inclination: ', ANIINCL)


#### interparticle properties, load lithology
path=os.getcwd()
sys.path.append(path)
filename='%s' % lithology
lith_module = importlib.import_module(filename)
module_dict = lith_module.__dict__
try:
    to_import = lith_module.__all__
except AttributeError:
    to_import = [name for name in module_dict if not name.startswith('_')]
globals().update({name: module_dict[name] for name in to_import})


#### define name of output file
output=GEOMETRY+'_'+Rpacking+'_'+NofParticles+'_'+lithology+'_'+str(SLOPEINCL)+'_'+ANIINCL+'_'+str(movAveInt)
print('output=',output)


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


#### changing output in O.engines (see loadTile.py)
interactionLaw.Key=(output)
recVTK.fileName=(output+'/'+output+'_')
recVTK.Key=(output)


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##### a few helpful functions

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


#### thickness of overburden
hOfOverb=Z-2*scaleFactor
print('thickness of overburden :',hOfOverb)


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


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##### preliminary cutting

#### cutting to right length
print('FIRST CUTTING TO RIGHT LENGTH (NO TOPO YET)')
erased_1stStep=0
nbSpheres=0

if GEOMETRY=='slopeStep':
    for b in O.bodies:
        if isinstance(b.shape,Sphere):
            if b.state.pos[0] < -(scaleFactor*forelandLengthRatio) or b.state.pos[0] > (scaleFactor/(tan(radians(SLOPEINCL))))+(scaleFactor*plateauLengthRatio) :
                erased_1stStep+=1
                O.bodies.erase(b.id)
            else:
                nbSpheres+=1

if GEOMETRY=='ridge'or GEOMETRY=='valley':
    for b in O.bodies:
        if isinstance(b.shape,Sphere):
            if abs(b.state.pos[0]) > (scaleFactor/(tan(radians(SLOPEINCL))))+scaleFactor :
                erased_1stStep+=1
                O.bodies.erase(b.id)
            else:
                nbSpheres+=1

print('number of spheres erased in 1st carving step: ',erased_1stStep)
print('number of spheres remaining after 1st carving step: ',nbSpheres)


#### counting weakness plane bonds after cutting to right length
if ANI==True:
    nbWPB=0
    for i in O.interactions:
        if i.phys.isOnJoint==True :
            nbWPB+=1
    
    print('number of weakness plane bonds after cutting to right length:',nbWPB)


#### get dimensions of the packing, again
dimensionsPack()
print('xinf=',xinf,' | yinf=',yinf,' | zinf=',zinf)
print('xsup=',xsup,' | ysup=',ysup,' | zsup=',zsup)
print('X=',X,' | Y=',Y,' | Z=',Z)


#### Identify indicator particles for the failure simulation

if GEOMETRY=='slopeStep':

    refPoint1=0
    refPoint2=0
    refPoint3=0
    refPoint4=0
    refPoint5=0

    Xref1=(X-(scaleFactor*forelandLengthRatio)-(scaleFactor*plateauLengthRatio))
    Yref1=0
    Zref1=scaleFactor-2.*Rmax

    Xref2=(3/4)*(X-(scaleFactor*forelandLengthRatio)-(scaleFactor*plateauLengthRatio))
    Yref2=0
    Zref2=(3/4)*scaleFactor-2.*Rmax

    Xref3=(1/2)*(X-(scaleFactor*forelandLengthRatio)-(scaleFactor*plateauLengthRatio))
    Yref3=0
    Zref3=(1/2)*scaleFactor-2.*Rmax

    Xref4=(1/4)*(X-(scaleFactor*forelandLengthRatio)-(scaleFactor*plateauLengthRatio))
    Yref4=0
    Zref4=(1/4)*scaleFactor-2.*Rmax

    Xref5=0
    Yref5=0
    Zref5=0-2.*Rmax

if GEOMETRY=='ridge':

    refPoint1=0
    refPoint2=0
    refPoint3=0
    refPoint4=0
    refPoint5=0

    Xref1=-((X/2)-scaleFactor)
    Yref1=0
    Zref1=0-2.*Rmax

    Xref2=-(1/2)*((X/2)-scaleFactor)
    Yref2=0
    Zref2=(1/2)*scaleFactor-2.*Rmax

    Xref3=0
    Yref3=0
    Zref3=scaleFactor-2.*Rmax

    Xref4=(1/2)*((X/2)-scaleFactor)
    Yref4=0
    Zref4=(1/2)*scaleFactor-2.*Rmax

    Xref5=((X/2)-scaleFactor)
    Yref5=0
    Zref5=0-2.*Rmax

if GEOMETRY=='valley':

    refPoint1=0
    refPoint2=0
    refPoint3=0
    refPoint4=0
    refPoint5=0

    Xref1=-((X/2)-scaleFactor)
    Yref1=0
    Zref1=scaleFactor-2.*Rmax

    Xref2=-(1/2)*((X/2)-scaleFactor)
    Yref2=0
    Zref2=(1/2)*scaleFactor-2.*Rmax

    Xref3=0
    Yref3=0
    Zref3=0-2.*Rmax

    Xref4=(1/2)*((X/2)-scaleFactor)
    Yref4=0
    Zref4=(1/2)*scaleFactor-2.*Rmax

    Xref5=((X/2)-scaleFactor)
    Yref5=0
    Zref5=scaleFactor-2.*Rmax

for o in O.bodies:
    if isinstance(o.shape,Sphere):
        if o.state.pos[0]>(Xref1-Rmax) and o.state.pos[0]<(Xref1+Rmax) and o.state.pos[1]>(Yref1-Rmax) and o.state.pos[1]<(Yref1+Rmax) and o.state.pos[2]>(Zref1-Rmax) and o.state.pos[2]<(Zref1+Rmax) :
            refPoint1=o.id
            px1=o.state.pos[0]
            py1=o.state.pos[1]
            pz1=o.state.pos[2]
        if o.state.pos[0]>(Xref2-Rmax) and o.state.pos[0]<(Xref2+Rmax) and o.state.pos[1]>(Yref2-Rmax) and o.state.pos[1]<(Yref2+Rmax) and o.state.pos[2]>(Zref2-Rmax) and o.state.pos[2]<(Zref2+Rmax) :
            refPoint2=o.id
            px2=o.state.pos[0]
            py2=o.state.pos[1]
            pz2=o.state.pos[2]
        if o.state.pos[0]>(Xref3-Rmax) and o.state.pos[0]<(Xref3+Rmax) and o.state.pos[1]>(Yref3-Rmax) and o.state.pos[1]<(Yref3+Rmax) and o.state.pos[2]>(Zref3-Rmax) and o.state.pos[2]<(Zref3+Rmax) :
            refPoint3=o.id
            px3=o.state.pos[0]
            py3=o.state.pos[1]
            pz3=o.state.pos[2]
        if o.state.pos[0]>(Xref4-Rmax) and o.state.pos[0]<(Xref4+Rmax) and o.state.pos[1]>(Yref4-Rmax) and o.state.pos[1]<(Yref4+Rmax) and o.state.pos[2]>(Zref4-Rmax) and o.state.pos[2]<(Zref4+Rmax) :
            refPoint4=o.id
            px4=o.state.pos[0]
            py4=o.state.pos[1]
            pz4=o.state.pos[2]
        if o.state.pos[0]>(Xref5-Rmax) and o.state.pos[0]<(Xref5+Rmax) and o.state.pos[1]>(Yref5-Rmax) and o.state.pos[1]<(Yref5+Rmax) and o.state.pos[2]>(Zref5-Rmax) and o.state.pos[2]<(Zref5+Rmax) :
            refPoint5=o.id
            px5=o.state.pos[0]
            py5=o.state.pos[1]
            pz5=o.state.pos[2]

print('refPoint1=',refPoint1,' | Xref1=',Xref1,', Yref1=',Yref1,', Zref1=',Zref1)
O.bodies[refPoint1].shape.color=(1,1,0)
print('refPoint2=',refPoint2,' | Xref2=',Xref2,', Yref2=',Yref2,', Zref2=',Zref2)
O.bodies[refPoint2].shape.color=(1,1,0)
print('refPoint3=',refPoint3,' | Xref3=',Xref3,', Yref3=',Yref3,', Zref3=',Zref3)
O.bodies[refPoint3].shape.color=(1,1,0)
print('refPoint4=',refPoint4,' | Xref4=',Xref4,', Yref4=',Yref4,', Zref4=',Zref4)
O.bodies[refPoint4].shape.color=(1,1,0)
print('refPoint5=',refPoint5,' | Xref5=',Xref5,', Yref5=',Yref5,', Zref5=',Zref5)
O.bodies[refPoint5].shape.color=(1,1,0)


#### define boundary conditions, again
print('DEFINE BOUNDARY CONDITIONS')
e=3*Rmean
baseBodies=[]
for o in O.bodies:
    if isinstance(o.shape,Sphere):

        ### these ensure the boundaries to not move in their normal direction ("roller" like BC) -> should we just use the walls? -> easier for model rotation
        ## front particles
        if o.state.pos[0]<(xinf+e):
            o.state.vel=Vector3(0.,0.,0.)
            o.state.angVel=Vector3(0.,0.,0.)
            o.state.blockedDOFs='x'
            o.shape.color=(0,1,0)
        ### back particles
        if o.state.pos[0]>(xsup-e):
            o.state.vel=Vector3(0.,0.,0.)
            o.state.angVel=Vector3(0.,0.,0.)
            o.state.blockedDOFs='x'
            o.shape.color=(0,1,0)
        ## left particles
        if o.state.pos[1]<(yinf+e):
            o.state.vel=Vector3(0.,0.,0.)
            o.state.angVel=Vector3(0.,0.,0.)
            o.state.blockedDOFs+='y'
            o.shape.color=(0,0,1)
        ## right particles
        if o.state.pos[1]>(ysup-e):
            o.state.vel=Vector3(0.,0.,0.)
            o.state.angVel=Vector3(0.,0.,0.)
            o.state.blockedDOFs+='y'
            o.shape.color=(0,0,1)
        ## ground particles
        if o.state.pos[2]<(zinf+e):
            o.state.vel=Vector3(0.,0.,0.)
            o.state.angVel=Vector3(0.,0.,0.)
            o.state.blockedDOFs='xyzXYZ'
            o.shape.color=(1,0,0)
            baseBodies.append(o.id)


#### add wall to all the two vertical sides that got cut
def wallMat(): return JCFpmMat(type=0,density=2000,young=0.01*YOUNG,poisson=ALPHA,tensileStrength=0,cohesion=0,frictionAngle=radians(FRICT))
O.bodies.append(wall((xinf,yinf+Y/2.,zinf+Z/2.),0,sense=0,color=None,material=wallMat))
O.bodies.append(wall((xsup,yinf+Y/2.,zinf+Z/2.),0,sense=0,color=None,material=wallMat))


#### compute the vertical stress
sigmaZbase=0
def vertStress():
    global sigmaZbase
    sigmaZbase=utils.sumForces(baseBodies,(0,0,1))/(X*Y)

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
##### engine functions redefined (record data)

#### myRecorder function 
V=0
coeff=1
def myRecorder():
    global V
    for o in O.bodies :                             # |
        if abs(np.linalg.norm(o.state.vel))>0.01:   # | It calculates the volume that moves faster than 0.01 m/s
            V+=1.333*pi*pow(o.shape.radius,3)       # |
    rId1=O.bodies[refPoint1]
    rId2=O.bodies[refPoint2]
    rId3=O.bodies[refPoint3]
    rId4=O.bodies[refPoint4]
    rId5=O.bodies[refPoint5]
    plot.addData(Ek=utils.kineticEnergy(),v1=abs(np.linalg.norm(rId1.state.vel)),v2=abs(np.linalg.norm(rId2.state.vel)),v3=abs(np.linalg.norm(rId3.state.vel)),v4=abs(np.linalg.norm(rId4.state.vel)),v5=abs(np.linalg.norm(rId5.state.vel)),g=newton.gravity[2],px1=rId1.state.pos[0]-px1,py1=rId1.state.pos[1]-py1,pz1=rId1.state.pos[2]-pz1,px2=rId2.state.pos[0]-px2,py2=rId2.state.pos[1]-py2,pz2=rId2.state.pos[2]-pz2,px3=rId3.state.pos[0]-px3,py3=rId3.state.pos[1]-py3,pz3=rId3.state.pos[2]-pz3,px4=rId4.state.pos[0]-px4,py4=rId4.state.pos[1]-py4,pz4=rId4.state.pos[2]-pz4,px5=rId5.state.pos[0]-px5,py5=rId5.state.pos[1]-py5,pz5=rId5.state.pos[2]-pz5,iterations=O.iter,t=O.realtime,redF=coeff,vol=V,tc=interactionLaw.nbTensCracks,sc=interactionLaw.nbShearCracks,mA_Ek=movAve_kinE)
    plot.saveDataTxt(output+'/'+output)
    V=0                                  #reset Volume calculation


#### record position-files
def positions():
    export.text(output+'/'+posDir+'/'+output+'-pos_'+str(O.iter))


#### moving average of kinetic energy function
# redefining numpy array
movAve_kinE=0.
dynList_kinE=np.empty(movAveSamplingNb)
dynList_kinE.fill(kinE_bound*1.1)

# redefine the moving average sampling step
movAveKinE.iterPeriod=movAveSamplingStep

# the function redefined
def movAverage_kinEnergy():
    global dynList_kinE, movAve_kinE
    kinE=utils.kineticEnergy()
    
    dynList_kinE=np.append(dynList_kinE,kinE)
    dynList_kinE=np.delete(dynList_kinE,0)
    movAve_kinE=np.mean(dynList_kinE)


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##### add anti-infinite-fall engine: takes out particles if necessary to avoid artefact of particles flying into infinity and messing up strength reduction 

#### loop that sets particles "non-dynamic" if they went trough bottom wall
def antiInfiniteFall():
    for o in O.bodies:                                  # (maybe slows down calculation significantly, comment if not necessary)
        if o.state.pos[2]<(zinf-(0.1*scaleFactor)):
            o.dynamic=False

#### adding anti-infinite-fall function to engines
O.engines=O.engines+[PyRunner(command='antiInfiniteFall()',label='antiInfinFall',dead=1)]

#### anti-infinite-fall intervals
antiInfinFall.dead=0
antiInfinFall.iterPeriod=1000


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##### simulation begins

#### stabilisation and a VTK recording

print('SIMULATION BEGINS')

stabilisation()

vertStress()
print('vertStress = ',sigmaZbase)

recVTK.dead=0
O.step()
recVTK.dead=1


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##### removing overburden

#### cutting off overburden
print('CUTTING OFF OVERBURDEN (NO TOPO YET)')
erased_2ndStep=0

for c in range(1,carvingSteps_overburden+1):
    #for c in range(1,3):
    print('Overburden layer No. ',c,' of ',carvingSteps_overburden, ' is carved out!')
    for b in O.bodies:
        if isinstance(b.shape,Sphere):
            if b.state.pos[2] > zsup-(hOfOverb/carvingSteps_overburden)*c :
                erased_2ndStep+=1
                O.bodies.erase(b.id)

    stabilisation()

    vertStress()
    print('vertStress = ',sigmaZbase)

    recVTK.dead=0
    O.step()
    recVTK.dead=1

nbSpheres=nbSpheres-erased_2ndStep
print('number of spheres erased in 2nd carving step: ',erased_2ndStep)
print('number of spheres remaining after 2nd carving step: ',nbSpheres)


#### get dimensions of the packing, again
dimensionsPack()
print('xinf=',xinf,' | yinf=',yinf,' | zinf=',zinf)
print('xsup=',xsup,' | ysup=',ysup,' | zsup=',zsup)
print('X=',X,' | Y=',Y,' | Z=',Z)


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##### topo cutting

#### start carving
print('CARVING OUT TOPOGRAPHY')
erased_3rdStep=0
#nbSpheres=0
#carvingSteps=10


for c in range(1,carvingSteps+1):
    #for c in range(1,3):
    print('Layer No. ',c,' of ',carvingSteps, ' is carved out!')

    if GEOMETRY=='slopeStep':
        for b in O.bodies:
            if isinstance(b.shape,Sphere):
                if (b.state.pos[0] < 0 and b.state.pos[2] > scaleFactor-(scaleFactor/carvingSteps)*c) or (b.state.pos[0] > 0 and b.state.pos[2] > tan(radians(SLOPEINCL))*b.state.pos[0] and b.state.pos[2] > scaleFactor-(scaleFactor/carvingSteps)*c) :
                    erased_3rdStep+=1
                    O.bodies.erase(b.id)

    if GEOMETRY=='ridge':
        for b in O.bodies:
            if isinstance(b.shape,Sphere):
                if (b.state.pos[2] > -tan(radians(SLOPEINCL))*abs(b.state.pos[0])+scaleFactor and b.state.pos[2] > scaleFactor-(scaleFactor/carvingSteps)*c) :
                    erased_3rdStep+=1
                    O.bodies.erase(b.id)

    if GEOMETRY=='valley':
        for b in O.bodies:
            if isinstance(b.shape,Sphere):
                if (b.state.pos[2] > tan(radians(SLOPEINCL))*abs(b.state.pos[0]) and b.state.pos[2] > scaleFactor-(scaleFactor/carvingSteps)*c) :
                    erased_3rdStep+=1
                    O.bodies.erase(b.id)

    stabilisation()

    vertStress()
    print('vertStress = ',sigmaZbase)

    recVTK.dead=0
    O.step()
    recVTK.dead=1

print('TOPOGRAPHY IS CARVED OUT')

nbSpheres=nbSpheres-erased_3rdStep
print('number of spheres erased in 3rd carving step: ',erased_3rdStep)

erased_total=erased_1stStep+erased_2ndStep+erased_3rdStep
print('total number of spheres erased: ',erased_total)

print('number of spheres remaining after 3rd carving step: ',nbSpheres)


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##### delete potential "floating particles"

##### loop through particles
#Kmin=1                                                              # minimum number of bonds to keep the particle in the simulation
#erased_floating=0
#for o in O.bodies:
    ##if not o : continue 
    #nbCont=0
    #for i in O.interactions.withBody(o.id) :
        #if not i.isReal : continue
        #if i.phys.isCohesive :
            #nbCont+=1
    #if nbCont<Kmin :
        #erased_floating+=1
        #O.bodies.erase(o.id)


##### another recording step and output

#recVTK.dead=0
#O.step()
#recVTK.dead=1

#print('FLOATING PARTICLES ARE REMOVED')

#nbSpheres=nbSpheres-erased_floating
#print('number of spheres erased when removing floating particles: ',erased_floating)
#erased_total=erased_1stStep+erased_2ndStep+erased_3rdStep+erased_floating
#print('total number of spheres erased: ',erased_total)
#print('number of spheres remaining after carving and removing floating particles ',nbSpheres)


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##### apply final strength of weakness plane bonds

#### the loop through the particles and bonds
if ANI==True:
    print('APPLY FINAL STRENGTH OF WEAKNESS PLANE BONDS')
    nbWPB=0
    for i in O.interactions:
        if i.phys.isOnJoint==True :
            # bodies
            O.bodies[i.id1].mat.jointTensileStrength=STENS
            O.bodies[i.id2].mat.jointTensileStrength=STENS
            O.bodies[i.id1].mat.jointCohesion=SCOH
            O.bodies[i.id2].mat.jointCohesion=SCOH
            O.bodies[i.id1].mat.jointFrictionAngle=radians(SFRIC)
            O.bodies[i.id2].mat.jointFrictionAngle=radians(SFRIC)
            O.bodies[i.id1].mat.jointDilationAngle=radians(SDIL)
            O.bodies[i.id2].mat.jointDilationAngle=radians(SDIL)
            # interactions
            i.phys.FnMax=STENS
            i.phys.FsMax=SCOH
            i.phys.tanFrictionAngle=tan(radians(SFRIC))
            i.phys.tanDilationAngle=tan(radians(SDIL))
            
            nbWPB+=1
    
    print('number of weakness plane bonds after topography is carved out:',nbWPB)
    
    ## stabilisation and a VTK recording
    stabilisation()

    vertStress()
    print('vertStress = ',sigmaZbase)

    recVTK.dead=0
    O.step()
    recVTK.dead=1


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##### add the strength reduction engine

#### new strength reduction function
print('STRENGTH REDUCTION BEGINS')
coeff=1
printInt=1                              #printing interval for stabilisation loop, printInt=1 prints out every sampling step
count=printInt

nbCracksRef=interactionLaw.nbTensCracks+interactionLaw.nbShearCracks

def strengthUpdate():
    global count, coeff, nbCracksRef
    nbCracks=interactionLaw.nbTensCracks+interactionLaw.nbShearCracks
    kinE=utils.kineticEnergy()

    if kinE>=kinE_bound or movAve_kinE>=kinE_bound or nbCracks>nbCracksRef:
        count-=1
        nbCracksRef=nbCracks
        if count==0:
            print('iter= ',O.iter,', cracking/failing!: nbCracks=', nbCracks, ', kinE=', round(kinE), ', movAve_kinE=', round(movAve_kinE))
            count=printInt

    elif kinE<kinE_bound and movAve_kinE<kinE_bound and nbCracks==nbCracksRef and coeff!=0:
        coeff*=0.9
        if coeff<1e-4 : coeff=0         #if coeff<0.001 : coeff=0 (earlier used bound)
        print('iter= ',O.iter,' || strength degradation! reductionFactor = ',coeff)
        for i in O.interactions:
            if i.phys.isOnJoint==True : continue
            #if coeff==0 : i.phys.isCohesive=False
            i.phys.FnMax*=0.9
            i.phys.FsMax*=0.9

    elif kinE<kinE_bound and movAve_kinE<kinE_bound and nbCracks==nbCracksRef and coeff==0:    # simulation stops when coeff==0 and simulation is stable
        print('simulation is stable and coeff==0, simulation stops')
        O.pause()


#### adding strength reduction function to engines
O.engines=O.engines+[PyRunner(command='strengthUpdate()',label='strengthRed',dead=1)]

#### anti-infinite-fall intervals
strengthRed.dead=1
strengthRed.iterPeriod=strengthUpdateInterval


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##### topography is cut! run strength reduction simulation

#### activate recorders
recPositions.dead=0
recVTK.dead=0

#### recording intervals
recPositions.iterPeriod=recINT
recVTK.iterPeriod=recINT

#### activate strength reduction
strengthRed.dead=0

#### start strength reduction simulation
O.run(int(ITERMAX),True)


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


#### moving vtk-files
for f in VTKfiles:
    if (f.startswith(output + '_cracks')) :
        shutil.move(os.path.abspath(f), cracksPath)
    elif (f.startswith(output + '_intrs')) :
        shutil.move(os.path.abspath(f), intrsPath)
    elif (f.startswith(output + '_spheres')) :
        shutil.move(os.path.abspath(f), spheresPath)

#### move cracks text-file
shutil.move(path + '/cracks_'+output+'.txt', cracksPath)
os.chdir(path)


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#end

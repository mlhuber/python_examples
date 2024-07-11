#!/usr/bin/python
# -*- coding: utf-8 -*-

import gts, os.path, locale

O=Omega()
from yade import utils,pack,ymport,export

utils.readParamsFromTable(MESH=None,R=None,noTableOk=True) # batch mode: define parameters from table
from yade.params.table import *

## open guis
#from yade import qt
#v=qt.Controller()
#v=qt.View()


#--------------------------------------------------------------------------------------------------------------------------------
##### set up parameters

#### controlling parameters
#mesh='tile_75-4-22.gts'      #name of gts mesh, e.g. 'slope60_long.gts'; don't forget to adjust the total volume calculation (totVol) when the input geometry is changing
mesh=MESH
fileStem=mesh.replace('.gts', '')
#sphereRad=0.020             #directly define mean sphere radius, e.g. 0.035, check size ratio on output
sphereRad=R
SphInC=20000               #spheres in cell


#### material definition
def sphereMat(): return FrictMat(density=2500,young=1e9,poisson=0.3,frictionAngle=radians(0))


#### finalisation criterion
comp_bound=0.7              #compacity as a finalisation criterion


#--------------------------------------------------------------------------------------------------------------------------------
##### create the pacling in yade

#### import mesh
locale.setlocale(locale.LC_ALL,'en_US.UTF-8')   #gts is locale-dependend !!
surface=gts.read(open(mesh))
print('closed? ', surface.is_closed())


#### generate packing
if surface.is_closed():
    pred=pack.inGtsSurface(surface)
    # get characteristic dimensions
    aabb=pred.aabb()
    dim=pred.dim()
    center=pred.center()
    ### packing functions
    #O.bodies.append(pack.regularHexa(pred,radius=sphereRad,gap=0.,color=(0.9,0.8,0.6),material=sphereMat))
    #O.bodies.append(pack.regularOrtho(pred,radius=sphereRad,gap=0.,color=(0.9,0.8,0.6),material=sphereMat))

    #O.bodies.append(pack.randomDensePack(pred,radius=sphereRad,rRelFuzz=0.333,spheresInCell=3000,memoizeDb='/tmp/gts-triax-packings.sqlite',returnSpherePack=False,color=(0.9,0.8,0.6),material=sphereMat))

    O.bodies.append(pack.randomDensePack(pred,radius=sphereRad,rRelFuzz=0.333,spheresInCell=SphInC,returnSpherePack=False,color=(0.9,0.8,0.6),material=sphereMat))


#### gtsSurface2Facets
O.bodies.append(pack.gtsSurface2Facets(surface,wire=True,material=sphereMat))           # for what is this good for ?


#--------------------------------------------------------------------------------------------------------------------------------
##### some basic stuff

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
            o.shape.color=(0.7,0.5,0.3)
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


#### calculate size ratio
sizeRatio=min(X,Y,Z)/(Rmean*2) #defines discretisation of smallest packing dimension (sizeRatio=meshLength/particleDiameter)
print('sizeRatio=',sizeRatio)


#### compute volume function
def volume():
    global packingVolume, volSpheres, comp
    packingVolume=0
    volSpheres=0

    packingVolume=X*Y*Z
    for o in O.bodies:
        if isinstance(o.shape,Sphere):
            volSpheres+=(4./3.)*pi*(o.shape.radius)**3.

volume()
print('packingVolume=',packingVolume)
print('volSpheres=',volSpheres)

comp=volSpheres/packingVolume
print('initial compacity=',comp)


#--------------------------------------------------------------------------------------------------------------------------------
##### define engine and time step

#### engines
O.engines=[ForceResetter(),
           InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Facet_Aabb()]),
           InteractionLoop(
               [Ig2_Sphere_Sphere_ScGeom(),Ig2_Facet_Sphere_ScGeom()],
               [Ip2_FrictMat_FrictMat_FrictPhys()],
               [Law2_ScGeom_FrictPhys_CundallStrack()]
           ),
           NewtonIntegrator(damping=0.5,label='newton')
]


#### define one time step
O.dt=0.1*utils.PWaveTimeStep()                     # for what is this good for ???
print('O.dt=',O.dt)
O.step()


#--------------------------------------------------------------------------------------------------------------------------------
##### run the packing alteration loop

#### dropping in the loop
while 1:
    volSpheres=0.
    packingVolume=0.
    for o in O.bodies:
        if isinstance(o.shape,Sphere):
            o.shape.radius*=1.001
    dimensionsPack()
    volume()
    comp=volSpheres/packingVolume
    O.run(100,True)
    unb=unbalancedForce()
    print('unbF:',unb,' compacity: ',comp,' iter:',O.iter)
    if comp>=comp_bound:
        print('final compacity=',comp)
        O.run(1000,True)
        break


#### final values
print('FINAL VALUES')
print('xinf=',xinf,' | yinf=',yinf,' | zinf=',zinf)
print('xsup=',xsup,' | ysup=',ysup,' | zsup=',zsup)
print('X=',X,' | Y=',Y,' | Z=',Z)
dimensionsSph()
print('Rmax=',Rmax,' | Rmin=',Rmin,' | Rmean=',Rmean)
print('packingVolume=',packingVolume)
print('volSpheres=',volSpheres)


#### export treated packing
export.text(fileStem+'_'+str(sphereRad)+'_'+str(int(nbSpheres))+'.spheres')


#--------------------------------------------------------------------------------------------------------------------------------
## end

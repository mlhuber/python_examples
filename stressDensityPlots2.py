# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import gaussian_kde



#### data import and selecting area of interest

## load csv-files from paraview as pandas dataframes

sphereDataFileName = [
    'slopeStep_0.02_89212_aniLith19_60_0_3000_spheresDataFromClip.csv',
    'slopeStep_0.02_89212_aniLith19_60_10_3000_spheresDataFromClip.csv',
    'slopeStep_0.02_89212_aniLith19_60_20_3000_spheresDataFromClip.csv',
    'slopeStep_0.02_89212_aniLith19_60_30_3000_spheresDataFromClip.csv',
    'slopeStep_0.02_89212_aniLith19_60_40_3000_spheresDataFromClip.csv',
    'slopeStep_0.02_89212_aniLith19_60_50_3000_spheresDataFromClip.csv',
    'slopeStep_0.02_89212_aniLith19_60_60_3000_spheresDataFromClip.csv',
    'slopeStep_0.02_89212_aniLith19_60_70_3000_spheresDataFromClip.csv',
    'slopeStep_0.02_89212_aniLith19_60_80_3000_spheresDataFromClip.csv',
    'slopeStep_0.02_89212_aniLith19_60_90_3000_spheresDataFromClip.csv',
    'slopeStep_0.02_89212_aniLith19_60_100_3000_spheresDataFromClip.csv',
    'slopeStep_0.02_89212_aniLith19_60_110_3000_spheresDataFromClip.csv',
    'slopeStep_0.02_89212_aniLith19_60_120_3000_spheresDataFromClip.csv',
    'slopeStep_0.02_89212_aniLith19_60_130_3000_spheresDataFromClip.csv',
    'slopeStep_0.02_89212_aniLith19_60_140_3000_spheresDataFromClip.csv',
    'slopeStep_0.02_89212_aniLith19_60_150_3000_spheresDataFromClip.csv',
    'slopeStep_0.02_89212_aniLith19_60_160_3000_spheresDataFromClip.csv',
    'slopeStep_0.02_89212_aniLith19_60_170_3000_spheresDataFromClip.csv',
    'slopeStep_0.02_89212_isoLith19_60_noANI_3000_spheresDataFromClip.csv'
    ]

totalNbOfFiles = len(sphereDataFileName)


aniAngleList = []
for i in range(totalNbOfFiles):
    extract=sphereDataFileName[i].rsplit('_', 3)
    aniAngleList.append(extract[1])

#print(aniAngleList)

for i in range(totalNbOfFiles):
    globals()[f'sphereData{i}'] = pd.read_csv(sphereDataFileName[i], sep=',')


## choose only spheres in certain area

#for i in range(totalNbOfFiles):
    #globals()[f'sphereData{i}'] =  globals()[f'sphereData{i}'].loc[(globals()[f'sphereData{i}']['Points:0'] >= -500) & \
        #(globals()[f'sphereData{i}']['Points:0'] <= 1250) & \
            #(globals()[f'sphereData{i}']['Points:2'] >= -500)]

#sphereData1 = sphereData1.loc[(sphereData1['Points:0'] >= -500) & (sphereData1['Points:0'] <= 1250) & (sphereData1['Points:2'] >= -500)]
#sphereData1 = sphereData1.loc[(sphereData1['Points:0'] >= 0) & (sphereData1['Points:0'] <= 1000) & (sphereData1['Points:2'] >= -500)]


## import list for compacity values (1 - porosity) from the loading of each anisotropy orientation configuration

compacityTable_fileName = 'compacity.csv'
compacityTable = pd.read_csv(compacityTable_fileName, sep=',')

compacity = compacityTable['compacity']

#print(compacity)
#print(compacity[5])


#### interpolation step for drawing contours
#sources:   |https://stackoverflow.com/questions/70646328/how-to-create-a-polar-density-plot-having-a-list-of-angles-values-in-degrees-in
#           |https://stackoverflow.com/questions/31051882/matplotlib-density-plot-in-polar-coordinates/31142319#31142319

### assign variables for interpolation

## angle-values put to radians

for i in range(totalNbOfFiles):
    globals()[f'angles{i}'] = np.deg2rad(globals()[f'sphereData{i}']['dirIII_angle '])

#angles1 = np.deg2rad(sphereData1['dirIII_angle '])


## set which stress value you want to plot

for i in range(totalNbOfFiles):
    globals()[f'radii{i}'] = globals()[f'sphereData{i}']['diffStress'] * compacity[i]
stressValue = 'diffStress'

#radii1 = sphereData1['diffStress']

#for i in range(totalNbOfFiles):
    #globals()[f'radii{i}'] = abs(globals()[f'sphereData{i}']['sigIII'])
#stressValue = 'sigma1'

#for i in range(totalNbOfFiles):
    #globals()[f'radii{i}'] = abs(globals()[f'sphereData{i}']['meanStress'])
#stressValue = 'meanStress'

#radii1 = abs(sphereData1['sigIII'])


### interpolation step

#N = 20 # number of points used to evaluate the interpolation
#N = 50 # number of points used to evaluate the interpolation
N = 100 # number of points used to evaluate the interpolation
interpFactor_divisor = 1
#interpFactor_divisor = 3

for i in range(totalNbOfFiles):
    globals()[f'interp{i}'] = gaussian_kde(np.vstack((globals()[f'angles{i}'], globals()[f'radii{i}'])), bw_method='scott')
    globals()[f'interp{i}'].set_bandwidth(bw_method=globals()[f'interp{i}'].factor / interpFactor_divisor)
    

#interp1 = gaussian_kde(np.vstack((angles1, radii1)), bw_method='scott',)
#interp1.set_bandwidth(bw_method=interp1.factor / interpFactor_divisor)


## define ranges for interpolation

angleRange = np.linspace(0, 0.5*np.pi, N)

#upperRbound = 8e7
upperRbound = 4e7
#upperRbound = 6e7
radiRange = np.linspace(0, upperRbound, N)


## to compute colors, first get the meshgrid of angles and radii with shape N x N x 2

mesh = np.stack(np.meshgrid(angleRange, radiRange), 0)

## then, compute N x N matrix of colors using `interp`, I think we need to reshape to accomodate Gaussian KDE API, but maybe this can be avoided

for i in range(totalNbOfFiles):
    globals()[f'Z{i}'] = globals()[f'interp{i}'](mesh.reshape(2, -1)).reshape(N, N)
    

##_______________________________________________________________________________________________

#### plotting

plt.rcParams.update({'figure.autolayout': True,'axes.labelsize':14,'legend.fontsize':14,'xtick.labelsize':12,'ytick.labelsize':12})    
#plt.rcParams.update({'legend.numpoints':1,'font.size': 16,'axes.labelsize':16,'xtick.major.pad':10,'ytick.major.pad':10,'legend.fontsize':16,'xtick.labelsize':12,'ytick.labelsize':12})
#lw=2
#ms=10
#levels = np.linspace(0, 8e-8, 6)
#levels = np.linspace(0, 5e-8, 6)
levels = np.linspace(0, 7e-8, 6)
#print(levels)
#levels = np.linspace(0, 1, 10)
alphaValue = 1

nbOfSubFigures = 6
nbOfFigures = (totalNbOfFiles//nbOfSubFigures)+1
#print(nbOfFigures)


for i in range(1,nbOfFigures+1):
    globals()[f'fig{i}'] = plt.figure(figsize=(20,10))
    for j in range(i*nbOfSubFigures-nbOfSubFigures,i*nbOfSubFigures):
        if j <= totalNbOfFiles-1:
            globals()[f'ax{j}'] = globals()[f'fig{i}'].add_subplot(2,3,(j+1)-(nbOfSubFigures*(i-1)),projection='polar')
            globals()[f'ax{j}'].set_thetamin(0)
            globals()[f'ax{j}'].set_thetamax(90)
            #globals()[f'ax{k}'].grid(False)
            ctf = globals()[f'ax{j}'].contourf(angleRange, radiRange, globals()[f'Z{j}'], alpha = alphaValue, levels=levels, extend='max', cmap = 'Greys')
            plt.colorbar(ctf)
            globals()[f'ax{j}'].set_title('α = ' + aniAngleList[j] + '°',fontweight='bold', size = 20)
            globals()[f'ax{j}'].set_ylabel('σ$_{1}$-σ$_{3}$ [Pa]')
            #globals()[f'ax{j}'].set_rticks([0, 2e7, 4e7, 6e7, 8e7])
            globals()[f'ax{j}'].set_rticks([0, 1e7, 2e7, 3e7, 4e7])
            globals()[f'ax{j}'].set_xticks(np.deg2rad([0,45,90]))
    globals()[f'fig{i}'].savefig('stressDensityPlots_' + stressValue + '_interPolFactor' + str(interpFactor_divisor) + '_' + 'fig' + str(i) + '.eps',dpi=400,format='eps',transparent=False)
    globals()[f'fig{i}'].savefig('stressDensityPlots_' + stressValue + '_interPolFactor' + str(interpFactor_divisor) + '_' + 'fig' + str(i) + '.png',dpi=400,format='png',transparent=False)    


##----------------------------------------------------------------------------------------------
#### show
plt.show()

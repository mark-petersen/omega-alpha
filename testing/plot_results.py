#!/usr/bin/env python
from __future__ import absolute_import, division, print_function, unicode_literals
"""
Plot vertical profiles from individual columns
Mark Petersen
September 2022
"""

############################## model files, run dirs
runDir = './'
initFileName = 'init.nc'
outputFileName = 'output.nc'

############################## domain and IC
nx = 32
ny = 4
import numpy as np
dc = 1.0
Lx = nx*dc
Ly = ny*dc*np.sqrt(3)/2
kx = 1
ky = 0
nVertLevels=1

#print('loading libraries...')
from datetime import date
import matplotlib.pyplot as plt
import xarray as xr
initDS = xr.open_dataset(runDir+initFileName)
outDS = xr.open_dataset(runDir+outputFileName)

nCells = initDS.dims['nCells']
N = int(np.sqrt(nCells))
xCell = initDS.variables['xCell']
yCell = initDS.variables['yCell']
nEdges = initDS.dims['nEdges']
xEdge = initDS.variables['xEdge']
yEdge = initDS.variables['yEdge']
angleEdge = initDS.variables['angleEdge']
nVertices = initDS.dims['nVertices']
xVertex = initDS.variables['xVertex']
yVertex = initDS.variables['yVertex']

k=0
iTime = 1
figdpi = 200
fig = plt.figure(figsize=(40,16))
varNames = ['normalVelocity', 'divergenceSol', 'relativeVorticitySol', 'relativeVorticityCellSol','del2GradDivVelocitySol','del2GradVortVelocitySol','del2VelocitySol',
            'normalVelocity', 'divergence', 'relativeVorticity','relativeVorticityCell','del2GradDivVelocityTendency','del2GradVortVelocityTendency','hmixDel2VelocityTendency']
nVars = len(varNames)
loc = ['edge','cell','vertex','cell','edge','edge','edge',
       'edge','cell','vertex','cell','edge','edge','edge']
f = ['in','in','in','in','in','in','in','in','out','out','out','out','out','out']
norm = [1.0,1e-5,1e-5,1e-5,4e-11,4e-11,4e-11,
        1.0,1e-5,1e-5,1e-5,4e-11,4e-11,4e-11]
err = np.zeros(nVars)
rms = np.zeros(nVars)

size = int(32./(N/16.))
for j in range(nVars):
    ax = plt.subplot(3,7,j+1)

    if f[j]=='in':
        var = initDS.variables[varNames[j]][:,k] / norm[j]
    elif f[j]=='out':
        var = outDS.variables[varNames[j]][iTime,:,k] / norm[j]
    if loc[j]=='cell':
        im = plt.scatter(xCell/1000,yCell/1000,c=var,s=size,marker='H',cmap=plt.cm.jet)
    elif loc[j]=='edge':
        im = plt.scatter(xEdge/1000,yEdge/1000,c=var,s=size,marker='s',cmap=plt.cm.jet)
    elif loc[j]=='vertex':
        im = plt.scatter(xVertex/1000,yVertex/1000,c=var,s=size,marker='^',cmap=plt.cm.jet)
    plt.colorbar(im, ax=ax)
    if j>=nVars/2:
        varSol = initDS.variables[varNames[j-int(nVars/2)]][:,k] / norm[j]
        diff = (var - varSol)
        #diff = varSol/var
        #print(varSol/var)
        err[j] = np.max(abs(diff)) # divide by max value
        rms[j] = np.sqrt(np.mean(diff**2))
        #print('nx={} '.format(np.sqrt(nCells))+'maxabs: {:9.2E}'.format(np.max(abs(diff))),varNames[j])
        plt.title('model: '+varNames[j])

        ax = plt.subplot(3,7,j+1+7)
        if loc[j]=='cell':
            im = plt.scatter(xCell/1000,yCell/1000,c=diff,s=size,marker='H',cmap=plt.cm.jet)
        elif loc[j]=='edge':
            im = plt.scatter(xEdge/1000,yEdge/1000,c=diff,s=size,marker='s',cmap=plt.cm.jet)
        elif loc[j]=='vertex':
            im = plt.scatter(xVertex/1000,yVertex/1000,c=diff,s=size,marker='^',cmap=plt.cm.jet)
        plt.colorbar(im, ax=ax)
        plt.title(varNames[j]+' {:9.2E}'.format(err[j]))
    else:
        plt.title('exact: '+varNames[j])
print('   {}, {:9.2E}, {:9.2E}, {:9.2E}, {:9.2E}, {:9.2e}, {:9.2e},     {:9.2E}, {:9.2E}, {:9.2E}, {:9.2E}, {:9.2e}, {:9.2E}' \
     .format(N,err[8],err[9], err[10],err[11],err[12],err[13], \
               rms[8],rms[9], rms[10],rms[11],rms[12],rms[13]))
figfile = 'plot_del2_nx{:04d}_t{}'.format(N,iTime)+'.png'
plt.savefig(figfile, bbox_inches='tight')
plt.close()

exit()
# compute k x grad nu from exact nu
# This is the Fortran code:
#      do iEdge = 1, nEdgesOwned
#         vertex1 = verticesOnEdge(1,iEdge)
#         vertex2 = verticesOnEdge(2,iEdge)
#         dvEdgeInv = 1.0_RKIND / dvEdge(iEdge)
#      del2GradVortVelocityTendency(k,iEdge) = - edgeMask(k,iEdge)*visc2 *&
#                   (relVort(k,vertex2)-relVort(k,vertex1))*dvEdgeInv

verticesOnEdge = initDS.variables['verticesOnEdge']
dvEdge = initDS.variables['dvEdge']
relativeVorticitySol = initDS.variables['relativeVorticitySol']
del2GradVortVelocitySol = initDS.variables['del2GradVortVelocitySol']
kxGradExactVort = np.zeros([nEdges,1])
#print(verticesOnEdge.shape)
#for iEdge in range(nEdges):
#     vertex1 = verticesOnEdge[iEdge,0]-1
#     vertex2 = verticesOnEdge[iEdge,1]-1
#     kxGradExactVort[iEdge,k] = - ( \
#         relativeVorticitySol[vertex2,k] - \
#         relativeVorticitySol[vertex1,k])/dvEdge[iEdge]
#print('after iEdge loop')
kxGradExactVortDiff = kxGradExactVort - del2GradVortVelocitySol
kxGradExactVortErr = np.max(abs(kxGradExactVortDiff/4e-10)) # divide by max value
kxGradExactVortRms = np.sqrt(np.mean(kxGradExactVortDiff**2))/4e-10

print('   {}, {:9.2E}, {:9.2E}, {:9.2E}, {:9.2E}, {:9.2E}, {:9.2E}, {:9.2E},     {:9.2E}, {:9.2E}, {:9.2E}, {:9.2E}, {:9.2E}, {:9.2E}, {:9.2E}' \
     .format(N, err[7],err[8],err[9],err[10],err[11],err[12],err[13], \
                rms[7],rms[8],rms[9],rms[10],rms[11],rms[12],rms[13]) )
#, kxGradExactVortRms))

figdpi = 200
fig = plt.figure(figsize=(10,16))
ax = plt.subplot(3,1,1)
im = plt.scatter(xEdge/1000,yEdge/1000,c=del2GradVortVelocitySol/4e-10,s=size,marker='s',cmap=plt.cm.jet)
plt.title('del2GradVortVelocitySol')
plt.colorbar(im, ax=ax)

ax = plt.subplot(3,1,2)
im = plt.scatter(xEdge/1000,yEdge/1000,c=kxGradExactVort/4e-10,s=size,marker='s',cmap=plt.cm.jet)
plt.title('kxGradExactVort')
plt.colorbar(im, ax=ax)

ax = plt.subplot(3,1,3)
im = plt.scatter(xEdge/1000,yEdge/1000,c=kxGradExactVortDiff/4e-10,s=size,marker='s',cmap=plt.cm.jet)
plt.title('kxGradExactVortDiff')
plt.colorbar(im, ax=ax)

figfile = 'plot_del2GradVort_nx{:04d}'.format(N)+'.png'
plt.savefig(figfile) #, bbox_inches='tight')
plt.close()

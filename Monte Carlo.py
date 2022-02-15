#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 24 11:02:41 2021

@author: feryantama
"""

from random import random, choices
from math import pi, sin, cos, log, acos
from matplotlib import pyplot as plt

import plotly.graph_objects as go
import plotly.io as pio
import time

# ============================================================================
def read_gamma_CS(filename,sep = '\t'):
    CS = {'Energy':[],'Coherent':[],'Incoherent':[],'Photoelectric':[],
          'Nuclear Pair Production':[],'Electron Pair Production':[],
          'Total w/ Coherent':[],'Total wo/ Coherent':[]}
    with open(filename) as f:
        for n,line in enumerate(f):
            if n > 2 and line!='\n':
                for i in range(0,8):
                    CS[list(CS.keys())[i]] += [float(line.split(sep=sep)[i])]
    return CS

# ============================================================================
def Interpolate(D,E,column):
    i=1
    while E-D['Energy'][i]>=0:
        i+=1
    j=i-1
    intpl = (E-D['Energy'][i])/(D['Energy'][j]-D['Energy'][i])*(D[column][j]-D[column][i])+D[column][i]
    return intpl

# ============================================================================
def Finite_Cylinder(P,C):
    H1,H2,r,a,b = P
    x,y,z       = C
    if (x-a)**2+(y-b)**2-r**2 <= 0 and z < H2 and z >= H1:
        p = -1
    else:
        p = 1
    return p

# ============================================================================
def print_out(Collection,n):
    for j in range(0,len(Collection[n])):
        if j==0:
            print('\n|================================== Particle '+str(n)+' ===============================|\n'+
                  ' hit\t'+" x\t\t\t"+"y\t\t\t"+"z\t\t\t"+"E\t\t\t"+"lambda\t\t"+" particle\t"+"\n"+
                  '|===|===========|============|===========|===========|=============|==========|')     
        print(j,'\t',
              "{:.3e}\t".format(Collection[n][j]['x']),
              "{:.3e}\t".format(Collection[n][j]['y']),
              "{:.3e}\t".format(Collection[n][j]['z']),
              "{:.3e}\t".format(Collection[n][j]['E']),
              "{:.3e}\t".format(Collection[n][j]['l']),
              (Collection[n][j]['p']),'\t',
              '\t  ',(Collection[n][j]['C']))
    Trace = [go.Scatter3d(x=[Collection[n][j]['x'] for j in range(0,len(Collection[n]))],
                          y=[Collection[n][j]['y'] for j in range(0,len(Collection[n]))],
                          z=[Collection[n][j]['z'] for j in range(0,len(Collection[n]))],
                          marker= dict(size=2),
                          line  = dict(color='darkblue',width=1))]
    return Trace

# ============================================================================
def Plot(P,BB,renderer = 'firefox',Trace = []):
    H1,H2,r,a,b = P
    
    import numpy as np
    pio.renderers.default = renderer
    
    xa = [r*cos(2*pi*i/100)+a for i in range(0,100)]
    ya = [r*sin(2*pi*i/100)+b for i in range(0,100)]
    Trace += [go.Scatter3d(x = xa+[None]+xa,
                          y = ya+[None]+xa,
                          z = [H1]*100+[None]+[H2]*100,
                          mode = 'lines',line = dict(color='black', width=1),
                          showlegend = False)]
    th, v = np.meshgrid(np.linspace(0,2*pi,100),np.linspace(H2,H1,50))
    Trace+= [go.Surface(x = r*np.cos(th),y = r*np.sin(th),z = v,
                        opacity=0.2,showlegend = False,showscale=False)]
    
    fig = go.Figure(Trace)
    fig.update_layout(scene=dict(xaxis=dict(range=BB[0]),yaxis=dict(range=BB[1]),
                                 zaxis=dict(range=BB[2])),
                      scene_aspectmode='cube',showlegend=False)
    fig.layout.scene.camera.projection.type = "orthographic"
    fig.show()
        
# ============================================================================
# [ Number of Source Point, Initial Energy, Geometry Parameter(Finite Cyl),
#   Particle type, Geometry Check Function, Cross Section, Density ]
def Volumetric_Source(N,E,P,L,particle,CS,rho):
    Collection  = []
    H1,H2,r,a,b = P
    while len(Collection)<N:
        pt,pr = [random()*2*pi,random()*r]
        x,y,z = [pr*cos(pt),pr*sin(pt),random()*(H2-H1)+H1]
        Check = Finite_Cylinder(P,[x,y,z])
        
        if Check <= 0:
            if particle == 'b+' or particle == 'b-':
                l   = -L*log(1-random())
            elif particle == 'g':
                l   = 1/(Interpolate(CS, E, 'Total wo/ Coherent')*rho)
            Collection += [[{'x':x,'y':y,'z':z,'E':E,'l':l,'C':'initial','p':particle}]]
    return Collection

# ============================================================================
def Beta_Interaction(Collection,n,P,CS,rho,L,**kwargs):
    ps = random()*2*pi
    pt = random()*2*pi
    H1,H2,r,a,b = P
    l           = Collection[n][-1]["l"]
    x,y,z= [Collection[n][-1]["x"]+Collection[n][-1]["l"]*sin(ps)*cos(pt),
            Collection[n][-1]["y"]+Collection[n][-1]["l"]*sin(ps)*sin(pt),
            Collection[n][-1]["z"]+Collection[n][-1]["l"]*cos(ps)]
    if Finite_Cylinder(P,[x,y,z]) <= 0:
        if Collection[n][-1]["p"] == 'b+':
            com = 'annihilation'
            Collection = Positron_Annihilation(Collection,n,P,CS,rho)
        else:
            com = 'stable'
    elif (z >= H2 and (x-a)**2+(y-b)**2-r**2 < 0):
        com = 'detected'
    else:
        com = 'leak'
    Collection[n] += [{'x':x,'y':y,'z':z,'E':0,'l':l,'C':com,
                       'p':Collection[n][-1]['p']}]
    return Collection

# ============================================================================
def Positron_Annihilation(Collection,n,P,CS,rho):
    ps          = random()*2*pi
    pt          = random()*2*pi
    H1,H2,r,a,b = P
    l1          = 1/(Interpolate(CS, 0.511, 'Total wo/ Coherent')*rho)    
    x1,y1,z1,e1 = [Collection[n][-1]['x'],
                   Collection[n][-1]['y'],
                   Collection[n][-1]['z'],0.511]
    
    x2,y2,z2 = [x1 + l1*sin(ps)*cos(pt),
                y1 + l1*sin(ps)*sin(pt),
                z1 + l1*cos(ps)]
    if Finite_Cylinder(P, [x2,y2,z2]) <= 0:
        com2  = 'Travel'
    elif (z2 >= H2 and (x2-a)**2+(y2-b)**2-r**2 < 0):
        com2  = 'detected'
    else:
        com2  = 'leak'

    x3,y3,z3 = [x1 - l1*sin(ps)*cos(pt),
                y1 - l1*sin(ps)*sin(pt),
                z1 - l1*cos(ps)]
    if Finite_Cylinder(P, [x3,y3,z3]) <= 0:
        com3  = 'Travel'
    elif (z3 >= H2 and (x3-a)**2+(y3-b)**2-r**2 < 0):
        com3  = 'detected'
    else:
        com3  = 'leak'
        
    Collection += [[{'x':x1,'y':y1,'z':z1,'E':e1,'l':l1,'C':'annihilation','p':'g'},
                    {'x':x2,'y':y2,'z':z2,'E':0,'l':0,'C':com2,'p':'g'}],
                   [{'x':x1,'y':y1,'z':z1,'E':e1,'l':l1,'C':'annihilation','p':'g'},
                    {'x':x3,'y':y3,'z':z3,'E':0,'l':0,'C':com3,'p':'g'}]]
    return Collection

# ============================================================================
def Gamma_Interaction(Collection,n,P,CS,rho,L):
    stop = False
    H1,H2,r,a,b = P
    while stop==False:
        ps = random()*2*pi
        pt = random()*2*pi
        if len(Collection[n]) == 1:
            e = Collection[n][-1]["E"]
        else:
            e = Compton(Collection, n, ps, pt)
            Collection[n][-1]['E'] = e
            Collection[n][-1]['l'] = 1/(Interpolate(CS,e,'Total wo/ Coherent')*rho)
        x,y,z= [Collection[n][-1]["x"]+Collection[n][-1]["l"]*sin(ps)*cos(pt),
                Collection[n][-1]["y"]+Collection[n][-1]["l"]*sin(ps)*sin(pt),
                Collection[n][-1]["z"]+Collection[n][-1]["l"]*cos(ps)]        
        if Finite_Cylinder(P,[x,y,z]) <= 0:
            CS_To = Interpolate(CS, e, 'Total wo/ Coherent')
            CS_PE = Interpolate(CS, e, 'Photoelectric')
            CS_IC = Interpolate(CS, e, 'Incoherent')
            CS_PP = Interpolate(CS, e, 'Nuclear Pair Production')
            chc   = choices(['Incoherent','Photoelectric','Nuclear Pair Production'],
                            weights=[CS_IC/CS_To,CS_PE/CS_To,CS_PP/CS_To],k=1)[0]
            Collection[n] += [{'x':x,'y':y,'z':z,'E':0,'l':0,'C':chc,'p':'g'}]
            if   chc == 'Photoelectric':
                stop = True
            elif chc == 'Nuclear Pair Production':
                Collection = Pair_Production(Collection,n,P,CS,rho,L)
                stop = True
            else:
                Collection[n][-1]['C'] = 'Travel'
        else:
            if (z >= H2 and (x-a)**2+(y-b)**2-r**2 < 0):
                com = 'detected'
            else:
                com = 'leak'
            Collection[n] += [{'x':x,'y':y,'z':z,'E':0,'l':0,'C':com,'p':'g'}]
            stop= True
    return Collection

# ============================================================================
def Compton(Collection,n,ps,pt):
    x1,y1,z1 = Collection[n][-2]['x'],Collection[n][-2]['y'],Collection[n][-2]['z']
    x2,y2,z2 = Collection[n][-1]['x'],Collection[n][-1]['y'],Collection[n][-1]['z']
    
    v1 = [x2-x1,y2-y1,z2-z1]
    m1 = sum([i**2 for i in v1])**0.5
    v2 = [sin(ps)*cos(pt),sin(ps)*sin(pt),cos(ps)]; 
    m2 = sum([i**2 for i in v2])**0.5
    E  = Collection[n][-2]['E']
    
    costh = sum([i*j for i,j in zip(v1,v2)])/(m1*m2)
    th    = acos(costh)
    e     = E/(1+E/(0.511**2)*(1-cos(th)))
    return e

# ============================================================================
def Pair_Production(Collection,n,P,CS,rho,L):
    ps          = random()*2*pi
    pt          = random()*2*pi
    H1,H2,r,a,b = P
    l1          = L    
    x1,y1,z1,e1 = [Collection[n][-1]['x'],
                   Collection[n][-1]['y'],
                   Collection[n][-1]['z'],
                   Collection[n][-2]['E']]
    
    x2,y2,z2 = [x1 + l1*sin(ps)*cos(pt),
                y1 + l1*sin(ps)*sin(pt),
                z1 + l1*cos(ps)]
    if Finite_Cylinder(P, [x2,y2,z2]) <= 0:
        com2  = 'stable'
    elif (z2 >= H2 and (x2-a)**2+(y2-b)**2-r**2 < 0):
        com2  = 'detected'
    else:
        com2  = 'leak'

    x3,y3,z3 = [x1 - l1*sin(ps)*cos(pt),
                y1 - l1*sin(ps)*sin(pt),
                z1 - l1*cos(ps)]
    if Finite_Cylinder(P, [x3,y3,z3]) <= 0:
        com3  = 'annihilation'
    elif (z3 >= H2 and (x3-a)**2+(y3-b)**2-r**2 < 0):
        com3  = 'detected'
    else:
        com3  = 'leak'
        
    Collection += [[{'x':x1,'y':y1,'z':z1,'E':e1/2,'l':l1,'C':'annihilation','p':'b-'},
                    {'x':x2,'y':y2,'z':z2,'E':0,'l':0,'C':com2,'p':'b-'}],
                   [{'x':x1,'y':y1,'z':z1,'E':e1/2,'l':l1,'C':'annihilation','p':'b+'},
                    {'x':x3,'y':y3,'z':z3,'E':0,'l':0,'C':com3,'p':'b+'}]]
    if com3 == 'annihilation':
        Collection = Positron_Annihilation(Collection, -1, P, CS, rho)
    return Collection

# MAIN PPROGRAM
# ============================================================================
def Main_Program(r,L,h,E,N,rho,particle,dirpath,plot = False):
    #start = time.time()
    #units in mm
    P = [0,h,r,0,0]

    CS  = read_gamma_CS(dirpath)

    Collection   = Volumetric_Source(N, E, P, L, particle, CS, rho)

    n=0; Trace = []
    while n<len(Collection):
        if (Collection[n][-1]['C'] == 'Travel' or
            Collection[n][-1]['C'] == 'initial'):
            if Collection[n][-1]['p'] == 'b+' or Collection[n][-1]['p'] == 'b-':
                Collection = Beta_Interaction(Collection, n, P, CS, rho, L)
            elif Collection[n][-1]['p'] == 'g':
                Collection = Gamma_Interaction(Collection, n, P, CS, rho, L)
        Trace += print_out(Collection,n)
        n+=1

    #print('Time Elapsed : ',time.time()-start,' s\n for simulating ',len(Collection),' particles')
    if plot == True:
        Plot(P,[[-r,r],[-r,r],[0,h]],Trace = Trace)
    return Collection

# ============================================================================
start = time.time()
#units in cm
L = 0.2              # Beta Mean Range [2]
H = 0.5              # Cylinder Height [0.5]
R = 0.9              # Radius of Cylinder [0.9]
E = 2                # Energy of beta/Gamma [2]
N = 1e5
particle = 'b+'
rho      = 1         # Medium density [g/cm3]
dirpath  = 'H2O_CS'
n = 50

# Check And Visualize
'''
Stop=False
while Stop==False:
    Collection = Main_Program(R,L,H,E,N,rho,particle,dirpath)
    for part in Collection:
        if part[-1]['C'] == 'detected':
            Stop=True
            
Trace = []
for n in range(0,len(Collection)):
    Trace += print_out(Collection,n)
P = [0,H,R,0,0]
Plot(P,[[-R,R],[-R,R],[0,H]],Trace = Trace)
'''
# Radius Variation and Height Variation
# ============================================================================

g_detected = [0]*n; g_energy=[]; err_g=[0]*n
b_detected = [0]*n; b_energy=[]; err_b=[0]*n
g_produced = [0]*n

#for m,h in enumerate([H*i for i in range(1,n+1)]):
#for m,r in enumerate([R*i for i in range(1,n+1)]):
for m,c in enumerate([0]*n):
#for m,d in enumerate([10**i for i in range(1,n+1)]):
    #Collection = Main_Program(R,L,h,E,N,rho,particle,dirpath)
    #Collection = Main_Program(r,L,H,E,N,rho,particle,dirpath)
    Collection = Main_Program(R,L,H,E,N,rho,particle,dirpath)
    #Collection = Main_Program(R,L,H,E,d,rho,particle,dirpath)
    for part in Collection:
        if part[-1]['C'] == 'detected':
            if part[-1]['p'] == 'g':
                g_detected[m] += 1
                g_energy      += [part[-2]['E']]
            elif part[-1]['p'] == particle :
                b_detected[m] += 1
                b_energy      += [part[-2]['E']]
        else:
            if part[-1]['p'] == 'g':
                g_produced[m] += 1
    
    g_detected[m] *= 1/N
    g_produced[m] *= 1/N
    b_detected[m] *= 1/N
    
plt.figure(figsize=[5,4],dpi=300)

#plt.plot([R*i for i in range(1,n+1)],g_detected,label='gamma detected/source',c='b')
#plt.plot([R*i for i in range(1,n+1)],g_produced,label='gamma produced/source',c='g')
#plt.plot([R*i for i in range(1,n+1)],b_detected,label=particle+' detected/source',c='r')

#plt.plot([H*i for i in range(1,n+1)],g_detected,label='gamma detected/source')
#plt.plot([H*i for i in range(1,n+1)],g_produced,label='gamma produced/source')
#plt.plot([H*i for i in range(1,n+1)],b_detected,label=particle+'detected/source')

plt.plot([i for i in range(1,n+1)],g_detected,label='gamma detected/source')
plt.plot([i for i in range(1,n+1)],g_produced,label='gamma produced/source')
plt.plot([i for i in range(1,n+1)],b_detected,label=particle+'detected/source')

#plt.plot([10**(i) for i in range(1,n+1)],g_detected,label='gamma detected/source')
#plt.plot([10**(i) for i in range(1,n+1)],g_produced,label='gamma produced/source')
#plt.plot([10**(i) for i in range(1,n+1)],b_detected,label=particle+'detected/source')

plt.ylabel('efficiency')
#plt.xlabel('height [cm]')
#plt.xlabel('radius [cm]')
plt.xlabel('trial')
#plt.xlabel('Particle Simulated')

#plt.xlim([H,n*H])
#plt.xlim([R,n*R])
plt.xlim(1,n+1)
#plt.xlim(10,10**n)

plt.ylim([0,max(g_detected+g_produced+b_detected)])
plt.legend()
plt.grid()
plt.show()

print('Time Elapsed : ',time.time()-start,' s\n for simulating ',len(Collection),' particles')
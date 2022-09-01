# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 11:19:49 2019

@author: andrew.ferguson
"""
import numpy as np
import csv
import os

path = os.getcwd()
basepath = os.path.dirname(path)

##Builds up the dictionaries of thermal parameters

k = {}
c = {}
rho = {}

filename = basepath+"/src/thermal_properties.csv"

with open(filename,'r') as f:
    csvreader = csv.DictReader(f)
    for row in csvreader:
        k[row['material']]= float(row['thermal_conductivity'])
        c[row['material']]= float(row['specific_heat'])
        rho[row['material']]= float(row['density'])

##
          
def test2d():    
    nx=101
    ny=101
    dx = 3e-6
    dy = 1e-6

    
    x_lim = (nx-1)*dx/2
    y_lim = (ny-1)*dy/2
    
    x = np.linspace(-x_lim,x_lim,nx)
    y = np.linspace(-y_lim,y_lim,ny)
    
    k_in = np.zeros((nx,ny))
    k_in[:,:] =k['water']
    
    Q_in = np.zeros((nx,ny))
    Q_in_val = 1
    Q_in[:,:] = Q_in_val
    
    v_in = np.zeros((nx,ny))
    c_in = np.zeros((nx,ny))
    rho_in = np.zeros((nx,ny))
    
    c_in[:,:]=c['water']
    rho_in[:,:]=rho['water']
    
    return [k_in,c_in,rho_in,Q_in,v_in,x,y,dx,dy,nx,ny]

def single_wire_membrane(v = 0e-3,wet=True):
    nx=401
    ny=151
    dx = 2.5e-6
    dy = 1e-6
    
    x_lim = (nx-1)*dx/2
    y_lim = (ny-1)*dy/2
    
    x = np.linspace(-x_lim,x_lim,nx)
    y = np.linspace(-y_lim,y_lim,ny)
    
    k_in = np.zeros((nx,ny))
    c_in = np.zeros((nx,ny))
    rho_in = np.zeros((nx,ny))
    Q_in = np.zeros((nx,ny))       
    
    v_max=-v
    
    if wet == True:
        k_in[:,50:101 ] = k['water']
        c_in[:,50:101] = c['water']
        rho_in[:,50:101] = rho['water']
    elif wet == False:
        k_in[:,50:101 ] = k['air']
        c_in[:,50:101] = c['air']
        rho_in[:,50:101] = rho['air']
    
    ##put some glass on top
    k_in[:,101:] = k['silicon_dioxide']
    c_in[:,101:] = c['silicon_dioxide']
    rho_in[:,101:] = rho['silicon_dioxide']
    
    
    ##silicon dioxide underneath
    k_in[:,48:50] = k['silicon_dioxide']
    c_in[:,48:50] = c['silicon_dioxide']
    rho_in[:,48:50] = rho['silicon_dioxide']
   
     
    ##silicon dioxide underneath
    k_in[:,0:48] = k['silicon_dioxide']
    c_in[:,0:48] = c['silicon_dioxide']
    rho_in[:,0:48] = rho['silicon_dioxide']

    ##and some air
    k_in[200-50:200+50,0:48] = k['air']
    c_in[200-50:200+50,0:48] = c['air']
    rho_in[200-50:200+50,0:48] = rho['air']

    ##and some silicon
    k_in[0:200-50,0:48] = k['silicon']
    c_in[0:200-50,0:48] = c['silicon']
    rho_in[0:200-50,0:48] = rho['silicon']

    ##and some silicon
    k_in[200+50:,0:48] = k['silicon']
    c_in[200+50:,0:48] = c['silicon']
    rho_in[200+50,0:48] = rho['silicon']
        
    heats = 8
    
    Q_in[(200-5):(200+6),49:50] = heats
      
  
    v_in = np.zeros((nx,ny))
    for yval in range(50,101):
        v_in[:,yval] = v_max - 4*v_max/50**2*(yval-75)**2
        
    return [k_in,c_in,rho_in,Q_in,v_in,x,y,dx,dy,nx,ny]


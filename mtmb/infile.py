#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 18:04:57 2021

@author: jiedeng
"""
#import os  
#import numpy as np
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--plumed","-pf",type=str,default='plumed.dat',help="template plumed.dat file")
parser.add_argument("--lammps","-lf",type=str,default='in.lammps',help="template in.lammps file")
parser.add_argument('--energy',"-e", type=float,nargs="+",help="energy in ev, order does not matter")
parser.add_argument('--volume',"-v", type=float,nargs="+",help="volume in A3, order does not matter")
parser.add_argument('--pressure',"-p", type=float,nargs="+",help="pressure in gpa, order does not matter")
parser.add_argument('--temperature',"-t", type=float,nargs="+",help="temperature in k, order does not matter")
parser.add_argument('--starting_from_melt',"-m", default=False, action='store_true',help="melt the system first ")

args   = parser.parse_args()

pd=open(args.plumed)
pd_new = open('plumed.dat','w')
max_e =  max(args.energy)
min_e =  min(args.energy)
max_t =  max(args.temperature)
min_t = min(args.temperature)
max_p =  max(args.pressure)
min_p =  min(args.pressure)
max_v =  max(args.volume)
min_v =  min(args.volume)

e_factor = 96.487
p_factor = 1000000000*6.02E-10*1000
for line in pd:

    if 'bf1: BF_LEGENDRE' in line:
        line  = 'bf1: BF_LEGENDRE ORDER=8 MINIMUM={0} MAXIMUM={1}\n'.format(round(min_e*e_factor), round(max_e*e_factor))
    if 'bf2: BF_LEGENDRE' in line:
        line = 'bf2: BF_LEGENDRE ORDER=8 MINIMUM={0} MAXIMUM={1}\n'.format(min_v, max_v)
    if 'MIN_TEMP' in line:
        line = ' MIN_TEMP={0}\n'.format(min_t)
    if 'MAX_TEMP' in line:
        line = ' MAX_TEMP={0}\n'.format(max_t)
    if 'MIN_PRESSURE' in line:
        line = ' MIN_PRESSURE={0}\n'.format(round(min_p*p_factor))
    if 'MAX_PRESSURE' in line:
        line = ' MAX_PRESSURE={0}\n'.format(round(max_p*p_factor))
    if ' PRESSURE=' in line:
        line = ' PRESSURE={0}\n'.format(round((max_p+min_p)/2))
    if ' TEMP=' in line:
        line = ' TEMP={0}\n'.format(round((min_t+max_t)/2))
    pd_new.writelines(line)
pd_new.close()

print('plumed.dat file generated')
lammps = open(args.lammps)
lammps_new = open('in.lammps','w')

for line in lammps:
    if 'variable        temperature equal' in line:
        line  = 'variable        temperature equal {0}\n'.format((min_t+max_t)/2)
    if 'variable        pressure equal' in line:
        line = 'variable        pressure equal {0}\n'.format((min_p+max_p)/2*10000)
    if args.starting_from_melt:
        if 'velocity       all create ' in line:
            tin = min_t + 3000
            line  = 'velocity        all create {0} 23456\n'.format(tin)
        if 'fix             1 all npt  temp ${temperature} ${temperature} ${tempDamp} iso' in line:
            line = 'fix             1 all npt  temp {2} {3} 0.01 iso {0} {1} 0.1\n'.format(min_p*10000,min_p*10000, tin,tin)
    lammps_new.writelines(line)

pd_new.close()
print('in.lammps file generated')

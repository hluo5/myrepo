#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  9 08:35:40 2021

@author: jiedeng
"""

from lammps_logfile import File
import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description="Plot contents from lammps log files")
parser.add_argument("--input_file",'-i', type=str,default='log.lammps', help="Lammps log file containing thermo output from lammps simulation.")
parser.add_argument("--ptarget",'-p', type=float, help="target p value, in GPa")
parser.add_argument("--ttarget",'-t', type=float, help="target t value")
parser.add_argument("--plot",'-pl', default=False, action="store_true",help="plot the result?")
parser.add_argument("--exe",'-e', default=True, action="store_false",help="plot the result?")
parser.add_argument("--infile",'-if', default='./infile.py',help="infile.py to modify the plumed.dat and in.lammps")
parser.add_argument("--plumed",'-pf', default='./plumed.dat',help="plumed.dat template file")
parser.add_argument("--lammps",'-lf', default='./mtmb.lammps',help="in.lammps template")
args = parser.parse_args()

try:
    log = File(args.input_file)
except:
    from subprocess import call
    call("cp {0} {1}".format(args.input_file, args.input_file+'tmp'), shell=True) # do not change in the original file, better for checking on running sinulation
    call("sed -i 's/style restartinfo set but has//' {0}".format(args.input_file+'tmp'), shell=True)
    log = File(args.input_file+'tmp')
    call("rm {0}".format(args.input_file+'tmp'), shell=True)


y = ['Press','Temp','PotEng','Volume']
ys  = [log.get(y,run_num=-1) for y in y]


def check(dat):
    for ele in dat:
        if str(ele).replace('.','').replace('e','').replace('d','').replace('+','').replace('-','').isdigit(): # data may be in int or float, may contain . , e , d , + , -
            pass
        else:
            return False
    return True

def select(dat):
    selected_idx = []
    for i  in range(len(dat)):
        if str(dat[i]).isdigit():
            selected_idx.append(i)
    return selected_idx


Step = log.get('Step',run_num=-1)
if not check(Step):
    print('**data messed up**')
    print('    Step col is:', Step[:10])
    print('**data messed up**')
    selected_idx = select(Step)
    
#    x = (x[selected_idx]).astype(float)
    ys = [(y[selected_idx]).astype(float) for y in ys]
    print('**Fixed**')

def misfit(P,T,Ptarget,Ttarget):
    mis_P  = ((P-Ptarget)/Ptarget)**2
    mis_T  = ((T-Ttarget)/Ttarget)**2
    mis = (mis_P + mis_T)
    return mis

mis=misfit(ys[0],ys[1],args.ptarget*1e4,args.ttarget)

#plt.plot(mis)
#plt.plot(ys[0])
#plt.hist(mis,bins=100)


#plt.plot(mis[:100])
ys.append(mis)
ys = np.array(ys)
ys_sorted=ys[:,ys[-1,:].argsort()]

if args.plot:
    fig,ax = plt.subplots(4,figsize=(8,10))
    ax[0].plot(ys_sorted[-1,:100],ys_sorted[0,:100],label=y[0])
    ax[0].grid()
    ax[0].legend()
    
    ax[1].plot(ys_sorted[-1,:100],ys_sorted[1,:100],label=y[1])
    ax[1].grid()
    ax[1].legend()
    
    ax[2].plot(ys_sorted[-1,:100],ys_sorted[2,:100],label=y[2])
    ax[2].grid()
    ax[2].legend()
    
    ax[3].plot(ys_sorted[-1,:100],ys_sorted[3,:100],label=y[3])
    ax[2].grid()
    ax[3].legend()
    plt.show()

# Calculate the average value for each row
ys_sorted_avg = np.mean(ys_sorted, axis=1)
# Print the average values instead of the first value
print("**'step','Press','Temp','PotEng','Volume', 'mis'")
print(ys_sorted_avg)



# Calculate the average of each row (axis 1 means calculating the mean across columns)
ys_sorted_avg_2 = np.mean(ys_sorted[2, :])  # Average of the 3rd row
ys_sorted_avg_3 = np.mean(ys_sorted[3, :])  # Average of the 4th row
# Now use these averages in the calculation for emin, emax, vmin, vmax
emin, emax, vmin, vmax = ys_sorted_avg_2 - 50, ys_sorted_avg_2 + 50, ys_sorted_avg_3 - 50, ys_sorted_avg_3 + 50
pmin,pmax,tmin,tmax = args.ptarget - 15, args.ptarget +15,  args.ttarget -1500, args.ttarget +1500
#print('v',vmin,vmax)
#print('ys_sorted[3,-1]',ys_sorted[3,-1])

#string = "python ../../infile.py -pf ../../plumed.dat -lf ../../in.lammps -e {0} {1} -v  {2} {3} -p {4} {5} -t {6} {7} -m".format(emin,emax,vmin,vmax,pmin,pmax,tmin,tmax)
string = "python {8} -pf {9} -lf {10} -e {0} {1} -v  {2} {3} -p {4} {5} -t {6} {7}".format(emin,emax,vmin,vmax,pmin,pmax,tmin,tmax, args.infile,args.plumed, args.lammps)

# Convert from eV to kJ/mol
conversion_factor = 96.485
emin_kJmol = emin * conversion_factor
emax_kJmol = emax * conversion_factor
# Print the values in kJ/mol
print(f"emin in kJ/mol: {emin_kJmol}")
print(f"emax in kJ/mol: {emax_kJmol}")


print(string)
print('=='*20)
print('change nn in in.lammps if necessary!')

if args.exe:
    print("execute the above command")
    from subprocess import call
    call(string, shell=True)

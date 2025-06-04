#!/usr/bin/env python3

from scipy.optimize import fsolve
import numpy as np

#correction_method 3
# Define the system of equations
def equations(variables):
	x1, x2, y1, y2, m1, m2, n1, n2 = variables
        
	xo = x1/(x1+y1+m1+n1)
	xmg = y1/(x1+y1+m1+n1)
	xsi = m1/(x1+y1+m1+n1)
	xfe = n1/(x1+y1+m1+n1)
	xmgo = y2/(y2+m2+n2)
	xsio2 = m2/(y2+m2+n2)
	xfeo = n2/(y2+m2+n2)

	a1 = -4.19473809
	b1 = 0.0
	c1 = 59.4760141
	a2 = 0.0
	b2 = 0.0
	c2 = 0.0
	a3 = -8.86496022
	b3 = 0.0
	c3 = 78.4827585
	a4 = -8.85508053
	b4 = 0.0
	c4 = 122.987517
	hh = -2.97544334
	ho = 0.0
	hmg = 0.0
	hsi = 3.11554783
	oo =  5.46470737
	omg = 0.0
	osi = -14.6661946
	mgmg = -580.932981
	mgsi = 31.5172795
	sisi = 23.9802402
	
	a21 = a2
	b21 = b2
	c21 = c2
	a31 = a3
	b31 = b3
	c31 = c3
	a41 = a4
	b41 = b4
	c41 = c4
	oo1 = oo
	omg1 = omg
	osi1 = osi
	mgmg1 = mgmg
	mgsi1 = mgsi
	sisi1 = sisi
	
	pres = 40
	temp = 3750

	lnGamma_O = -oo1*1873/temp*np.log(1-xo) - (omg1*1873/temp*xmg*(1+np.log(1-xmg)/xmg-1/(1-xo)) + osi1*1873/temp*xsi*(1+np.log(1-xsi)/xsi-1/(1-xo))) + (omg1*1873/temp*xmg**2*xo*(1/(1-xo)+1/(1-xmg)+xo/(2*(1-xo)**2)-1) + osi1*1873/temp*xsi**2*xo*(1/(1-xo)+1/(1-xsi)+xo/(2*(1-xo)**2)-1))
 
# Define your ten equations here
	eq1 = x1+x2-34000
	eq2 = y1+y2-12000
	eq3 = m1+m2-10714
	eq4 = n1+n2-11000
	eq5 = y2+2*m2+n2-x2-0
	eq6 = a21 + b21/temp + c21*pres/temp - lnGamma_O - np.log(xfe*xo/xfeo) - 0.0
	eq7 = a31 + b31/temp + c31*pres/temp - lnGamma_O - (-mgmg1*1873/temp*np.log(1-xmg) - (omg1*1873/temp*xo*(1+np.log(1-xo)/xo-1/(1-xmg)) + mgsi1*1873/temp*xsi*(1+np.log(1-xsi)/xsi-1/(1-xmg))) + (omg1*1873/temp*xo**2*xmg*(1/(1-xmg)+1/(1-xo)+xmg/(2*(1-xmg)**2)-1) + mgsi1*1873/temp*xsi**2*xmg*(1/(1-xmg)+1/(1-xsi)+xmg/(2*(1-xmg)**2)-1))) - np.log(xmg*xo/xmgo) - 0.0
	eq8 = a41 + b41/temp + c41*pres/temp - 2*lnGamma_O - (-sisi1*1873/temp*np.log(1-xsi) - (osi1*1873/temp*xo*(1+np.log(1-xo)/xo-1/(1-xsi)) + mgsi1*1873/temp*xmg*(1+np.log(1-xmg)/xmg-1/(1-xsi))) + (osi1*1873/temp*xo**2*xsi*(1/(1-xsi)+1/(1-xo)+xsi/(2*(1-xsi)**2)-1) + mgsi1*1873/temp*xmg**2*xsi*(1/(1-xsi)+1/(1-xmg)+xsi/(2*(1-xsi)**2)-1))) - np.log(xsi*xo**2/xsio2) - 0.0

	return [eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8]

# Initial guess
initial_guess = np.array([1086, 32913, 19, 11980, 858, 9855, 9778, 1221])

# Solve the system of equations
solution = fsolve(equations, initial_guess)

print("Solution:", solution)



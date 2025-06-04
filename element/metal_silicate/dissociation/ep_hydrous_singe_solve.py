#!/usr/bin/env python3

from scipy.optimize import fsolve
import numpy as np

#correction_method 3
# Define the system of equations
def equations(variables):
	p1, p2, x1, x2, y1, y2, m1, m2, n1, n2 = variables
        
	xh = p1/(p1+x1+y1+m1+n1)
	xo = x1/(p1+x1+y1+m1+n1)
	xmg = y1/(p1+x1+y1+m1+n1)
	xsi = m1/(p1+x1+y1+m1+n1)
	xfe = n1/(p1+x1+y1+m1+n1)
	xh2o = p2/2/(p2/2+y2+m2+n2)
	xmgo = y2/(p2/2+y2+m2+n2)
	xsio2 = m2/(p2/2+y2+m2+n2)
	xfeo = n2/(p2/2+y2+m2+n2)

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
	
	pres = 54
	temp = 3350

	lnGamma_O = -oo*1873/temp*np.log(1-xo) - (ho*1873/temp*xh*(1+np.log(1-xh)/xh-1/(1-xo)) + omg*1873/temp*xmg*(1+np.log(1-xmg)/xmg-1/(1-xo)) + osi*1873/temp*xsi*(1+np.log(1-xsi)/xsi-1/(1-xo))) + (ho*1873/temp*xh**2*xo*(1/(1-xo)+1/(1-xh)+xo/(2*(1-xo)**2)-1) + omg*1873/temp*xmg**2*xo*(1/(1-xo)+1/(1-xmg)+xo/(2*(1-xo)**2)-1) + osi*1873/temp*xsi**2*xo*(1/(1-xo)+1/(1-xsi)+xo/(2*(1-xo)**2)-1))
    
# Define your ten equations here
	eq1 = p1+p2-5000
	eq2 = x1+x2-36500
	eq3 = y1+y2-12000
	eq4 = m1+m2-10714
	eq5 = n1+n2-11000
	eq6 = p2/2+y2+2*m2+n2-x2-0
	eq7 = a1 + b1/temp + c1*pres/temp - lnGamma_O - 2*(-hh*1873/temp*np.log(1-xh) - (ho*1873/temp*xo*(1+np.log(1-xo)/xo-1/(1-xh)) + hmg*1873/temp*xmg*(1+np.log(1-xmg)/xmg-1/(1-xh)) + hsi*1873/temp*xsi*(1+np.log(1-xsi)/xsi-1/(1-xsi))) + (ho*1873/temp*xo**2*xh*(1/(1-xh)+1/(1-xo)+xh/(2*(1-xh)**2)-1) + hmg*1873/temp*xmg**2*xh*(1/(1-xh)+1/(1-xmg)+xh/(2*(1-xh)**2)-1) + hsi*1873/temp*xsi**2*xh*(1/(1-xh)+1/(1-xsi)+xh/(2*(1-xh)**2)-1))) - np.log(xh**2*xo/xh2o) - 0.0
	eq8 = a2 + b2/temp + c2*pres/temp - lnGamma_O - np.log(xfe*xo/xfeo) - 0.0
	eq9 = a3 + b3/temp + c3*pres/temp - lnGamma_O - (-mgmg*1873/temp*np.log(1-xmg) - (hmg*1873/temp*xh*(1+np.log(1-xh)/xh-1/(1-xmg)) + omg*1873/temp*xo*(1+np.log(1-xo)/xo-1/(1-xmg)) + mgsi*1873/temp*xsi*(1+np.log(1-xsi)/xsi-1/(1-xmg))) + (hmg*1873/temp*xh**2*xmg*(1/(1-xmg)+1/(1-xh)+xmg/(2*(1-xmg)**2)-1) + omg*1873/temp*xo**2*xmg*(1/(1-xmg)+1/(1-xo)+xmg/(2*(1-xmg)**2)-1) + mgsi*1873/temp*xsi**2*xmg*(1/(1-xmg)+1/(1-xsi)+xmg/(2*(1-xmg)**2)-1))) - np.log(xmg*xo/xmgo) - 0.0
	eq10 = a4 + b4/temp + c4*pres/temp - 2*lnGamma_O - (-sisi*1873/temp*np.log(1-xsi) - (hsi*1873/temp*xh*(1+np.log(1-xh)/xh-1/(1-xsi)) + osi*1873/temp*xo*(1+np.log(1-xo)/xo-1/(1-xsi)) + mgsi*1873/temp*xmg*(1+np.log(1-xmg)/xmg-1/(1-xsi))) + (hsi*1873/temp*xh**2*xsi*(1/(1-xsi)+1/(1-xh)+xsi/(2*(1-xsi)**2)-1) + osi*1873/temp*xo**2*xsi*(1/(1-xsi)+1/(1-xo)+xsi/(2*(1-xsi)**2)-1) + mgsi*1873/temp*xmg**2*xsi*(1/(1-xsi)+1/(1-xmg)+xsi/(2*(1-xsi)**2)-1))) - np.log(xsi*xo**2/xsio2) - 0.0

	return [eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10]

# Initial guess
initial_guess = np.array([1522, 1477, 968, 34031, 50, 11919, 502, 10411, 9849, 1250])

# Solve the system of equations
solution = fsolve(equations, initial_guess)

print("Solution:", solution)



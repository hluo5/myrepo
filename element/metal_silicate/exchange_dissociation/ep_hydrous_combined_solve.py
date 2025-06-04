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

	a1 = 0.0
	b1 = -10503.1021
	c1 = 0.0
	a2 = 0.0
	b2 = -3871.90415
	c2 = 0.0
	a3 = -10.6391133
	b3 = 0.0
	c3 = 79.6122262
	a4 = -3.44034886
	b4 = -8805.31458
	c4 = 0.0
	hh = 0.0
	ho = -3.67745982
	hmg = 0.0
	hsi = 11.3747330
	oo = 0.0
	omg = -25.6828142
	osi = -9.73803053
	mgmg = 0.0
	mgsi = 0.0
	sisi = 0.0
	
	pres = 50
	temp = 3500

	lnGamma_O = -oo*1873/temp*np.log(1-xo) - (ho*1873/temp*xh*(1+np.log(1-xh)/xh-1/(1-xo)) + omg*1873/temp*xmg*(1+np.log(1-xmg)/xmg-1/(1-xo)) + osi*1873/temp*xsi*(1+np.log(1-xsi)/xsi-1/(1-xo))) + (ho*1873/temp*xh**2*xo*(1/(1-xo)+1/(1-xh)+xo/(2*(1-xo)**2)-1) + omg*1873/temp*xmg**2*xo*(1/(1-xo)+1/(1-xmg)+xo/(2*(1-xo)**2)-1) + osi*1873/temp*xsi**2*xo*(1/(1-xo)+1/(1-xsi)+xo/(2*(1-xo)**2)-1))
    
# Define your ten equations here
	eq1 = p1+p2-8000
	eq2 = x1+x2-38000
	eq3 = y1+y2-12000
	eq4 = m1+m2-11000
	eq5 = n1+n2-11000
	eq6 = p2/2+y2+2*m2+n2-x2-0
	eq7 = a1 + b1/temp + c1*pres/temp - 2*(-hh*1873/temp*np.log(1-xh) - (ho*1873/temp*xo*(1+np.log(1-xo)/xo-1/(1-xh)) + hmg*1873/temp*xmg*(1+np.log(1-xmg)/xmg-1/(1-xh)) + hsi*1873/temp*xsi*(1+np.log(1-xsi)/xsi-1/(1-xsi))) + (ho*1873/temp*xo**2*xh*(1/(1-xh)+1/(1-xo)+xh/(2*(1-xh)**2)-1) + hmg*1873/temp*xmg**2*xh*(1/(1-xh)+1/(1-xmg)+xh/(2*(1-xh)**2)-1) + hsi*1873/temp*xsi**2*xh*(1/(1-xh)+1/(1-xsi)+xh/(2*(1-xh)**2)-1))) - np.log(xh**2*xfeo/xh2o/xfe) - 0.0
	eq8 = a2 + b2/temp + c2*pres/temp - lnGamma_O - np.log(xfe*xo/xfeo) - 0.0
	eq9 = a3 + b3/temp + c3*pres/temp - lnGamma_O - (-mgmg*1873/temp*np.log(1-xmg) - (hmg*1873/temp*xh*(1+np.log(1-xh)/xh-1/(1-xmg)) + omg*1873/temp*xo*(1+np.log(1-xo)/xo-1/(1-xmg)) + mgsi*1873/temp*xsi*(1+np.log(1-xsi)/xsi-1/(1-xmg))) + (hmg*1873/temp*xh**2*xmg*(1/(1-xmg)+1/(1-xh)+xmg/(2*(1-xmg)**2)-1) + omg*1873/temp*xo**2*xmg*(1/(1-xmg)+1/(1-xo)+xmg/(2*(1-xmg)**2)-1) + mgsi*1873/temp*xsi**2*xmg*(1/(1-xmg)+1/(1-xsi)+xmg/(2*(1-xmg)**2)-1))) - np.log(xmg*xo/xmgo) - 0.0
	eq10 = a4 + b4/temp + c4*pres/temp - (-sisi*1873/temp*np.log(1-xsi) - (hsi*1873/temp*xh*(1+np.log(1-xh)/xh-1/(1-xsi)) + osi*1873/temp*xo*(1+np.log(1-xo)/xo-1/(1-xsi)) + mgsi*1873/temp*xmg*(1+np.log(1-xmg)/xmg-1/(1-xsi))) + (hsi*1873/temp*xh**2*xsi*(1/(1-xsi)+1/(1-xh)+xsi/(2*(1-xsi)**2)-1) + osi*1873/temp*xo**2*xsi*(1/(1-xsi)+1/(1-xo)+xsi/(2*(1-xsi)**2)-1) + mgsi*1873/temp*xmg**2*xsi*(1/(1-xsi)+1/(1-xmg)+xsi/(2*(1-xsi)**2)-1))) - np.log(xsi*xfeo**2/xsio2/xfe**2) - 0.0

	return [eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10]

# Initial guess
initial_guess = np.array([1200, 800, 1500, 33500, 500, 11500, 1000, 10000, 8000, 3000])

# Solve the system of equations
solution = fsolve(equations, initial_guess)

print("Solution:", solution)



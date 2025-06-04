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
	
	pres = 50
	temp = 3500

	lnGamma_O = -oo1*1873/temp*np.log(1-xo) - (omg1*1873/temp*xmg*(1+np.log(1-xmg)/xmg-1/(1-xo)) + osi1*1873/temp*xsi*(1+np.log(1-xsi)/xsi-1/(1-xo))) + (omg1*1873/temp*xmg**2*xo*(1/(1-xo)+1/(1-xmg)+xo/(2*(1-xo)**2)-1) + osi1*1873/temp*xsi**2*xo*(1/(1-xo)+1/(1-xsi)+xo/(2*(1-xo)**2)-1))
 
# Define your ten equations here
	eq1 = x1+x2-34000
	eq2 = y1+y2-12000
	eq3 = m1+m2-11000
	eq4 = n1+n2-11000
	eq5 = y2+2*m2+n2-x2-0
	eq6 = a21 + b21/temp + c21*pres/temp - lnGamma_O - np.log(xfe*xo/xfeo) - 0.0
	eq7 = a31 + b31/temp + c31*pres/temp - lnGamma_O - (-mgmg1*1873/temp*np.log(1-xmg) - (omg1*1873/temp*xo*(1+np.log(1-xo)/xo-1/(1-xmg)) + mgsi1*1873/temp*xsi*(1+np.log(1-xsi)/xsi-1/(1-xmg))) + (omg1*1873/temp*xo**2*xmg*(1/(1-xmg)+1/(1-xo)+xmg/(2*(1-xmg)**2)-1) + mgsi1*1873/temp*xsi**2*xmg*(1/(1-xmg)+1/(1-xsi)+xmg/(2*(1-xmg)**2)-1))) - np.log(xmg*xo/xmgo) - 0.0
	eq8 = a41 + b41/temp + c41*pres/temp - (-sisi1*1873/temp*np.log(1-xsi) - (osi1*1873/temp*xo*(1+np.log(1-xo)/xo-1/(1-xsi)) + mgsi1*1873/temp*xmg*(1+np.log(1-xmg)/xmg-1/(1-xsi))) + (osi1*1873/temp*xo**2*xsi*(1/(1-xsi)+1/(1-xo)+xsi/(2*(1-xsi)**2)-1) + mgsi1*1873/temp*xmg**2*xsi*(1/(1-xsi)+1/(1-xmg)+xsi/(2*(1-xsi)**2)-1))) - np.log(xsi*xfeo**2/xsio2/xfe**2) - 0.0

	return [eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8]

# Initial guess
initial_guess = np.array([800, 32200, 30, 11970, 1500, 9500, 9000, 2000])

# Solve the system of equations
solution = fsolve(equations, initial_guess)

print("Solution:", solution)



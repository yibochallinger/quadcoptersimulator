'''calculate thrust and torque using blade element theory and conversation of momentum'''
from constant import *
import numpy as np
import math

def thrust_and_torque(v_z, v2): #motor_omega, v_z, disk plane velocity
	

	rad = 3./4.*blade_length
	#local blade element setting angle
	theta = math.atan(pitch/2.0/math.pi/rad) 
	#initial guess of inflow and swirl factor
	a = .1
	b = .01
	finished = False
	total = 1
	chord = blade_width
	while(finished is False):
		#axial velocity
		v0 = v_z*(1+a)
		phi = math.atan2(abs(v0),abs(v2))  
		if (v_z>0):
			alpha = theta-phi
		if (v_z<0):
			alpha = theta+phi
		cl = 6.2*alpha
		cd = .008-.003*cl+.01*cl**2
		#local velocity at blade
		vlocal = math.sqrt(v0**2+v2**2)
		#thrust grading per blade 
		DtDr = 0.5*rho*vlocal**2.0*chord*(cl*math.cos(phi)-cd*math.sin(phi))

		#torque grading per blade
		DqDr = 0.5*rho*vlocal**2.0*chord*rad*(cd*math.cos(phi)+cl*math.sin(phi))

		#momentum check on inflow and swirl factors
		tem1 = DtDr/(4.0*math.pi*rad*rho*v_z**2*(1+a))
		anew = 0.5*(a+tem1)

		#check for convergence
		if(abs(anew-a)<1.0e-5):
			finished = True
		a = anew

		total =+1
		#check to see if iteration stuck
		if (total>500):
			finished=True
	if (v2>0):
		thrust = DtDr*blade_length
		torque = DqDr*blade_length
	else:
		thrust = -DtDr*blade_length
		torque = -DqDr*blade_length
	return (thrust,torque)

print thrust_and_torque(0.001,100)
print thrust_and_torque(-0.001,100)


				





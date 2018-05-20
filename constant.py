constants = {
	'g':9.80655, #m/s^2
	#'k':1,  #thrust coefficient
	'Km':.0335, #motor constant
	'R':.07, #ohm
	'lambda_m': 0.00276,
	'N_P': 7,
	'k_theta':51.8, #stiffness
	'blade_c': 0,

	#body dimension
	'body_mass':2.8,  #kg
	'body_arm':.36,   #m
	#'body_thickness':.1, #m

	#rotor dimension
	'motor_mass':.075,   #kg
	'motor_radius':.029, #m
	'blade_mass':.021, #kg
	'blade_width':.02, #m 
	'blade_length':.23, #m 
	'pitch':.23, #m
	'rho':1.225, #kg/m^3
	'Cd':.01	
}

constants['S'] = 2.*constants['blade_length']*constants['blade_width']
constants['b'] = .5*constants['rho']*constants['Cd']*constants['S']*constants['blade_length']**3.
globals().update(constants)
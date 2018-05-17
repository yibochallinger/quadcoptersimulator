from sympy import *
from sympy.physics.mechanics import *
from sympy.functions import sign
from pydy.system import System
import numpy as np
from scipy.integrate import odeint
import math
from constant import *
import os

def square_inertia_xx_yy(m, l, z):
    return 1./12. * m*(l**2+z**2)

def square_inertia_zz(m, l):
    return 1./6. * m*l**2

def ball_inertia_xx_yy_zz(m, r):
    return 2./5. * m*r**2

def tube_inertia_xx_zz(m, r, l):
    return 1./4. * m*r**2 + 1./12. *m*l**2

def tube_inertia_yy(m, r):
    return 1./2. *m*r**2

def vector_in_frame(F, inlist):
    return F.x*inlist[0]+F.y*inlist[1]+F.z*inlist[2]

def quat_321_roll(quat):
    qr = quat[0]
    qi = quat[1]
    qj = quat[2]
    qk = quat[3]

    norm = sqrt(qr**2+qi**2+qj**2+qk**2)
    qr /= norm
    qi /= norm
    qj /= norm
    qk /= norm

    return math.atan2(2*(qr*qi+qj*qk), 1-2*(qi**2+qj**2))

def quat_321_pitch(quat):
    qr = quat[0]
    qi = quat[1]
    qj = quat[2]
    qk = quat[3]

    norm = sqrt(qr**2+qi**2+qj**2+qk**2)
    qr /= norm
    qi /= norm
    qj /= norm
    qk /= norm

    return math.asin(2*(qr*qj-qk*qi))
 
N = ReferenceFrame('N')
O = Point('O')
O.set_vel(N, 0)
gravity = g*N.z

#define states 
pos = dynamicsymbols('pos(0:3)')                               #position
pos_dot = dynamicsymbols('pos(0:3)',1)                         
vel = dynamicsymbols('vel(0:3)')                               #velocity
quat = dynamicsymbols('quat(0:4)')                             #quaternion
omega = dynamicsymbols('omega(0:3)')                           #angular velocity
motor_theta = dynamicsymbols('motor_theta(0:4)')               #motor angle
motor_theta_dot = dynamicsymbols('motor_theta(0:4)',1)
motor_omega = dynamicsymbols('motor_omega(0:4)')               #motor angular velocity
blade_theta = dynamicsymbols('blade_theta(0:8)')              #blade deflection angle
blade_theta_dot = dynamicsymbols('blade_theta(0:8)',1)
blade_omega = dynamicsymbols('blade_omega(0:8)')
motor_voltage = dynamicsymbols('motor_voltage(0:4)')


thrust_func = dynamicsymbols('thrust(0:8)')
torque_func = dynamicsymbols('torque(0:8)')

#define copter body frame

#body_inertia_xx_yy = square_inertia_xx_yy(body_mass, body_arm, body_thickness)
#body_inertia_zz = square_inertia_zz(body_mass, body_arm)

body_frame = N.orientnew('body_frame', 'Quaternion', quat)
body_frame.set_ang_vel(N, vector_in_frame(body_frame, omega))
body_inertia = inertia(body_frame, 0.042, 0.062 , 0.092)
body_masscenter = O.locatenew('body_masscenter', vector_in_frame(N, pos))
body_masscenter.set_vel(N, vector_in_frame(N, vel))
body_body = RigidBody('body', body_masscenter, body_frame, body_mass, (body_inertia, body_masscenter))

#copter eqns
bodies = []
forces = []
kdes = []

kdes.extend([pos_dot[i]-vel[i] for i in range(3)])
kdes.extend(kinematic_equations(omega, quat, 'quaternion'))
forces.append((body_masscenter, body_mass*gravity))
bodies.append(body_body)

#add motor frames

#motor_inertia_xx_yy_zz = ball_inertia_xx_yy_zz(motor_mass, motor_radius)

motor_frame =[]
motor_inertia = []
motor_masscenter = []
motor_body =[]
motorloc = [1/np.sqrt(2)*(body_frame.x+body_frame.y), -1/np.sqrt(2)*(body_frame.x+body_frame.y), 1/np.sqrt(2)*(body_frame.x-body_frame.y), 1/np.sqrt(2)*(-body_frame.x+body_frame.y)]
direction = [-1, -1, 1, 1]

for i in range (4):
    motor_frame.append(body_frame.orientnew('motor_frame[%u]' % (i,), 'Axis', [direction[i]*motor_theta[i], body_frame.z]))
    motor_inertia.append(inertia(motor_frame[i], 0.000031,  0.000031, 0.000063))
    motor_masscenter.append(body_masscenter.locatenew('motor_masscenter[%u]' % (i,), body_arm*motorloc[i]))
    motor_masscenter[i].v2pt_theory(body_masscenter, N, body_frame)
    motor_body.append(RigidBody('motor_body[%u]' % (i,), motor_masscenter[i], motor_frame[i], motor_mass, (motor_inertia[i], motor_masscenter[i])))
    forces.append((motor_frame[i], direction[i]*(motor_voltage[i]-Km*motor_omega[i])*Km/R*motor_frame[i].z))  #motor torque
    forces.append((body_frame, -direction[i]*(motor_voltage[i]-Km*motor_omega[i])*Km/R*motor_frame[i].z))  #body reaction
    forces.append((motor_masscenter[i], motor_mass*gravity))
    bodies.append(motor_body[i])
    kdes.append(motor_omega[i]-motor_theta_dot[i])
    #print motor_frame[i].ang_vel_in(N).to_matrix(body_frame)[2]

#add blade frames
#blade_inertia_xx_zz = tube_inertia_xx_zz(blade_mass, blade_width, blade_length)
#blade_inertia_yy = tube_inertia_yy(blade_mass, blade_width)

blade_frame =[]
blade_inertia = []
blade_masscenter = []
blade_body = []
force_center = []
blade_tip = []
blade_v_z = []
blade_v_2 = []
blade_thrust = []
blade_torque = []
for i in range(4):
    for j in range(2):
        blade_frame.append(motor_frame[i].orientnew('blade_frame[%u]' % (2*i+j,), 'Body', [pi*j, 0, direction[i]*blade_theta[2*i+j]], '321'))
        blade_inertia.append(inertia(blade_frame[2*i+j], 0.000091, 7e-7, 0.000091))
        blade_masscenter.append(motor_masscenter[i].locatenew('blade_masscenter[%u]' % (2*i+j,), -direction[i]*blade_length/2.0*blade_frame[2*i+j].y))
        blade_masscenter[2*i+j].v2pt_theory(motor_masscenter[i], N, blade_frame[2*i+j])
        #print (blade_masscenter[2*i+j].pos_from(motor_masscenter[i]).to_matrix(motor_frame[i]))
        blade_body.append(RigidBody('blade_body[%u]' % (2*i+j,), blade_masscenter[2*i+j], blade_frame[2*i+j], blade_mass, (blade_inertia[2*i+j], blade_masscenter[2*i+j])))
        force_center.append(motor_masscenter[i].locatenew('force_center[%u]' % (2*i+j,), -direction[i]*3.0/4.0*blade_length*blade_frame[2*i+j].y))
        force_center[2*i+j].v2pt_theory(blade_masscenter[-1], N, blade_frame[-1])
        blade_tip.append(motor_masscenter[i].locatenew('blade_tip[%u]' % (2*i+j,), blade_length*blade_frame[2*i+j].y))
        blade_tip[2*i+j].v2pt_theory(blade_masscenter[-1], N, blade_frame[-1])
        blade_v_z.append(force_center[i].vel(N).to_matrix(blade_frame[i])[2])
        blade_v_2.append(force_center[i].vel(N).to_matrix(blade_frame[i])[0])
        blade_thrust.append(Function('get_blade_thrust')(blade_v_z[-1], blade_v_2[-1]))
        blade_torque.append(Function('get_blade_torque')(blade_v_z[-1], blade_v_2[-1]))

        forces.append((force_center[2*i+j], -blade_thrust[-1]*blade_frame[2*i+j].z))  #thrust on blades
        forces.append((blade_frame[2*i+j], -direction[i]*k_theta*blade_theta[2*i+j]*blade_frame[2*i+j].x)) #blade stiffness torque
        forces.append((motor_frame[i], direction[i]*k_theta*blade_theta[2*i+j]*blade_frame[2*i+j].x)) #reaction torque of blade stiffness
        forces.append((blade_frame[2*i+j], direction[i]*blade_torque[-1]*blade_frame[2*i+j].z)) #torque on blades
        forces.append((blade_masscenter[2*i+j], blade_mass*gravity))
        bodies.append(blade_body[2*i+j])
        kdes.append(blade_omega[2*i+j]-blade_theta_dot[2*i+j])

        print blade_frame[-1].z.to_matrix(body_frame)[2]

q_sym = quat+pos+motor_theta+blade_theta
u_sym = omega+vel+motor_omega+blade_omega
in_sym = motor_voltage

print q_sym+u_sym

KM = KanesMethod(N, q_ind=q_sym, u_ind=u_sym, kd_eqs=kdes)

# Per https://github.com/sympy/sympy/issues/10945
try:
    (fr, frstar) = KM.kanes_equations(bodies=bodies, loads=forces)
except TypeError:
    (fr, frstar) = KM.kanes_equations(forces, bodies)

kdd = KM.kindiffdict()
mm = KM.mass_matrix_full.xreplace(kdd)
fo = KM.forcing_full.xreplace(kdd)

def count_expression(count_expr, within_expr):
    if hasattr(within_expr, "__getitem__"):
        return sum(map(lambda x: count_expression(count_expr, x), within_expr))
    else:
        return within_expr.count(count_expr)

def extractSubexpressions(inexprs, prefix='X', threshold=1, prev_subx=[]):
    # Perform common subexpression extraction on inexprs and existing subexpressions
    subexprs, outexprs = cse(inexprs+[x[1] for x in prev_subx], symbols=numbered_symbols('__TMP__'), order='none')

    subexprs = list(zip([x[0] for x in prev_subx], outexprs[len(inexprs):]))+subexprs
    outexprs = outexprs[:len(inexprs)]

    done = False
    while not done:
        done = True

        for i in reversed(range(len(subexprs))):
            from sympy.logic.boolalg import Boolean
            ops_saved = (count_expression(subexprs[i][0], [[x[1] for x in subexprs], outexprs])-1)*subexprs[i][1].count_ops()
            if ops_saved < threshold or isinstance(subexprs[i][1], Boolean):
                sub = dict([subexprs.pop(i)])
                subexprs = [(x[0],x[1].xreplace(sub)) for x in subexprs]
                outexprs = [x.xreplace(sub) for x in outexprs]
                done = False
                break

    for i in range(len(subexprs)):
        newSym = Symbol('%s_%u' % (prefix,i))
        sub = {subexprs[i][0]:newSym}
        subexprs[i] = (newSym,subexprs[i][1])
        subexprs = [(x[0],x[1].xreplace(sub)) for x in subexprs]
        outexprs = [x.xreplace(sub) for x in outexprs]

    outexprs = [x.as_mutable() if type(x) is ImmutableDenseMatrix else x for x in outexprs]

    return tuple(outexprs+[subexprs])

mm, fo, subx = extractSubexpressions([mm, fo], 'subx')

import os
if not os.path.exists('output'):
    os.makedirs('output')

with open('output/derivation.srepr', 'wb') as f:
    f.truncate()
    f.write(srepr({'mm':mm, 'fo':fo, 'subx':subx, 'inputs':in_sym, 'states':q_sym+u_sym, 'constants':constants, 'blade_v_z':Matrix(blade_v_z).xreplace(kdd), 'blade_v_2':Matrix(blade_v_2).xreplace(kdd)}))
    print('wrote output/derivation.srepr')


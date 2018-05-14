from highfidelitymodel import *
from pydy.viz.shapes import Cylinder, Sphere
import pydy.viz
from pydy.viz.visualization_frame import VisualizationFrame
from pydy.viz.scene import Scene

body_shape = Sphere(color='black', radius=0.1)
motor_shape = Sphere(color='black', radius=0.1)
arm_shape = Cylinder(radius=0.08, length=body_arm, color='blue')
blade_shape = Cylinder(radius=0.08, length=blade_length, color='red')
blade_tip_shape = Sphere(color='black', radius=.08)

body_viz_frame = VisualizationFrame(N, body_masscenter, body_shape)
motor_viz_frame=[]
arm_center=[]
arm_viz_frame=[]
blade_tip_viz_frame=[]
blade_center=[]
blade_viz_frame=[]

for i in range(4):
	motor_viz_frame.append(VisualizationFrame(N, motor_masscenter[i], motor_shape))
	arm_center.append(Point('a_c[%u]' % (i,)))
	arm_center[i].set_pos(body_masscenter, body_arm / 2. * motorloc[i])
	arm_viz_frame.append(VisualizationFrame('arm[%u]' % (i,), body_frame.orientnew('arm%u'%(i,), 'Axis', [i*pi/2,body_frame.z]), arm_center[i], arm_shape))
for i in range(4):
	for j in range(2):
		blade_tip_viz_frame.append(VisualizationFrame(N, blade_tip[2*i+j], blade_tip_shape))
		blade_center.append(Point('p_c[%u]' % (2*i+j,)))
		blade_center[2*i+j].set_pos(motor_masscenter[i], blade_length / 2. * blade_frame[2*i+j].y)
		blade_viz_frame.append(VisualizationFrame('blade[%u]' % (2*i+j,), blade_frame[2*i+j], blade_center[2*i+j], blade_shape))

scene = Scene(N, O)
scene.visualization_frames = [body_viz_frame] + arm_viz_frame + motor_viz_frame + blade_viz_frame
print scene.visualization_frames
scene.states_symbols = q_sym + u_sym
scene.states_trajectories = y
scene.constants = {}
scene.display()
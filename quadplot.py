import matplotlib.pyplot as plt
import csv
t=[]
v1=[]
v2=[]
v3=[]
v4=[]
mot_omega1=[]
mot_omega2=[]
mot_omega3=[]
mot_omega4=[]
mot_omega1_real=[]
mot_omega2_real=[]
mot_omega3_real=[]
mot_omega4_real=[]
w_x=[]
w_y=[]
w_z=[]
w_x_real=[]
w_y_real=[]
w_z_real=[]

with open('data.txt', 'r') as csvfile:
		plots = csv.reader(csvfile, delimiter=',')
		for row in plots:
			t.append(row[0])
			v1.append(row[1])
			v2.append(row[2])
			v3.append(row[3])
			v4.append(row[4])
			mot_omega1_real.append(row[11])
			mot_omega2_real.append(row[12])
			mot_omega3_real.append(row[13])
			mot_omega4_real.append(row[14])
			mot_omega1.append(row[15])
			mot_omega2.append(row[16])
			mot_omega3.append(row[17])
			mot_omega4.append(row[18])
			w_x.append(row[5])
			w_y.append(row[6])
			w_z.append(row[7])
			w_x_real.append(row[8])
			w_y_real.append(row[9])
			w_z_real.append(row[10])

plt.subplot(311)
plt.title('motor input voltages')
plt.ylabel('V')
plt.plot(t,v1, label='motor1 V')
plt.plot(t,v2, label='motor2 V')
plt.plot(t,v3, label='motor2 V')
plt.plot(t,v4, label='motor3 V')
plt.legend()

plt.subplot(312)
plt.title('motor speeds')
plt.ylabel('rad/s')
plt.plot(t,mot_omega1, label='motor1 omega')
plt.plot(t,mot_omega2, label='motor2 omega')
plt.plot(t,mot_omega3, label='motor2 omega')
plt.plot(t,mot_omega4, label='motor3 omega')
plt.plot(t,mot_omega1_real, label='motor1 omega')
plt.plot(t,mot_omega2_real, label='motor2 omega')
plt.plot(t,mot_omega3_real, label='motor2 omega')
plt.plot(t,mot_omega4_real, label='motor3 omega')
#plt.legend()

plt.subplot(313)
plt.title('')
plt.xlabel('seconds')
plt.ylabel('rad/s')
#plt.plot(t,w_x, label='ang_vel.x')
#plt.plot(t,w_x_real, label='x gyro')
#plt.plot(t,w_y, label='ang_vel.y')
#plt.plot(t,w_y_real, label='y gyro')
plt.plot(t,w_z, label='ang_vel.z')
plt.plot(t,w_z_real, label='z gyro')
plt.legend()
plt.suptitle('Flight comparison plots')
#plt.ylim(0,1)
plt.show()
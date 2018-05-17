#pragma once

#include <math.h>
#include <iostream>
using namespace std;

static inline void get_blade_thrust_torque(double v_z, double v2, double& ret_thrust, double& ret_torque) {

    //cout << v_z << " " << v2 << endl;
    double rad = 3.0/4.0 * blade_length;
    //local blade element setting angle
    double theta = atan(pitch/2.0/M_PI/rad); 
    double chord = blade_width;
    double phi = atan2(fabs(v_z),fabs(v2));
    //cout << phi << endl;
    double alpha =  v_z > 0 ? theta-phi: theta+phi;
    double cl = 6.2*alpha;
    double cd = .008-.003*cl+.01*pow(cl,2);
    //local velocity at blade
    double vlocal = sqrt(pow(v_z,2)+pow(v2,2));

    //thrust grading per bladeś
    double DtDr = 0.5*rho*pow(vlocal,2)*chord*(cl*cos(phi)-cd*sin(phi));
    //torque grading per blade
    double DqDr = 0.5*rho*pow(vlocal,2)*chord*rad*(cd*cos(phi)+cl*sin(phi));
    ret_thrust = v2 > 0 ? DtDr*blade_length: -DtDr*blade_length;//DtDr*blade_length;
    ret_torque = v2 > 0 ? DqDr*blade_length: -DqDr*blade_length;

    //cout << ret_thrust << " " << ret_torque << endl << endl;
}

static inline double get_blade_thrust(double v_z, double v2) {
    double thrust, torque;
    get_blade_thrust_torque(v_z, v2, thrust, torque);
    return thrust;
}

static inline double get_blade_torque(double v_z, double v2) {
    double thrust, torque;
    get_blade_thrust_torque(v_z, v2, thrust, torque);
    return torque;
}

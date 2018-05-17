#pragma once

#include <math.h>
#include <iostream>
using namespace std;

static inline void get_blade_thrust_torque(double v_z, double v2, double& ret_thrust, double& ret_torque) {

    //cout << v_z << " " << v2 << endl;
    double rad = 3.0/4.0 * blade_length;
    //local blade element setting angle
    double theta = atan(pitch/2.0/M_PI/rad); 
    //initial guess of inflow and swirl factor
    double a = 0.1;
    int total = 1;
    double chord = blade_width;
    double DtDr;
    double DqDr;
    while(true){
        //axial velocity
        double v0 = v_z*(1+a);
        double phi = atan2(fabs(v0),fabs(v2));
//         cout << phi << endl;
        double alpha =  v0 > 0 ? theta-phi: theta+phi;
        double cl = 6.2*alpha;
        double cd = .008-.003*cl+.01*pow(cl,2);
        //local velocity at blade
        double vlocal = sqrt(pow(v0,2)+pow(v2,2));

        //thrust grading per bladeś
        DtDr = 0.5*rho*pow(vlocal,2)*chord*(cl*cos(phi)-cd*sin(phi));
        //torque grading per blade
        DqDr = 0.5*rho*pow(vlocal,2)*chord*rad*(cd*cos(phi)+cl*sin(phi));
        //momentum check on inflow and swirl factors
        double tem1 = DtDr/(4.0*M_PI*rad*rho*pow(v_z,2)*(1+a));
        double anew = 0.5*(a+tem1);
        
        total = total+1;
        //check to see if iteration stuck or convergence
        if (total>500 || abs(anew-a)<1.0e-5) {
            break;
        }
        a = anew;
    }

    ret_thrust = v2 > 0 ? DtDr*blade_length: -DtDr*blade_length;//DtDr*blade_length;
    ret_torque = v2 > 0 ? DqDr*blade_length: -DqDr*blade_length;

    //cout << ret_thrust << " " << ret_torque << endl << endl;
}

static inline double get_blade_thrust(double v_z, double v2) {
    double thrust, torque;
    get_blade_thrust_torque(-v_z, v2, thrust, torque);
    return thrust;
}

static inline double get_blade_torque(double v_z, double v2) {
    double thrust, torque;
    get_blade_thrust_torque(-v_z, v2, thrust, torque);
    return torque;
}

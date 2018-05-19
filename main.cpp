#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <math.h>
#include "output/quadcoptereqns.h"
#include <specifieds.h.backup>

#define SQ(x) ((x)*(x))

using namespace Eigen;
using namespace std;

static Matrix<double,37,1> states;
static Matrix<double,4,1> inputs;


Matrix<double,37,1> compute_xdot(Matrix<double,37,1> states, Matrix<double,4,1> inputs) {
    Matrix<double,37,37> mm = get_mm(states, inputs);
    Matrix<double,37,1> fo = get_fo(states, inputs);
    return mm.llt().solve(fo);;
}

void wrap_2pi(double& x) {
    if (x >= 2*M_PI) {
        x -= 2*M_PI;
    } else if (x < 0) {
        x += 2*M_PI;
    }
}

int main() {
    states << 0, 0, 0, 1, // quat 0 1 2 3
              0, 0, 0, // pos 4 5 6
              0, 0, 0, 0, // motor angle 7 8 9 10
              0, 0, 0, 0, 0, 0, 0, 0, //blade angle 11 12 13 14 15 16 17 18
              0, 0, 0, // ang vel 19 20 21
              0, 0, 0.01, // vel 22 23 24
              380, 380, 380, 380, //motor ang vel 25 26 27 28
              0, 0, 0, 0, 0, 0, 0, 0; // blade ang vel 29 30 31 32 33 34 35 36


    inputs << 9.4, 
              9.4, 
              9.4, 
              9.4; // motor voltage

    double t = 0;
    ofstream datafile("data.txt");
    get_blade_thrust(0.1, 50);
    for (int a=0; a<10; a++) {
        Matrix<double,37,1> xdot;
        for (int i=0; i<1000; i++) {
            double dt = 0.001;
            t += dt;

            xdot = compute_xdot(states,inputs);

            Matrix<double,37,1> k1 = dt * compute_xdot(states,inputs);
            Matrix<double,37,1> k2 = dt * compute_xdot(states+k1/2,inputs);
            Matrix<double,37,1> k3 = dt * compute_xdot(states+k2/2,inputs);
            Matrix<double,37,1> k4 = dt * compute_xdot(states+k3,inputs);

            Matrix<double,37,1> xdelta = (k1+2*k2+2*k3+k4)/6;

            if (!xdelta.allFinite()) {
                cout << "FAIL2" << endl << "state: " << endl << states << endl << "input: " << endl << inputs << endl << endl;
                return 0;
            }

            states += xdelta;

            double norm = sqrt(SQ(states(0,0))+SQ(states(1,0))+SQ(states(2,0))+SQ(states(3,0)));
            states(0,0) /= norm;
            states(1,0) /= norm;
            states(2,0) /= norm;
            states(3,0) /= norm;

            wrap_2pi(states(7,0));
            wrap_2pi(states(8,0));
            wrap_2pi(states(9,0));
            wrap_2pi(states(10,0));
        }


        cout << "xdot:" << endl << xdot << endl << endl;


        Matrix<double, 8,1> blade_v_z = get_blade_v_z(states, inputs);
        Matrix<double, 8,1> blade_v_2 = get_blade_v_2(states, inputs);

        cout << "blade_v_2" << endl << blade_v_2 << endl << endl;
        cout << "blade_v_z" << endl << blade_v_z << endl << endl;
        Matrix<double, 8,1> blade_thrust;
        Matrix<double, 8,1> blade_torque;
        for (int j=0; j<8; j++) {
            blade_thrust(j,0) = get_blade_thrust(blade_v_z(j,0), blade_v_2(j,0));
            blade_torque(j,0) = get_blade_torque(blade_v_z(j,0), blade_v_2(j,0));
        }

         cout << "blade_thrust" << endl << blade_thrust  << endl << endl << "blade_torque" << endl << blade_torque << endl << endl;

        double qr = states(0,0);
        double qi = states(1,0);
        double qj = states(2,0);
        double qk = states(3,0);

        //double norm = sqrt(SQ(states(0,0))+SQ(states(1,0))+SQ(states(2,0))+SQ(states(3,0)));

        double pitch = 180*asin(2*(qr*qj-qk*qi))/M_PI;
    
        cout << "       pitch: " << pitch << endl;

        cout << "       quat0: " << states(0,0) << endl
             << "       quat1: " << states(1,0) << endl
             << "       quat2: " << states(2,0) << endl
             << "       quat3: " << states(3,0) << endl
             << "        pos0: " << states(4,0) << endl
             << "        pos1: " << states(5,0) << endl
             << "        pos2: " << states(6,0) << endl
             << "motor_angle0: " << states(7,0) << endl
             << "motor_angle1: " << states(8,0) << endl
             << "motor_angle2: " << states(9,0) << endl
             << "motor_angle3: " << states(10,0) << endl
             << "      omega0: " << states(19,0) << endl
             << "      omega1: " << states(20,0) << endl
             << "      omega2: " << states(21,0) << endl
             << "        vel0: " << states(22,0) << endl
             << "        vel1: " << states(23,0) << endl
             << "        vel2: " << states(24,0) << endl
             << "motor_omega0: " << states(25,0) << endl
             << "motor_omega1: " << states(26,0) << endl
             << "motor_omega2: " << states(27,0) << endl
             << "motor_omega3: " << states(28,0) << endl
             << "blade_theta0: " << states(11,0) << endl
             << "blade_theta1: " << states(12,0) << endl
             << "blade_theta2: " << states(13,0) << endl
             << "blade_theta3: " << states(14,0) << endl
             << "blade_theta4: " << states(15,0) << endl
             << "blade_theta5: " << states(16,0) << endl
             << "blade_theta6: " << states(17,0) << endl
             << "blade_theta7: " << states(18,0) << endl
             ;

        datafile << "       quat0: " << states(0,0) << endl
                 << "       quat1: " << states(1,0) << endl
                 << "       quat2: " << states(2,0) << endl
                 << "       quat3: " << states(3,0) << endl
                 << "        pos0: " << states(4,0) << endl
                 << "        pos1: " << states(5,0) << endl
                 << "        pos2: " << states(6,0) << endl
                 << "motor_angle0: " << states(7,0) << endl
                 << "motor_angle1: " << states(8,0) << endl
                 << "motor_angle2: " << states(9,0) << endl
                 << "motor_angle3: " << states(10,0) << endl
                 << "      omega0: " << states(19,0) << endl
                 << "      omega1: " << states(20,0) << endl
                 << "      omega2: " << states(21,0) << endl
                 << "        vel0: " << states(22,0) << endl
                 << "        vel1: " << states(23,0) << endl
                 << "        vel2: " << states(24,0) << endl
                 << "motor_omega0: " << states(25,0) << endl
                 << "motor_omega1: " << states(26,0) << endl
                 << "motor_omega2: " << states(27,0) << endl
                 << "motor_omega3: " << states(28,0) << endl
                 << "blade_theta0: " << states(11,0) << endl
                 << "blade_theta1: " << states(12,0) << endl
                 << "blade_theta2: " << states(13,0) << endl
                 << "blade_theta3: " << states(14,0) << endl
                 << "blade_theta4: " << states(15,0) << endl
                 << "blade_theta5: " << states(16,0) << endl
                 << "blade_theta6: " << states(17,0) << endl
                 << "blade_theta7: " << states(18,0) << endl
                 ;
            cout << "           t: " << t << endl;
    }

//         cout << "fo:" << endl << fo << endl << endl << "mm:" << endl << mm << endl;
    datafile.close();
    return 0;
}

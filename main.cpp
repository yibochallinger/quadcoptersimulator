#include <Eigen/Dense>
#include <iostream>
#include <math.h>
#include "output/quadcoptereqns.h"

#define SQ(x) ((x)*(x))

using namespace Eigen;
using namespace std;

static Matrix<double,37,1> states;
static Matrix<double,4,1> inputs;



int main() {
    states << 0, 0, 0, 1, // quat 0 1 2 3
              0, 0, 0, // pos 4 5 6
              0, 0, 0, 0, // motor angle 7 8 9 10
              -.01, -.01, -.01, -.01, -.01, -.01, -.01, -.01, //blade angle 11 12 13 14 15 16 17 18
              0, 0, 0, // ang vel 19 20 21
              0, 0, 0.01, // vel 22 23 24
              0, 0, 0, 0, //motor ang vel 25 26 27 28
              0, 0, 0, 0, 0, 0, 0, 0; // blade ang vel 29 30 31 32 33 34 35 36

    inputs << 0, 
              0, 
              0, 
              0; // motor voltage

    double t = 0;
    for (int a=0; a<1; a++) {
        for (int i=0; i<1; i++) {
            double dt = 0.0000001;
            t += dt;

            Matrix<double,37,37> mm = get_mm(states, inputs);
            Matrix<double,37,1> fo = get_fo(states, inputs);

            if (!mm.allFinite() || !fo.allFinite()) {
                cout << "FAIL1" << endl << "state: " << endl << states << endl << "input: " << endl << inputs << endl << endl;
                return 0;
            }

            Matrix<double,37,1> xdot = mm.fullPivLu().solve(fo);
            cout<<xdot*dt<<endl<<endl;

            states += xdot * dt;

            if (!xdot.allFinite()) {
                cout << "FAIL2" << endl << "state: " << endl << states << endl << "input: " << endl << inputs << endl << endl;
                return 0;
            }

            double norm = sqrt(SQ(states(0,0))+SQ(states(1,0))+SQ(states(2,0))+SQ(states(3,0)));
            states(0,0) /= norm;
            states(1,0) /= norm;
            states(2,0) /= norm;
            states(3,0) /= norm;
        }

        double qr = states(0,0);
        double qi = states(1,0);
        double qj = states(2,0);
        double qk = states(3,0);

        //double norm = sqrt(SQ(states(0,0))+SQ(states(1,0))+SQ(states(2,0))+SQ(states(3,0)));

        double pitch = 180*asin(2*(qr*qj-qk*qi))/M_PI;

        cout << "           t: " << t << endl;

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
             
    }

//         cout << "fo:" << endl << fo << endl << endl << "mm:" << endl << mm << endl;


    return 0;
}
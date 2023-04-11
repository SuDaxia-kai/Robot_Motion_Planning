#include "trajectory_generator_waypoint.h"
#include <stdio.h>
#include <ros/ros.h>
#include <ros/console.h>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;    
using namespace Eigen;

TrajectoryGeneratorWaypoint::TrajectoryGeneratorWaypoint(){}
TrajectoryGeneratorWaypoint::~TrajectoryGeneratorWaypoint(){}

//define factorial function, input i, output i!
int TrajectoryGeneratorWaypoint::Factorial(int x)
{
    int fac = 1;
    for(int i = x; i > 0; i--)
        fac = fac * i;
    return fac;
}


Eigen::MatrixXd Generate_Tm(double t, int n_order, int r)
{
    MatrixXd Tm = MatrixXd::Zero(1, n_order+1);
    int F = 1;
    for(int i = r; i <= n_order; i++)
    {
        for(int j = i; j > i - r; j--)
            F *= j;
        Tm(i) = F * pow(t, i-r); 
    }
    return Tm;
}

/*

    STEP 2: Learn the "Closed-form solution to minimum snap" in L5, then finish this PolyQPGeneration function

    variable declaration: input       const int d_order,                    // the order of derivative
                                      const Eigen::MatrixXd &Path,          // waypoints coordinates (3d)
                                      const Eigen::MatrixXd &Vel,           // boundary velocity
                                      const Eigen::MatrixXd &Acc,           // boundary acceleration
                                      const Eigen::VectorXd &Time)          // time allocation in each segment
                          output      MatrixXd PolyCoeff(m, 3 * p_num1d);   // position(x,y,z), so we need (3 * p_num1d) coefficients

*/

Eigen::MatrixXd TrajectoryGeneratorWaypoint::PolyQPGeneration(
            const int d_order,                    // the order of derivative
            const Eigen::MatrixXd &Path,          // waypoints coordinates (3d)
            const Eigen::MatrixXd &Vel,           // boundary velocity
            const Eigen::MatrixXd &Acc,           // boundary acceleration
            const Eigen::VectorXd &Time)          // time allocation in each segment
{
    // enforce initial and final velocity and accleration, for higher order derivatives, just assume them be 0;
    int p_order   = 2 * d_order - 1;              // the order of polynomial
    int p_num1d   = p_order + 1;                  // the number of variables in each segment

    int m = Time.size();                          // the number of segments
    MatrixXd PolyCoeff = MatrixXd::Zero(m, 3 * p_num1d);           // position(x,y,z), so we need (3 * p_num1d) coefficients
    VectorXd Px(p_num1d * m), Py(p_num1d * m), Pz(p_num1d * m);

    /*
        Method1: Using OOQP to solve this problem 
        Matrix A and D are the constraint matrixs
    */
    MatrixXd A = MatrixXd::Zero(2 * d_order - 5 + 5 * m, m*p_num1d);
    MatrixXd D = MatrixXd::Zero(2 * d_order - 5 + 5 * m, 3); 
    /** 
        Method2: Close-form solution for this problem
        whose form is [ D_F, D_P ]'
        D_F is 2*d_order + (k - 1)
        D_P is 3 * (k-1)
    */
    MatrixXd M = MatrixXd::Zero(2 * d_order * m, m * p_num1d);
    MatrixXd C_T = MatrixXd::Zero(2 * d_order + (m - 1) * 2 * d_order, 2 * d_order + m - 1 + 3 * (m - 1));  // C is the selection matrix
    MatrixXd DF = MatrixXd::Zero(2 * d_order + m - 1, 3);
    MatrixXd DP = MatrixXd::Zero(3 * (m - 1), 3);                     // It is the solution of close-form
    MatrixXd DD = MatrixXd::Zero(2 * d_order + m - 1 + 3 * (m - 1), 3);  

    /*
        public parameter:
        matrix Q
    */
    MatrixXd Q = MatrixXd::Zero(m * p_num1d, m * p_num1d);

    /*   Produce Mapping Matrix M to the entire trajectory, M is a mapping matrix that maps polynomial coefficients to derivatives.   */
    for( auto seg = 0; seg < m; seg++ )
    {
        M.block( seg * 2 * d_order, seg * p_num1d, 2 * d_order, p_num1d ) << Generate_Tm(0, p_order, 0),
                                                                             Generate_Tm(0, p_order, 1),
                                                                             Generate_Tm(0, p_order, 2),
                                                                             Generate_Tm(0, p_order, 3),
                                                                             Generate_Tm(Time[seg], p_order, 0),
                                                                             Generate_Tm(Time[seg], p_order, 1),
                                                                             Generate_Tm(Time[seg], p_order, 2),
                                                                             Generate_Tm(Time[seg], p_order, 3);
    }
    cout << "M IS " << endl << M << endl;
    /*  Produce the fix constrain matrix DF */
    DF.block(0, 0, 2*d_order, 3) << Path.row(0), Vel.row(0), Acc.row(0), MatrixXd::Zero(1, 3),
                                    Path.row(Path.rows()-1), Vel.row(Vel.rows()-1), Acc.row(Acc.rows()-1), MatrixXd::Zero(1, 3);
    for( auto seg_point = 0; seg_point < m - 1; seg_point++ )
    {
        DF.row(seg_point + 2*d_order) = Path.row(seg_point + 1);
    }

    /*  Produce selection matrix C  */
    C_T.block(0, 0, d_order, d_order) = MatrixXd::Identity(d_order, d_order);
    C_T.block((2*m-1) * d_order, d_order, d_order, d_order) = MatrixXd::Identity(d_order, d_order);
    for( auto num_cp = 0; num_cp < (m - 1); num_cp++ )
    {
        C_T(2*num_cp*d_order + d_order, 2 * d_order + num_cp) = 1;                              // p
        C_T(2*num_cp*d_order + d_order + 1, 2 * d_order + (m - 1) + num_cp * 3) = 1;            // v
        C_T(2*num_cp*d_order + d_order + 2, 2 * d_order + (m - 1) + num_cp * 3 + 1) = 1;        // a
        C_T(2*num_cp*d_order + d_order + 3, 2 * d_order + (m - 1) + num_cp * 3 + 2) = 1;        // j

        C_T((num_cp + 1) * 2*d_order, 2 * d_order + num_cp) = 1;                                // p
        C_T((num_cp + 1) * 2*d_order + 1, 2 * d_order + (m - 1) + num_cp * 3) = 1;              // v
        C_T((num_cp + 1) * 2*d_order + 2, 2 * d_order + (m - 1) + num_cp * 3 + 1) = 1;          // a
        C_T((num_cp + 1) * 2*d_order + 3, 2 * d_order + (m - 1) + num_cp * 3 + 2) = 1;          // j
    }


    /*  Produce matrix Q */
    for( auto num_Q = 0; num_Q < m; num_Q++ )
    {
        MatrixXd QJ = MatrixXd::Zero(p_num1d, p_num1d);
        for( auto row = (d_order - 1); row < p_num1d; row++ )
        {
            for( auto line = row; line < p_num1d; line++ )
            {
                // refer to the blog
                QJ(row, line) = QJ(line, row) = row*(row-1)*(row-2)*(row-3)*line*(line-1)*(line-2)*(line-3)*pow(Time(num_Q), row+line-7)/row+line-7;
            }
        }
        Q.block(num_Q*p_num1d, num_Q*p_num1d, p_num1d, p_num1d) = QJ;
    }
    // cout << " Q matrix is " << endl << Q << endl;

    /*  Produce the RFP, Rpp and R  */
    MatrixXd R = C_T.transpose() * M.inverse().transpose() * Q * M.inverse() * C_T;
    MatrixXd RFP = R.block(0, 2 * d_order + m - 1, 2 * d_order + m - 1, 3 * (m - 1));
    MatrixXd RPP = R.block(2 * d_order + m - 1, 2 * d_order + m - 1, 3 * (m - 1), 3 * (m - 1));

    DP = -RPP.inverse() * RFP.transpose() * DF;
    DD << DF,
          DP;
    cout << "DD is " << endl << DD << endl;
    MatrixXd P = M.inverse() * C_T * DD;
    MatrixXd temp_p = MatrixXd::Zero(1, 3*p_num1d);
    for( auto seg = 0; seg < m; seg++ )
    {
        temp_p << P.block(seg * p_num1d, 0, p_num1d, 1).transpose(),
                  P.block(seg * p_num1d, 1, p_num1d, 1).transpose(),
                  P.block(seg * p_num1d, 2, p_num1d, 1).transpose();
        PolyCoeff.row(seg) = temp_p;
    }
    /*  Method1: Using OOQP to solve this problem  */
    // A.block(0, 0, d_order, p_num1d) << Generate_Tm(0, p_order, 0),
    //                                    Generate_Tm(0, p_order, 1),
    //                                    Generate_Tm(0, p_order, 2),
    //                                    Generate_Tm(0, p_order, 3);
    
                                         
    // A.block(d_order, A.cols() - p_num1d, d_order, p_num1d) << Generate_Tm(Time(-1), p_order, 0),
    //                                                           Generate_Tm(Time(-1), p_order, 1),
    //                                                           Generate_Tm(Time(-1), p_order, 2),
    //                                                           Generate_Tm(Time(-1), p_order, 3);

    // D.block(0, 0, 2*d_order, 3) << Path(0),Vel(0),Acc(0), MatrixXd::Zero(1, 3),
    //                                Path(-1),Vel(-1),Acc(-1), MatrixXd::Zero(1, 3);

    
    // // Waypoint Is Firm
    // for(int i = 0; i < m - 1; i++)
    // {
    //     A.block(2*d_order + i, i * p_num1d, 0, p_num1d) = Generate_Tm(Time(i), p_order, 0);
    //     D.row(i + 2*d_order) = Path.row(i + 1);
    // }

    // // Position continuity constrain between each 2 segments
    // for(int i = 0; i < m - 1; i++)
    // {
    //     A.block((2*d_order + m - 1), i * p_num1d, 0, 2 * p_num1d) << Generate_Tm(Time(i), p_order, 0),
    //                                                                  -Generate_Tm(0, p_order, 0);
    //     D.row(i + 2*d_order + m - 1) = MatrixXd::Zero(1, 3);
    // }

    // // Velocity continuity constrain between each 2 segments
    // for(int i = 0; i < m - 1; i++)
    // {
    //     A.block((2*d_order + 2*(m - 1)), i * p_num1d, 0, 2 * p_num1d) << Generate_Tm(Time(i), p_order, 1),
    //                                                                      -Generate_Tm(0, p_order, 1);
    //     D.row(i + 2*d_order + 2*(m - 1)) = MatrixXd::Zero(1, 3);
    // }

    // // acceleration continuity constrain between each 2 segments
    // for(int i = 0; i < m - 1; i++)
    // {
    //     A.block((2*d_order + 3*(m - 1)), i * p_num1d, 0, 2 * p_num1d) << Generate_Tm(Time(i), p_order, 2),
    //                                                                      -Generate_Tm(0, p_order, 2);
    //     D.row(i + 2*d_order + 3*(m - 1)) = MatrixXd::Zero(1, 3);
    // }

    // // jerk continuity constrain between each 2 segments
    // for(int i = 0; i < m - 1; i++)
    // {
    //     A.block((2*d_order + 4*(m - 1)), i * p_num1d, 0, 2 * p_num1d) << Generate_Tm(Time(i), p_order, 3),
    //                                                                      -Generate_Tm(0, p_order, 3);
    //     D.row(i + 2*d_order + 4*(m - 1)) = MatrixXd::Zero(1, 3);
    // }
    
    

    /*   Produce the dereivatives in X, Y and Z axis directly.  */
    












    /*   Produce the Minimum Snap cost function, the Hessian Matrix   */




    return PolyCoeff;
}

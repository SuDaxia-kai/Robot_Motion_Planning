#include <hw_tool.h>

using namespace std;
using namespace Eigen;

void Homeworktool::initGridMap(double _resolution, Vector3d global_xyz_l, Vector3d global_xyz_u, int max_x_id, int max_y_id, int max_z_id)
{   
    gl_xl = global_xyz_l(0);
    gl_yl = global_xyz_l(1);
    gl_zl = global_xyz_l(2);

    gl_xu = global_xyz_u(0);
    gl_yu = global_xyz_u(1);
    gl_zu = global_xyz_u(2);
    
    GLX_SIZE = max_x_id;
    GLY_SIZE = max_y_id;
    GLZ_SIZE = max_z_id;
    GLYZ_SIZE  = GLY_SIZE * GLZ_SIZE;
    GLXYZ_SIZE = GLX_SIZE * GLYZ_SIZE;

    resolution = _resolution;
    inv_resolution = 1.0 / _resolution;    

    data = new uint8_t[GLXYZ_SIZE];
    memset(data, 0, GLXYZ_SIZE * sizeof(uint8_t));
}

void Homeworktool::setObs(const double coord_x, const double coord_y, const double coord_z)
{   
    if( coord_x < gl_xl  || coord_y < gl_yl  || coord_z <  gl_zl || 
        coord_x >= gl_xu || coord_y >= gl_yu || coord_z >= gl_zu )
        return;

    int idx_x = static_cast<int>( (coord_x - gl_xl) * inv_resolution);
    int idx_y = static_cast<int>( (coord_y - gl_yl) * inv_resolution);
    int idx_z = static_cast<int>( (coord_z - gl_zl) * inv_resolution);      
    
    data[idx_x * GLYZ_SIZE + idx_y * GLZ_SIZE + idx_z] = 1;
}

bool Homeworktool::isObsFree(const double coord_x, const double coord_y, const double coord_z)
{
    Vector3d pt;
    Vector3i idx;
    
    pt(0) = coord_x;
    pt(1) = coord_y;
    pt(2) = coord_z;
    idx = coord2gridIndex(pt);

    int idx_x = idx(0);
    int idx_y = idx(1);
    int idx_z = idx(2);

    return (idx_x >= 0 && idx_x < GLX_SIZE && idx_y >= 0 && idx_y < GLY_SIZE && idx_z >= 0 && idx_z < GLZ_SIZE && 
           (data[idx_x * GLYZ_SIZE + idx_y * GLZ_SIZE + idx_z] < 1));
}

Vector3d Homeworktool::gridIndex2coord(const Vector3i & index) 
{
    Vector3d pt;

    pt(0) = ((double)index(0) + 0.5) * resolution + gl_xl;
    pt(1) = ((double)index(1) + 0.5) * resolution + gl_yl;
    pt(2) = ((double)index(2) + 0.5) * resolution + gl_zl;

    return pt;
}

Vector3i Homeworktool::coord2gridIndex(const Vector3d & pt) 
{
    Vector3i idx;
    idx <<  min( max( int( (pt(0) - gl_xl) * inv_resolution), 0), GLX_SIZE - 1),
            min( max( int( (pt(1) - gl_yl) * inv_resolution), 0), GLY_SIZE - 1),
            min( max( int( (pt(2) - gl_zl) * inv_resolution), 0), GLZ_SIZE - 1);                  
  
    return idx;
}

Eigen::Vector3d Homeworktool::coordRounding(const Eigen::Vector3d & coord)
{
    return gridIndex2coord(coord2gridIndex(coord));
}

double Homeworktool::OptimalBVP(Eigen::Vector3d _start_position,Eigen::Vector3d _start_velocity,Eigen::Vector3d _target_position)
{
    double optimal_cost = 100000; // this just to initial the optimal_cost, you can delete it 
    /*
    

    STEP 2: go to the hw_tool.cpp and finish the function Homeworktool::OptimalBVP
    the solving process has been given in the document

    because the final point of trajectory is the start point of OBVP, so we input the pos,vel to the OBVP

    after finish Homeworktool::OptimalBVP, the Trajctory_Cost will record the optimal cost of this trajectory


    */

    /**
     *  利用Matlab推导出来的多项式
     *  (- 3*v0x^2 - 3*v0y^2 - 3*v0z^2)*T^2 + 
     *  (12*pfx*v0x - 12*p0y*v0y - 12*p0z*v0z - 12*p0x*v0x + 12*pfy*v0y + 12*pfz*v0z)*T 
     *  -9*p0x^2 + 18*p0x*pfx - 9*p0y^2 + 18*p0y*pfy - 9*p0z^2 + 18*p0z*pfz - 9*pfx^2 - 9*pfy^2 - 9*pfz^2) 
    */

    Vector3d c;
    c << 9*(pow(_start_position[0], 2) + pow(_start_position[1], 2) + pow(_start_position[2], 2) + pow(_target_position[0], 0) + pow(_target_position[1], 1) + pow(_target_position[2], 2) + 18.0*double((_start_position.transpose() * _target_position))),
         -12*(double((_target_position.transpose() * _start_velocity)) - double((_start_position.transpose() * _start_velocity))), 
         -3*(_start_velocity.squaredNorm());
    Matrix4d M = Matrix4d::Zero(4, 4);
    M(1, 0) = M(2, 1) = M(3, 2) = 1.0;
    M(0, 3) = c[0];
    M(1, 3) = c[1];
    M(2, 3) = c[2];

    Matrix<complex<double>, Dynamic, Dynamic> M_eig = M.eigenvalues();
    int iter = 0, real_flag = false;
    for(; iter < M.rows(); iter++)
    {
        if(fabs(M_eig(iter, 0).imag()-0.0) < 1e-8 && M_eig(iter, 0).real() > 0)
        {
            real_flag = true;
            break;
        }
    }
    if(real_flag == true)
    {
        double T = M_eig(iter, 0).real();
        Vector3d Dp = _target_position - _start_position - _start_velocity * T;
        optimal_cost = T + 3*Dp.squaredNorm()/pow(T, 3.0);
    }
    return optimal_cost;
}

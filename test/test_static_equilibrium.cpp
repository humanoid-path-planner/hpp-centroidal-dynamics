#include <vector>
#include <Eigen/Dense>
#include <Eigen/src/Core/util/Macros.h>

#include "cdd/cddmp.h"
#include "cdd/setoper.h"
#include "cdd/cddtypes.h"
#include "cdd/cdd.h"

#include <robust-equilibrium-lib/static_equilibrium.hh>

using namespace robust_equilibrium;
using namespace Eigen;

typedef Eigen::AngleAxis<value_type> angle_axis_t;


void init_library()
{
  dd_set_global_constants();dd_debug = false;
}

void release_library()
{
  //dd_free_global_constants();
}

bool uniform(Cref_matrixXX lower_bounds, Cref_matrixXX upper_bounds, Ref_matrixXX out)
{

  assert(lower_bounds.rows()==out.rows());
  assert(upper_bounds.rows()==out.rows());
  assert(lower_bounds.cols()==out.cols());
  assert(upper_bounds.cols()==out.cols());
  for(int i=0; i<out.rows(); i++)
    for(int j=0; j<out.cols(); j++)
      out(i,j) = (rand()/ value_type(RAND_MAX))*(upper_bounds(i,j)-lower_bounds(i,j)) + lower_bounds(i,j);
}

bool euler_matrix(double roll, double pitch, double yaw, Ref_rotation R)
{
  R = (angle_axis_t(roll, Vector3::UnitX())
       * angle_axis_t(pitch, Vector3::UnitY())
       * angle_axis_t(yaw, Vector3::UnitZ())).toRotationMatrix();
}

bool generate_rectangle_contacts(double lx, double ly, Cref_vector3 pos, Cref_vector3 rpy,
                                 Ref_matrix43 p, Ref_matrix43 N)
{
  // compute rotation matrix
  Quaternion<value_type> q = angle_axis_t(rpy(0), Vector3::UnitX())
                             * angle_axis_t(rpy(1), Vector3::UnitY())
                             * angle_axis_t(rpy(2), Vector3::UnitZ());
  //euler_matrix(rpy[0], rpy[1], rpy[2], R);
  // contact points in local frame
  p << lx,  ly, 0,
       lx, -ly, 0,
      -lx, -ly, 0,
      -lx,  ly, 0;
  // contact points in world frame
  p.row(0) = pos + q*p.row(0);
  p.row(1) = pos + q*p.row(1);
  p.row(2) = pos + q*p.row(2);
  p.row(3) = pos + q*p.row(3);
  // normal direction in local frame
  RVector3 n;
  n << 0, 0, 1;
  // normal directions in world frame
  n = q*n;
  N << n, n, n, n;
  return true;
}

int main()
{
  int ret = 0;
  double mu = 0.3;  // friction coefficient
  unsigned int generatorsPerContact = 4;
  unsigned int N_CONTACTS = 8;
  unsigned int N_POINTS = 50;   // number of com positions to test
  double MIN_FEET_DISTANCE = 0.3;
  double LX = 0.5*0.2172;        // half foot size in x direction
  double LY = 0.5*0.138;         // half foot size in y direction
  RVector3 CONTACT_POINT_LOWER_BOUNDS, CONTACT_POINT_UPPER_BOUNDS;
  CONTACT_POINT_LOWER_BOUNDS <<-0.5, -0.5,  0.0;
  CONTACT_POINT_UPPER_BOUNDS << 0.5,  0.5,  1.0;
  double gamma = atan(mu);   // half friction cone angle
  RVector3 RPY_LOWER_BOUNDS, RPY_UPPER_BOUNDS;
  RPY_LOWER_BOUNDS << -2*gamma, -2*gamma, -M_PI;
  RPY_UPPER_BOUNDS << +2*gamma, +2*gamma, +M_PI;
  double X_MARG = 0.07;
  double Y_MARG = 0.07;

  StaticEquilibrium solver(generatorsPerContact, SOLVER_LP_QPOASES);
  MatrixXX contact_pos = MatrixXX::Zero(N_CONTACTS, 3);
  MatrixXX contact_rpy = MatrixXX::Zero(N_CONTACTS, 3);
  MatrixXX p = MatrixXX::Zero(4*N_CONTACTS,3); // contact points
  MatrixXX N = MatrixXX::Zero(4*N_CONTACTS,3); // contact normals
  VectorX frictionCoefficients(N_CONTACTS);
  frictionCoefficients.fill(mu);

  // Generate contact positions and orientations
  bool collision;
  for(int i=0; i<N_CONTACTS; i++)
  {
    while(true) // generate contact position
    {
      uniform(CONTACT_POINT_LOWER_BOUNDS, CONTACT_POINT_UPPER_BOUNDS, contact_pos.row(i));
      collision = false;
      for(int j=0; j<i-1; j++)
        if((contact_pos.row(i)-contact_pos.row(j)).norm() < MIN_FEET_DISTANCE)
          collision = true;
      if(collision==false)
        break;
    }

    // generate contact orientation
    uniform(RPY_LOWER_BOUNDS, RPY_UPPER_BOUNDS, contact_rpy.row(i));
    generate_rectangle_contacts(LX, LY, contact_pos.row(i), contact_rpy.row(i),
                                p.middleRows<4>(i*4), N.middleRows<4>(i*4));
  }

  solver.setNewContacts(p, N, frictionCoefficients, STATIC_EQUILIBRIUM_ALGORITHM_PP);


  RVector2 com_LB, com_UB;
  com_LB(0) = p.col(0).minCoeff()-X_MARG;
  com_UB(0) = p.col(0).maxCoeff()+X_MARG;
  com_LB(1) = p.col(1).minCoeff()-Y_MARG;
  com_UB(1) = p.col(1).maxCoeff()+Y_MARG;
  MatrixXX com(N_POINTS,2);
  for(int i=0; i<N_POINTS; i++)
  {
    uniform(com_LB, com_UB, com.row(i));
    double robustness = solver.checkRobustEquilibrium(com.row(i), 0.0);
    printf("CoM position (%.3f, %.3f) robustness %.3f\n", com(i,0), com(i,1), robustness);
  }

  return ret;
}

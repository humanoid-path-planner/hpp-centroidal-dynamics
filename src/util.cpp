#ifndef _ROBUST_EQUILIBRIUM_LIB_CONFIG_HH
#define _ROBUST_EQUILIBRIUM_LIB_CONFIG_HH

#include <Eigen/Dense>
#include <Eigen/src/Core/util/Macros.h>

#include "cdd/cddmp.h"
#include "cdd/setoper.h"
#include "cdd/cddtypes.h"
#include "cdd/cdd.h"

namespace robust_equilibrium
{

//#define USE_FLOAT 1;
#ifdef USE_FLOAT
typedef float value_type;
#else
typedef double value_type;
#endif

typedef Eigen::Matrix <value_type, 2, 1>                                            Vector2;
typedef Eigen::Matrix <value_type, 1, 2>                                            RVector2;
typedef Eigen::Matrix <value_type, 3, 1>                                            Vector3;
typedef Eigen::Matrix <value_type, 1, 3>                                            RVector3;
typedef Eigen::Matrix <value_type, Eigen::Dynamic, 1>                               VectorX;
typedef Eigen::Matrix <value_type, 1, Eigen::Dynamic>                               RVectorX;
typedef Eigen::Matrix <value_type, 3, 3, Eigen::RowMajor>                           Rotation;
typedef Eigen::Matrix <value_type, 3, 3, Eigen::RowMajor>                           Matrix3;
typedef Eigen::Matrix <value_type, Eigen::Dynamic, 3, Eigen::RowMajor>              MatrixX3;
typedef Eigen::Matrix <value_type, 4, 3, Eigen::RowMajor>                           Matrix43;
typedef Eigen::Matrix <value_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXX;

//define Eigen ref if available
#if EIGEN_VERSION_AT_LEAST(3,2,0)
typedef Eigen::Ref<Vector2>     Ref_vector2;
typedef Eigen::Ref<Vector3>     Ref_vector3;
typedef Eigen::Ref<VectorX>     Ref_vectorX;
typedef Eigen::Ref<Rotation>    Ref_rotation;
typedef Eigen::Ref<MatrixX3>    Ref_matrixX3;
typedef Eigen::Ref<Matrix43>    Ref_matrix43;
typedef Eigen::Ref<MatrixXX>    Ref_matrixXX;

typedef const Eigen::Ref<const Vector2>     & Cref_vector2;
typedef const Eigen::Ref<const Vector3>     & Cref_vector3;
typedef const Eigen::Ref<const VectorX>     & Cref_vectorX;
typedef const Eigen::Ref<const Rotation>    & Cref_rotation;
typedef const Eigen::Ref<const MatrixX3>    & Cref_matrixX3;
typedef const Eigen::Ref<const Matrix43>    & Cref_matrix43;
typedef const Eigen::Ref<const MatrixXX>    & Cref_matrixXX;
#else
typedef vector2_t   & Ref_vector2;
typedef vector3_t   & Ref_vector3;
typedef vector_t    & Ref_vector;
typedef rotation_t  & Ref_rotation;
typedef T_rotation_t& Ref_matrixX3;
typedef T_rotation_t& Ref_matrix43;
typedef matrix_t    & Ref_matrixXX;

typedef const vector2_t   & Cref_vector2;
typedef const vector3_t   & Cref_vector3;
typedef const vector_t    & Cref_vector;
typedef const rotation_t  & Cref_rotation;
typedef const T_rotation_t& Cref_matrixX3;
typedef const T_rotation_t& Cref_matrix43;
typedef const matrix_t    & Cref_matrixXX;
#endif

typedef Eigen::AngleAxis<value_type> angle_axis_t;


dd_MatrixPtr cone_span_eigen_to_cdd(Cref_matrixXX input)
{
  dd_debug = false;
  dd_MatrixPtr M=NULL;
  dd_rowrange i;
  dd_colrange j;
  dd_rowrange m_input = (dd_rowrange)(input.rows());
  dd_colrange d_input = (dd_colrange)(input.cols()+1);
  dd_RepresentationType rep=dd_Generator;
  mytype value;
  dd_NumberType NT = dd_Real;
  dd_init(value);

  M=dd_CreateMatrix(m_input, d_input);
  M->representation=rep;
  M->numbtype=NT;

  for (i = 0; i < input.rows(); i++)
  {
    dd_set_d(value, 0);
    dd_set(M->matrix[i][0],value);
    for (j = 1; j < d_input; j++)
    {
      dd_set_d(value, input(i,j-1));
      dd_set(M->matrix[i][j],value);
    }
  }
  dd_clear(value);
  return M;
}


void init_library()
{
  dd_set_global_constants();dd_debug = false;
}

void release_library()
{
  //dd_free_global_constants();
}

void uniform(Cref_matrixXX lower_bounds, Cref_matrixXX upper_bounds, Ref_matrixXX out)
{

  assert(lower_bounds.rows()==out.rows());
  assert(upper_bounds.rows()==out.rows());
  assert(lower_bounds.cols()==out.cols());
  assert(upper_bounds.cols()==out.cols());
  for(int i=0; i<out.rows(); i++)
    for(int j=0; j<out.cols(); j++)
      out(i,j) = (rand()/ value_type(RAND_MAX))*(upper_bounds(i,j)-lower_bounds(i,j)) + lower_bounds(i,j);
}

void euler_matrix(double roll, double pitch, double yaw, Ref_rotation R)
{
  const int i = 0;
  const int j = 1;
  const int k = 2;
  double si = sin(roll);
  double sj = sin(pitch);
  double sk = sin(yaw);
  double ci = cos(roll);
  double cj = cos(pitch);
  double ck = cos(yaw);
  double cc = ci*ck;
  double cs = ci*sk;
  double sc = si*ck;
  double ss = si*sk;

  R(i, i) = cj*ck;
  R(i, j) = sj*sc-cs;
  R(i, k) = sj*cc+ss;
  R(j, i) = cj*sk;
  R(j, j) = sj*ss+cc;
  R(j, k) = sj*cs-sc;
  R(k, i) = -sj;
  R(k, j) = cj*si;
  R(k, k) = cj*ci;
//  R = (angle_axis_t(roll, Vector3::UnitX())
//       * angle_axis_t(pitch, Vector3::UnitY())
//       * angle_axis_t(yaw, Vector3::UnitZ())).toRotationMatrix();
}

bool generate_rectangle_contacts(double lx, double ly, Cref_vector3 pos, Cref_vector3 rpy,
                                 Ref_matrix43 p, Ref_matrix43 N)
{
  // compute rotation matrix
  Rotation R;
  euler_matrix(rpy(0), rpy(1), rpy(2), R);
  // contact points in local frame
  p << lx,  ly, 0,
       lx, -ly, 0,
      -lx, -ly, 0,
      -lx,  ly, 0;
  // contact points in world frame
  p.row(0) = pos + (R*p.row(0).transpose());
  p.row(1) = pos + (R*p.row(1).transpose());
  p.row(2) = pos + (R*p.row(2).transpose());
  p.row(3) = pos + (R*p.row(3).transpose());
  // normal direction in local frame
  RVector3 n;
  n << 0, 0, 1;
  // normal directions in world frame
  n = (R*n.transpose()).transpose();
  N << n, n, n, n;
  return true;
}

Rotation crossMatrix(Cref_vector3 x)
{
    Rotation res = Rotation::Zero();
    res(0,1) = - x(2); res(0,2) =   x(1);
    res(1,0) =   x(2); res(1,2) = - x(0);
    res(2,0) = - x(1); res(2,1) =   x(0);
    return res;
}

} //namespace robust_equilibrium

#endif //_ROBUST_EQUILIBRIUM_LIB_CONFIG_HH

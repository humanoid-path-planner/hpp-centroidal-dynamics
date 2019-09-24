/*
 * Copyright 2015, LAAS-CNRS
 * Author: Andrea Del Prete
 */

#include <ctime>
#include <hpp/centroidal-dynamics/util.hh>
#include <ctime>

namespace centroidal_dynamics {

dd_MatrixPtr cone_span_eigen_to_cdd(Cref_matrixXX input, const bool canonicalize) {
  dd_debug = false;
  dd_MatrixPtr M = NULL;
  dd_rowrange i;
  dd_colrange j;
  dd_rowrange m_input = (dd_rowrange)(input.rows());
  dd_colrange d_input = (dd_colrange)(input.cols() + 1);
  dd_RepresentationType rep = dd_Generator;
  mytype value;
  dd_NumberType NT = dd_Real;
  dd_init(value);
  M = dd_CreateMatrix(m_input, d_input);
  M->representation = rep;
  M->numbtype = NT;
  for (i = 0; i < input.rows(); i++) {
    dd_set_d(value, 0);
    dd_set(M->matrix[i][0], value);
    for (j = 1; j < d_input; j++) {
      dd_set_d(value, input(i, j - 1));
      dd_set(M->matrix[i][j], value);
    }
  }
  dd_clear(value);
  if (canonicalize) {
    dd_ErrorType error = dd_NoError;
    dd_rowset redset, impl_linset;
    dd_rowindex newpos;
    dd_MatrixCanonicalize(&M, &impl_linset, &redset, &newpos, &error);
    set_free(redset);
    set_free(impl_linset);
    free(newpos);
  }
  return M;
}

void init_cdd_library() {
  dd_set_global_constants();
  dd_debug = false;
}

void release_cdd_library() {
  // dd_free_global_constants();
}

void uniform3(Cref_vector3 lower_bounds, Cref_vector3 upper_bounds, Ref_vector3 out) {
  assert(lower_bounds.rows() == out.rows());
  assert(upper_bounds.rows() == out.rows());
  assert(lower_bounds.cols() == out.cols());
  assert(upper_bounds.cols() == out.cols());
  for (int i = 0; i < out.rows(); i++)
    for (int j = 0; j < out.cols(); j++)
      out(i, j) = (rand() / value_type(RAND_MAX)) * (upper_bounds(i, j) - lower_bounds(i, j)) + lower_bounds(i, j);
}

void uniform(Cref_matrixXX lower_bounds, Cref_matrixXX upper_bounds, Ref_matrixXX out) {
  assert(lower_bounds.rows() == out.rows());
  assert(upper_bounds.rows() == out.rows());
  assert(lower_bounds.cols() == out.cols());
  assert(upper_bounds.cols() == out.cols());
  for (int i = 0; i < out.rows(); i++)
    for (int j = 0; j < out.cols(); j++)
      out(i, j) = (rand() / value_type(RAND_MAX)) * (upper_bounds(i, j) - lower_bounds(i, j)) + lower_bounds(i, j);
}

void euler_matrix(double roll, double pitch, double yaw, Ref_rotation R) {
  const int i = 0;
  const int j = 1;
  const int k = 2;
  double si = sin(roll);
  double sj = sin(pitch);
  double sk = sin(yaw);
  double ci = cos(roll);
  double cj = cos(pitch);
  double ck = cos(yaw);
  double cc = ci * ck;
  double cs = ci * sk;
  double sc = si * ck;
  double ss = si * sk;

  R(i, i) = cj * ck;
  R(i, j) = sj * sc - cs;
  R(i, k) = sj * cc + ss;
  R(j, i) = cj * sk;
  R(j, j) = sj * ss + cc;
  R(j, k) = sj * cs - sc;
  R(k, i) = -sj;
  R(k, j) = cj * si;
  R(k, k) = cj * ci;
  //  R = (angle_axis_t(roll, Vector3::UnitX())
  //       * angle_axis_t(pitch, Vector3::UnitY())
  //       * angle_axis_t(yaw, Vector3::UnitZ())).toRotationMatrix();
}

bool generate_rectangle_contacts(double lx, double ly, Cref_vector3 pos, Cref_vector3 rpy, Ref_matrix43 p,
                                 Ref_matrix43 N) {
  // compute rotation matrix
  Rotation R;
  euler_matrix(rpy(0), rpy(1), rpy(2), R);
  // contact points in local frame
  p << lx, ly, 0, lx, -ly, 0, -lx, -ly, 0, -lx, ly, 0;
  // contact points in world frame
  p.row(0) = pos + (R * p.row(0).transpose());
  p.row(1) = pos + (R * p.row(1).transpose());
  p.row(2) = pos + (R * p.row(2).transpose());
  p.row(3) = pos + (R * p.row(3).transpose());
  // normal direction in local frame
  RVector3 n;
  n << 0, 0, 1;
  // normal directions in world frame
  n = (R * n.transpose()).transpose();
  N << n, n, n, n;
  return true;
}

Rotation crossMatrix(Cref_vector3 x) {
  Rotation res = Rotation::Zero();
  res(0, 1) = -x(2);
  res(0, 2) = x(1);
  res(1, 0) = x(2);
  res(1, 2) = -x(0);
  res(2, 0) = -x(1);
  res(2, 1) = x(0);
  return res;
}

std::string getDateAndTimeAsString() {
  time_t rawtime;
  struct tm* timeinfo;
  char buffer[80];
  time(&rawtime);
  timeinfo = localtime(&rawtime);
  strftime(buffer, 80, "%Y%m%d_%I%M%S", timeinfo);
  return std::string(buffer);
}
/*
int fact(const int n)
{
    assert(n>=0);
    int res = 1;
    for (int i=2 ; i <= n ; ++i)
       res *= i;
    return res;
}

int choosenk(const int n, const int k)
{
    return fact(n) / (fact(k) * fact(n - k));
}*/

/* is this faster ?
value_type choosenk(const int n, const int k)
{
    if(k>n/2)
        return nchoosek(n,n-k);
    else if(k==1)
        return n;
    else
    {
        double c = 1;
        for(int i = 1;i<=k;i++)
            c *= (((double)n-k+i)/((double)i));
        return std::round(c);
    }
}*/

value_type nchoosek(const int n, const int k) {
  if (k > n / 2)
    return nchoosek(n, n - k);
  else if (k == 1)
    return n;
  else {
    value_type c = 1;
    for (int i = 1; i <= k; i++) c *= (((value_type)n - k + i) / ((value_type)i));
    return round(c);
  }
}

}  // namespace centroidal_dynamics

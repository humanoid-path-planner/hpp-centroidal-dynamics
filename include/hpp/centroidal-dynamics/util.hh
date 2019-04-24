/*
 * Copyright 2015, LAAS-CNRS
 * Author: Andrea Del Prete
 */
#ifndef HPP_CENTROIDAL_DYNAMICS_UTIL_HH
#define HPP_CENTROIDAL_DYNAMICS_UTIL_HH

#include <iostream>
#include <fstream>
#include <cmath>
#include <cassert>

#include <Eigen/Dense>
#include <Eigen/src/Core/util/Macros.h>

#include "cddmp.h"
#include "setoper.h"
#include "cddtypes.h"
#include "cdd.h"

namespace centroidal_dynamics
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
  typedef Eigen::Matrix <value_type, 6, 1>                                            Vector6;
  typedef Eigen::Matrix <value_type, Eigen::Dynamic, 1>                               VectorX;
  typedef Eigen::Matrix <value_type, 1, Eigen::Dynamic>                               RVectorX;
  typedef Eigen::Matrix <value_type, 3, 3, Eigen::RowMajor>                           Rotation;
  typedef Eigen::Matrix <value_type, Eigen::Dynamic, 2, Eigen::RowMajor>              MatrixX2;
  typedef Eigen::Matrix <value_type, 3, 3, Eigen::RowMajor>                           Matrix3;
  typedef Eigen::Matrix <value_type, Eigen::Dynamic, 3, Eigen::RowMajor>              MatrixX3;
  typedef Eigen::Matrix <value_type, 3, Eigen::Dynamic, Eigen::RowMajor>              Matrix3X;
  typedef Eigen::Matrix <value_type, 4, 3, Eigen::RowMajor>                           Matrix43;
  typedef Eigen::Matrix <value_type, 6, Eigen::Dynamic, Eigen::RowMajor>              Matrix6X;
  typedef Eigen::Matrix <value_type, 6, 2, Eigen::RowMajor>                           Matrix62;
  typedef Eigen::Matrix <value_type, 6, 3, Eigen::RowMajor>                           Matrix63;
  typedef Eigen::Matrix <value_type, Eigen::Dynamic, 6, Eigen::RowMajor>              MatrixX6;
  typedef Eigen::Matrix <value_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXX;

  typedef Eigen::Ref<Vector2>     Ref_vector2;
  typedef Eigen::Ref<Vector3>     Ref_vector3;
  typedef Eigen::Ref<VectorX>     Ref_vectorX;
  typedef Eigen::Ref<Rotation>    Ref_rotation;
  typedef Eigen::Ref<MatrixX3>    Ref_matrixX3;
  typedef Eigen::Ref<Matrix43>    Ref_matrix43;
  typedef Eigen::Ref<Matrix6X>    Ref_matrix6X;
  typedef Eigen::Ref<MatrixXX>    Ref_matrixXX;

  typedef const Eigen::Ref<const Vector2>     & Cref_vector2;
  typedef const Eigen::Ref<const Vector3>     & Cref_vector3;
  typedef const Eigen::Ref<const Vector6>     & Cref_vector6;
  typedef const Eigen::Ref<const VectorX>     & Cref_vectorX;
  typedef const Eigen::Ref<const Rotation>    & Cref_rotation;
  typedef const Eigen::Ref<const MatrixX3>    & Cref_matrixX3;
  typedef const Eigen::Ref<const Matrix43>    & Cref_matrix43;
  typedef const Eigen::Ref<const Matrix6X>    & Cref_matrix6X;
  typedef const Eigen::Ref<const Matrix63>    & Cref_matrix63;
  typedef const Eigen::Ref<const MatrixXX>    & Cref_matrixXX;

  /**Column major definitions for compatibility with classical eigen use**/
  typedef Eigen::Matrix <value_type, Eigen::Dynamic, 3, Eigen::ColMajor>              MatrixX3ColMajor;
  typedef Eigen::Matrix <value_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> MatrixXXColMajor;
  typedef const Eigen::Ref<const MatrixX3ColMajor>  & Cref_matrixX3ColMajor;
  typedef Eigen::Ref<MatrixXXColMajor>              &  ref_matrixXXColMajor;

  /**
   * Write the specified matrix to a binary file with the specified name.
   */
  template<class Matrix>
  bool writeMatrixToFile(const std::string &filename, const Matrix& matrix)
  {
    std::ofstream out(filename.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
    if(!out.is_open())
      return false;
    typename Matrix::Index rows=matrix.rows(), cols=matrix.cols();
    out.write((char*) (&rows), sizeof(typename Matrix::Index));
    out.write((char*) (&cols), sizeof(typename Matrix::Index));
    out.write((char*) matrix.data(), rows*cols*sizeof(typename Matrix::Scalar) );
    out.close();
    return true;
  }

  /**
   * Read a matrix from the specified input binary file.
   */
  template<class Matrix>
  bool readMatrixFromFile(const std::string &filename, Matrix& matrix)
  {
    std::ifstream in(filename.c_str(), std::ios::in | std::ios::binary);
    if(!in.is_open())
      return false;
    typename Matrix::Index rows=0, cols=0;
    in.read((char*) (&rows),sizeof(typename Matrix::Index));
    in.read((char*) (&cols),sizeof(typename Matrix::Index));
    matrix.resize(rows, cols);
    in.read( (char *) matrix.data() , rows*cols*sizeof(typename Matrix::Scalar) );
    in.close();
    return true;
  }

  /**
 * Convert the specified list of rays from Eigen to cdd format.
 * @param input The mXn input Eigen matrix contains m rays of dimension n.
 * @param bool whether to remove redundant inequalities
 * @return The mX(n+1) output cdd matrix, which contains an additional column,
 * the first one, with all zeros.
 */
  dd_MatrixPtr cone_span_eigen_to_cdd(Cref_matrixXX input, const bool canonicalize=false);

  /**
 * Compute the cross-product skew-symmetric matrix associated to the specified vector.
 */
  Rotation crossMatrix(Cref_vector3 x);


  void init_cdd_library();

  void release_cdd_library();

  // in some distribution the conversion Ref_matrixXX to Ref_vector3 does not compile
  void uniform3(Cref_vector3 lower_bounds, Cref_vector3 upper_bounds, Ref_vector3 out);
  void uniform(Cref_matrixXX lower_bounds, Cref_matrixXX upper_bounds, Ref_matrixXX out);

  void euler_matrix(double roll, double pitch, double yaw, Ref_rotation R);

  bool generate_rectangle_contacts(double lx, double ly, Cref_vector3 pos, Cref_vector3 rpy,
                                   Ref_matrix43 p, Ref_matrix43 N);

  std::string getDateAndTimeAsString();

  /**
 * Computes a binomal coefficient
 * @return  n!/((n–k)! k!).
 */
  value_type nchoosek(const int n, const int k);

  template < typename DerivedV, typename DerivedU>
  void doCombs(Eigen::Matrix<typename DerivedU::Scalar,1,Eigen::Dynamic>& running,
               int& running_i, int& running_j, Eigen::PlainObjectBase<DerivedU> & U, const Eigen::MatrixBase<DerivedV> & V,
               int offset, int k)
  {
      int N = (int)(V.size());
      if(k==0)
      {
          U.row(running_i) = running;
          running_i++;
          return;
      }
      for (int i = offset; i <= N - k; ++i)
      {
          running(running_j) = V(i);
          running_j++;
          doCombs(running, running_i, running_j, U, V, i+1,k-1);
          running_j--;
      }
  }

  /**
 * Computes a matrix C containing all possible combinations of the elements of vector v taken k at a time.
 * Matrix C has k columns and n!/((n–k)! k!) rows, where n is length(v).
 * @param V  n-long vector of elements
 * @param k  size of sub-set to consider
 * @param U  result matrix
 * @return nchoosek by k long matrix where each row is a unique k-size
 * the first one, with all zeros.
 */
template <typename DerivedV, typename DerivedU>
void nchoosek(const Eigen::MatrixBase<DerivedV>& V, const int k, Eigen::PlainObjectBase<DerivedU>& U) {
  using namespace Eigen;
  if (V.size() == 0) {
    U.resize(0, k);
    return;
  }
  assert((V.cols() == 1 || V.rows() == 1) && "V must be a vector");
  U.resize(nchoosek((int)(V.size()), k), k);
  int running_i = 0;
  int running_j = 0;
  Matrix<typename DerivedU::Scalar, 1, Dynamic> running(1, k);
  doCombs(running, running_i, running_j, U, V, 0, k);
}
}  // namespace centroidal_dynamics

#endif  // HPP_CENTROIDAL_DYNAMICS_UTIL_HH

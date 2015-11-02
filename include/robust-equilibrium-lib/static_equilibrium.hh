/*
 * Copyright 2015, LAAS-CNRS
 * Author: Andrea Del Prete
 */

#ifndef ROBUST_EQUILIBRIUM_LIB_STATIC_EQUILIBRIUM_H
#define ROBUST_EQUILIBRIUM_LIB_STATIC_EQUILIBRIUM_H

#include <Eigen/Dense>
#include <robust-equilibrium-lib/config.hh>
#include <robust-equilibrium-lib/util.hh>
#include <robust-equilibrium-lib/solver_LP_abstract.hh>

namespace robust_equilibrium
{

enum ROBUST_EQUILIBRIUM_DLLAPI StaticEquilibriumAlgorithm
{
  STATIC_EQUILIBRIUM_ALGORITHM_LP,
  STATIC_EQUILIBRIUM_ALGORITHM_LP2,
  STATIC_EQUILIBRIUM_ALGORITHM_DLP,
  STATIC_EQUILIBRIUM_ALGORITHM_PP,
  STATIC_EQUILIBRIUM_ALGORITHM_IP,
  STATIC_EQUILIBRIUM_ALGORITHM_DIP
};

class ROBUST_EQUILIBRIUM_DLLAPI StaticEquilibrium
{
private:
  static bool m_is_cdd_initialized;

  std::string                 m_name;
  StaticEquilibriumAlgorithm  m_algorithm;
  Solver_LP_abstract*         m_solver;
  SolverLP                    m_solver_type;

  unsigned int m_generatorsPerContact;
  double m_mass;
  Vector3 m_gravity;

   /** Tangent directions for all contacts (numberOfContacts*generatorsPerContact X 3) */
  MatrixX3 m_T1, m_T2;
  /** Matrix mapping contact forces to gravito-inertial wrench (6 X 3*numberOfContacts) */
  Matrix6X m_A;
  /** Lists of contact generators (3 X numberOfContacts*generatorsPerContact) */
  Matrix3X m_G;
  /** Gravito-inertial wrench generators (6 X numberOfContacts*generatorsPerContact) */
  Matrix6X m_G_centr;
  /** Inequality matrix defining the gravito-inertial wrench cone H w <= h */
  MatrixXX m_H;
  /** Inequality vector defining the gravito-inertial wrench cone H w <= h */
  VectorX m_h;
  /** Inequality matrix defining the CoM support polygon HD com + Hd <= h */
  MatrixX2 m_HD;
  VectorX  m_Hd;
  /** Matrix and vector mapping 2d com position to GIW */
  Matrix62 m_D;
  Vector6 m_d;

  bool computePolytopeProjection(Cref_matrix6X v);

public:
  StaticEquilibrium(std::string name, double mass, unsigned int generatorsPerContact, SolverLP solver_type);

  std::string getName(){ return m_name; }

  bool setNewContacts(Cref_matrixX3 contactPoints, Cref_matrixX3 contactNormals,
                      Cref_vectorX frictionCoefficients, StaticEquilibriumAlgorithm alg);

  bool computeEquilibriumRobustness(Cref_vector2 com, double &robustness);

  bool checkRobustEquilibrium(Cref_vector2 com, bool &equilibrium, double e_max=0.0);

  /** Compute the extremum CoM position over the line a*x + a0 that is in robust equilibrium.
   * This is equivalent to following the following LP:
   *     find          c, b
   *     maximize      c
   *     subject to    G b = D (a c + a0) + d
   *                   b  >= b0
   *   where
   *     b         are the m coefficients of the contact force generators (f = G b)
   *     b0        is the m-dimensional vector of identical values that are proportional to e_max
   *     c         is the 1d line parameter
   *     G         is the 6xm matrix whose columns are the gravito-inertial wrench generators
   *     D         is the 6x2 matrix mapping the CoM position in gravito-inertial wrench
   *     d         is the 6d vector containing the gravity part of the gravito-inertial wrench
   * @param a 2d vector representing the line direction
   * @param a0 2d vector representing an arbitrary point over the line
   * @param e_max Desired robustness in terms of the maximum force error tolerated by the system
   * @return True if the operation succeeded, false otherwise.
  */
  bool findExtremumOverLine(Cref_vector2 a, Cref_vector2 a0, double e_max, Ref_vector2 com);

  bool findExtremumInDirection(Cref_vector2 direction, Ref_vector2 com, double e_max=0.0);

};

} // end namespace robust_equilibrium

#endif

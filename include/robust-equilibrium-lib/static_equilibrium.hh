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

  StaticEquilibriumAlgorithm  m_algorithm;
  Solver_LP_abstract*         m_solver;

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
  StaticEquilibrium(double mass, unsigned int generatorsPerContact, SolverLP solver_type);

  bool setNewContacts(Cref_matrixX3 contactPoints, Cref_matrixX3 contactNormals,
                      Cref_vectorX frictionCoefficients, StaticEquilibriumAlgorithm alg);

  double computeEquilibriumRobustness(Cref_vector2 com);

  bool checkRobustEquilibrium(Cref_vector2 com, double e_max=0.0);

  double findExtremumOverLine(Cref_vector2 a, double b, double e_max=0.0);

  double findExtremumInDirection(Cref_vector2 direction, double e_max=0.0);

};

} // end namespace robust_equilibrium

#endif

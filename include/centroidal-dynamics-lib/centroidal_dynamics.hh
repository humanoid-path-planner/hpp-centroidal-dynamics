/*
 * Copyright 2015, LAAS-CNRS
 * Author: Andrea Del Prete
 */

#ifndef CENTROIDAL_DYNAMICS_LIB_STATIC_EQUILIBRIUM_H
#define CENTROIDAL_DYNAMICS_LIB_STATIC_EQUILIBRIUM_H

#include <Eigen/Dense>
#include <centroidal-dynamics-lib/config.hh>
#include <centroidal-dynamics-lib/util.hh>
#include <centroidal-dynamics-lib/solver_LP_abstract.hh>

namespace centroidal_dynamics
{

enum CENTROIDAL_DYNAMICS_DLLAPI EquilibriumAlgorithm
{
  EQUILIBRIUM_ALGORITHM_LP,  /// primal LP formulation
  EQUILIBRIUM_ALGORITHM_LP2, /// another primal LP formulation
  EQUILIBRIUM_ALGORITHM_DLP, /// dual LP formulation
  EQUILIBRIUM_ALGORITHM_PP,  /// polytope projection algorithm
  EQUILIBRIUM_ALGORITHM_IP,  /// incremental projection algorithm based on primal LP formulation
  EQUILIBRIUM_ALGORITHM_DIP  /// incremental projection algorithm based on dual LP formulation
};

class CENTROIDAL_DYNAMICS_DLLAPI Equilibrium
{
private:
  static bool m_is_cdd_initialized;   /// true if cdd lib has been initialized, false otherwise

  std::string                 m_name;         /// name of this object
  EquilibriumAlgorithm  m_algorithm;    /// current algorithm used
  Solver_LP_abstract*         m_solver;       /// LP solver
  SolverLP                    m_solver_type;  /// type of LP solver

  unsigned int  m_generatorsPerContact; /// number of generators to approximate the friction cone per contact point
  double        m_mass;                 /// mass of the system
  Vector3       m_gravity;              /// gravity vector

  /** Gravito-inertial wrench generators (6 X numberOfContacts*generatorsPerContact) */
  Matrix6X m_G_centr;

  /** Inequality matrix and vector defining the gravito-inertial wrench cone H w <= h */
  MatrixXX m_H;
  VectorX m_h;
  /** False if a numerical instability appeared in the computation H and h*/
  bool m_is_cdd_stable;
  /** EQUILIBRIUM_ALGORITHM_PP: If double description fails,
    * indicate the max number of attempts to compute the cone by introducing
    * a small pertubation of the system */
  const unsigned max_num_cdd_trials;
  /** whether to remove redundant inequalities when computing double description matrices*/
  const bool canonicalize_cdd_matrix;

  /** Inequality matrix and vector defining the CoM support polygon HD com + Hd <= h */
  MatrixX3 m_HD;
  VectorX  m_Hd;

  /** Matrix and vector mapping 2d com position to GIW */
  Matrix63 m_D;
  Vector6 m_d;

  /** Coefficient used for converting the robustness measure in Newtons */
  double m_b0_to_emax_coefficient;

  bool computePolytopeProjection(Cref_matrix6X v);
  bool computeGenerators(Cref_matrixX3 contactPoints, Cref_matrixX3 contactNormals,
                         double frictionCoefficient, const bool perturbate = false);

  /**
   * @brief Given the smallest coefficient of the contact force generators it computes
   * the minimum norm of force error necessary to have a contact force on
   * the associated friction cone boundaries.
   * @param b0 Minimum coefficient of the contact force generators.
   * @return Minimum norm of the force error necessary to result in a contact force being
   * on the friction cone boundaries.
   */
  double convert_b0_to_emax(double b0);

  double convert_emax_to_b0(double emax);

public:

  /**
   * @brief Equilibrium constructor.
   * @param name Name of the object.
   * @param mass Mass of the system for which to test equilibrium.
   * @param generatorsPerContact Number of generators used to approximate the friction cone per contact point.
   * @param solver_type Type of LP solver to use.
   * @param useWarmStart Whether the LP solver can warm start the resolution.
   * @param max_num_cdd_trials indicate the max number of attempts to compute the cone by introducing
   * @param canonicalize_cdd_matrix whether to remove redundant inequalities when computing double description matrices
   * a small pertubation of the system
   */
  Equilibrium(std::string name, double mass, unsigned int generatorsPerContact,
                    SolverLP solver_type = SOLVER_LP_QPOASES, bool useWarmStart=true, const unsigned int max_num_cdd_trials=0,
                    const bool canonicalize_cdd_matrix=false);

  /**
   * @brief Returns the useWarmStart flag.
   * @return True if the LP solver is allowed to use warm start, false otherwise.
   */
  bool useWarmStart(){ return m_solver->getUseWarmStart(); }

  /**
   * @brief Specifies whether the LP solver is allowed to use warm start.
   * @param uws True if the LP solver is allowed to use warm start, false otherwise.
   */
  void setUseWarmStart(bool uws){ m_solver->setUseWarmStart(uws); }

  /**
   * @brief Get the name of this object.
   * @return The name of this object.
   */
  std::string getName(){ return m_name; }

  EquilibriumAlgorithm getAlgorithm(){ return m_algorithm; }

  /**
   * @brief Specify a new set of contacts.
   * All 3d vectors are expressed in a reference frame having the z axis aligned with gravity.
   * In other words the gravity vecotr is (0, 0, -9.81).
   * @param contactPoints List of N 3d contact points as an Nx3 matrix.
   * @param contactNormals List of N 3d contact normal directions as an Nx3 matrix.
   * @param frictionCoefficient The contact friction coefficient.
   * @param alg Algorithm to use for testing equilibrium.
   * @return True if the operation succeeded, false otherwise.
   */
  bool setNewContacts(Cref_matrixX3 contactPoints, Cref_matrixX3 contactNormals,
                      double frictionCoefficient, EquilibriumAlgorithm alg);

  /**
   * @brief Compute a measure of the robustness of the equilibrium of the specified com position.
   * This amounts to solving the following LP:
   *       find          b, b0
   *       maximize      b0
   *       subject to    G b = D c + d
   *                     b >= b0
   *  where:
   *     b         are the coefficient of the contact force generators (f = G b)
   *     b0        is a parameter proportional to the robustness measure
   *     c         is the specified CoM position
   *     G         is the 6xm matrix whose columns are the gravito-inertial wrench generators
   *     D         is the 6x3 matrix mapping the CoM position in gravito-inertial wrench
   *     d         is the 6d vector containing the gravity part of the gravito-inertial wrench
   * @param com The 3d center of mass position to test.
   * @param robustness The computed measure of robustness.
   * @return The status of the LP solver.
   * @note If the system is in force closure the status will be LP_STATUS_UNBOUNDED, meaning that the
   * system can reach infinite robustness. This is due to the fact that we are not considering
   * any upper limit for the friction cones.
   */
  LP_status computeEquilibriumRobustness(Cref_vector3 com, double &robustness);

  /**
   * @brief Compute a measure of the robustness of the equilibrium of the specified com position.
   * This amounts to solving the following LP:
   *       find          b, b0
   *       maximize      b0
   *       subject to    G b = D c + d
   *                     b >= b0
   *  where:
   *     b         are the coefficient of the contact force generators (f = G b)
   *     b0        is a parameter proportional to the robustness measure
   *     c         is the specified CoM position
   *     G         is the 6xm matrix whose columns are the gravito-inertial wrench generators
   *     D         is the 6x3 matrix mapping the CoM position in gravito-inertial wrench
   *     d         is the 6d vector containing the gravity part of the gravito-inertial wrench
   * @param com The 3d center of mass position to test.
   * @param acc The 3d acceleration of the CoM.
   * @param robustness The computed measure of robustness.
   * @return The status of the LP solver.
   * @note If the system is in force closure the status will be LP_STATUS_UNBOUNDED, meaning that the
   * system can reach infinite robustness. This is due to the fact that we are not considering
   * any upper limit for the friction cones.
   */
  LP_status computeEquilibriumRobustness(Cref_vector3 com, Cref_vector3 acc, double &robustness);

  /**
   * @brief Check whether the specified com position is in robust equilibrium.
   * This amounts to solving the following feasibility LP:
   *       find          b
   *       minimize      1
   *       subject to    G b = D c + d
   *                     b >= b0
   *  where:
   *     b         are the coefficient of the contact force generators (f = G b)
   *     b0        is a parameter proportional to the specified robustness measure
   *     c         is the specified CoM position
   *     G         is the 6xm matrix whose columns are the gravito-inertial wrench generators
   *     D         is the 6x3 matrix mapping the CoM position in gravito-inertial wrench
   *     d         is the 6d vector containing the gravity part of the gravito-inertial wrench
   * @param com The 3d center of mass position to test.
   * @param equilibrium True if com is in robust equilibrium, false otherwise.
   * @param e_max Desired robustness level.
   * @return The status of the LP solver.
   */
  LP_status checkRobustEquilibrium(Cref_vector3 com, bool &equilibrium, double e_max=0.0);

  /**
   * @brief Compute the extremum CoM position over the line a*x + a0 that is in robust equilibrium.
   * This amounts to solving the following LP:
   *     find          c, b
   *     maximize      c
   *     subject to    G b = D (a c + a0) + d
   *                   b  >= b0
   *   where:
   *     b         are the m coefficients of the contact force generators (f = G b)
   *     b0        is an m-dimensional vector of identical values that are proportional to e_max
   *     c         is the 1d line parameter
   *     G         is the 6xm matrix whose columns are the gravito-inertial wrench generators
   *     D         is the 6x3 matrix mapping the CoM position in gravito-inertial wrench
   *     d         is the 6d vector containing the gravity part of the gravito-inertial wrench
   * @param a 2d vector representing the line direction
   * @param a0 2d vector representing an arbitrary point over the line
   * @param e_max Desired robustness in terms of the maximum force error tolerated by the system
   * @return The status of the LP solver.
   * @note If the system is in force closure the status will be LP_STATUS_UNBOUNDED, meaning that the
   * system can reach infinite robustness. This is due to the fact that we are not considering
   * any upper limit for the friction cones.
  */
  LP_status findExtremumOverLine(Cref_vector3 a, Cref_vector3 a0, double e_max, Ref_vector3 com);

  /**
   * @brief Find the extremum com position that is in robust equilibrium in the specified direction.
   * This amounts to solving the following LP:
   *     find          c, b
   *     maximize      a^T c
   *     subject to    G b = D c + d
   *                   b  >= b0
   * where:
   *     a         is the specified 2d direction
   *     b         are the m coefficients of the contact force generators (f = G b)
   *     b0        is an m-dimensional vector of identical values that are proportional to e_max
   *     c         is the 3d com position
   *     G         is the 6xm matrix whose columns are the gravito-inertial wrench generators
   *     D         is the 6x3 matrix mapping the CoM position in gravito-inertial wrench
   *     d         is the 6d vector containing the gravity part of the gravito-inertial wrench
   * @param direction Desired 3d direction.
   * @param com Output 3d com position.
   * @param e_max Desired robustness level.
   * @return The status of the LP solver.
   * @note If the system is in force closure the status will be LP_STATUS_UNBOUNDED, meaning that the
   * system can reach infinite robustness. This is due to the fact that we are not considering
   * any upper limit for the friction cones.
   */
  LP_status findExtremumInDirection(Cref_vector3 direction, Ref_vector3 com, double e_max=0.0);

  /**
   * @brief Retrieve the inequalities that define the admissible wrenchs
   * for the current contact set.
   * @param H reference to the H matrix to initialize
   * @param h reference to the h vector to initialize
   * @return The status of the inequalities. If the inequalities are not defined
   * due to numerical instabilities, will send appropriate error message,
   * and return LP_STATUS_ERROR. If they are not defined because no
   * contact has been defined, will return LP_STATUS_INFEASIBLE
   */
  LP_status getPolytopeInequalities(MatrixXX& H, VectorX& h) const;

  /**
   * @brief findMaximumAcceleration Find the maximal acceleration along a given direction
          find          b, alpha0
          maximize      alpha0
          subject to    -h <= [-G  (Hv)] [b a0]^T   <= -h
                        0       <= [b a0]^T <= Inf


          b         are the coefficient of the contact force generators (f = V b)
          v         is the vector3 defining the direction of the motion
          alpha0    is the maximal amplitude of the acceleration, for the given direction v
          c         is the CoM position
          G         is the matrix whose columns are the gravito-inertial wrench generators
          A         is [-G  (Hv)]
   * @param A
   * @param h
   * @param alpha0
   * @return The status of the LP solver.
   */
  LP_status findMaximumAcceleration(Cref_matrixXX A, Cref_vector6 h, double& alpha0);

  /**
   * @brief checkAdmissibleAcceleration return true if the given acceleration is admissible for the given contacts
          find          b
          subject to    G b = Ha + h
                        0       <= b <= Inf
          b         are the coefficient of the contact force generators (f = V b)
          a         is the vector3 defining the acceleration
          G         is the matrix whose columns are the gravito-inertial wrench generators
          h and H come from polytope inequalities
   * @param G
   * @param H
   * @param h
   * @param a
   * @return true if the acceleration is admissible, false otherwise
   */
  bool checkAdmissibleAcceleration(Cref_matrixXX G, Cref_matrixXX H, Cref_vector6 h, Cref_vector3 a );

};

} // end namespace centroidal_dynamics

#endif

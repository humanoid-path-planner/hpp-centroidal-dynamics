#ifndef ROBUST_EQUILIBRIUM_LIB_STATIC_EQUILIBRIUM_H
#define ROBUST_EQUILIBRIUM_LIB_STATIC_EQUILIBRIUM_H

#include <Eigen/Dense>
#include <robust-equilibrium-lib/config.hh>
#include <robust-equilibrium-lib/solver_LP_abstract.hh>

namespace robust_equilibrium
{

enum ROBUST_EQUILIBRIUM_DLLAPI StaticEquilibriumAlgorithm
{
  STATIC_EQUILIBRIUM_ALGORITHM_LP,
  STATIC_EQUILIBRIUM_ALGORITHM_PP,
  STATIC_EQUILIBRIUM_ALGORITHM_IP,
  STATIC_EQUILIBRIUM_ALGORITHM_DIP
};

class ROBUST_EQUILIBRIUM_DLLAPI StaticEquilibrium
{
private:
  StaticEquilibriumAlgorithm  m_algorithm;
  SolverLPAbstract*           m_solver;

public:
  StaticEquilibrium(unsigned int generatorsPerContact, SolverLP solver_type);

  bool setNewContacts(Cref_matrixX3 contactPoints, Cref_matrixX3 contactNormals,
                      Cref_vectorX frictionCoefficients, StaticEquilibriumAlgorithm alg);

  double computeEquilibriumRobustness(Cref_vector2 com);

  double checkRobustEquilibrium(Cref_vector2 com, double e_max=0.0);

  double findExtremumOverLine(Cref_vector2 a, double b, double e_max=0.0);

  double findExtremumInDirection(Cref_vector2 direction, double e_max=0.0);

};

} // end namespace robust_equilibrium

#endif

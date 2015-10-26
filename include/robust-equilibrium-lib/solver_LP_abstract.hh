#ifndef ROBUST_EQUILIBRIUM_LIB_SOLVER_LP_ABSTRACT_HH
#define ROBUST_EQUILIBRIUM_LIB_SOLVER_LP_ABSTRACT_HH

#include <Eigen/Dense>
#include <robust-equilibrium-lib/config.hh>

namespace robust_equilibrium
{

enum ROBUST_EQUILIBRIUM_DLLAPI SolverLP
{
  SOLVER_LP_QPOASES
};

class ROBUST_EQUILIBRIUM_DLLAPI SolverLPAbstract
{
private:
  std::string m_name;

public:

  SolverLPAbstract(std::string name);

  /** Solve the linear program
   *  minimize    c' x
   *  subject to  lb <= A x <= ub
   */
  virtual bool solve(Cref_vectorX c, Cref_matrixXX A, Cref_vectorX lb, Cref_vectorX ub) = 0;

};

} // end namespace robust_equilibrium

#endif //ROBUST_EQUILIBRIUM_LIB_SOLVER_LP_ABSTRACT_HH

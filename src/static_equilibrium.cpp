#include <robust-equilibrium-lib/static_equilibrium.hh>

namespace robust_equilibrium
{


StaticEquilibrium::StaticEquilibrium(unsigned int generatorsPerContact, SolverLP solver_type)
{

}

bool StaticEquilibrium::setNewContacts(Cref_matrixX3 contactPoints, Cref_matrixX3 contactNormals,
                      Cref_vectorX frictionCoefficients, StaticEquilibriumAlgorithm alg)
{
  return true;
}

double StaticEquilibrium::computeEquilibriumRobustness(Cref_vector2 com)
{
  return 0.0;
}

double StaticEquilibrium::checkRobustEquilibrium(Cref_vector2 com, double e_max)
{
  return 0.0;
}

double StaticEquilibrium::findExtremumOverLine(Cref_vector2 a, double b, double e_max)
{
  return 0.0;
}

double StaticEquilibrium::findExtremumInDirection(Cref_vector2 direction, double e_max)
{
  return 0.0;
}


} // end namespace robust_equilibrium

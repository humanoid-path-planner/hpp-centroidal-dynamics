#include "centroidal-dynamics-lib/centroidal_dynamics.hh"

#include <boost/python.hpp>

namespace centroidal_dynamics
{
using namespace boost::python;

BOOST_PYTHON_MODULE(centroidal_dynamics)
{
    /** BEGIN enum types **/
  #ifdef CLP_FOUND
    enum_<SolverLP>("SolverLP")
            .value("SOLVER_LP_QPOASES", SOLVER_LP_QPOASES)
            .value("SOLVER_LP_CLP", SOLVER_LP_CLP)
            .export_values();
  #else
    enum_<SolverLP>("SolverLP")
            .value("SOLVER_LP_QPOASES", SOLVER_LP_QPOASES)
            .export_values();
  #endif


    enum_<EquilibriumAlgorithm>("EquilibriumAlgorithm")
            .value("EQUILIBRIUM_ALGORITHM_LP", EQUILIBRIUM_ALGORITHM_LP)
            .value("EQUILIBRIUM_ALGORITHM_LP2", EQUILIBRIUM_ALGORITHM_LP2)
            .value("EQUILIBRIUM_ALGORITHM_DLP", EQUILIBRIUM_ALGORITHM_DLP)
            .value("EQUILIBRIUM_ALGORITHM_PP", EQUILIBRIUM_ALGORITHM_PP)
            .value("EQUILIBRIUM_ALGORITHM_IP", EQUILIBRIUM_ALGORITHM_IP)
            .value("EQUILIBRIUM_ALGORITHM_DIP", EQUILIBRIUM_ALGORITHM_DIP)
            .export_values();

    /** END enum types **/
    class_<Equilibrium>("Equilibrium", init<std::string, double, unsigned int, optional <SolverLP, bool, const unsigned int, const bool> >())
            .def("useWarmStart", &Equilibrium::useWarmStart)
            .def("setUseWarmStart", &Equilibrium::setUseWarmStart)
            .def("getName", &Equilibrium::getName)
            .def("getAlgorithm", &Equilibrium::getAlgorithm)
    ;
}

} // namespace centroidal_dynamics

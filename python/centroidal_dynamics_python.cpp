#include "hpp/centroidal-dynamics/centroidal_dynamics.hh"
#include "hpp/centroidal-dynamics/util.hh"

#include <eigenpy/memory.hpp>
#include <eigenpy/eigenpy.hpp>

#include <boost/python.hpp>

EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(centroidal_dynamics::Equilibrium)

namespace centroidal_dynamics
{
using namespace boost::python;

boost::python::tuple wrapComputeQuasiEquilibriumRobustness(Equilibrium& self, const Vector3& com)
{
    double robustness;
    LP_status status = self.computeEquilibriumRobustness(com, robustness);
    return boost::python::make_tuple(status, robustness);
}

boost::python::tuple wrapComputeEquilibriumRobustness(Equilibrium& self, const Vector3& com, const Vector3& acc)
{
    double robustness;
    LP_status status = self.computeEquilibriumRobustness(com, acc, robustness);
    return boost::python::make_tuple(status, robustness);
}

boost::python::tuple wrapCheckRobustEquilibrium(Equilibrium& self, const Vector3& com)
{
    bool robustness;
    LP_status status = self.checkRobustEquilibrium(com, robustness);
    return boost::python::make_tuple(status, robustness);
}

bool wrapSetNewContacts(Equilibrium& self, const MatrixX3ColMajor& contactPoints, const MatrixX3ColMajor&  contactNormals,
                        const double frictionCoefficient, const EquilibriumAlgorithm alg)
{
    return self.setNewContacts(contactPoints, contactNormals, frictionCoefficient, alg);
}

bool wrapSetNewContactsFull(Equilibrium& self, const MatrixX3ColMajor& contactPoints, const MatrixX3ColMajor&  contactNormals,
                      const double frictionCoefficient, const EquilibriumAlgorithm alg)
{
    return self.setNewContacts(contactPoints, contactNormals, frictionCoefficient, alg);
}

boost::python::tuple wrapGetPolytopeInequalities(Equilibrium& self)
{
    MatrixXX H;
    VectorX h;
    H = MatrixXX::Zero(6,6);
    h = VectorX::Zero(6,1);
    self.getPolytopeInequalities(H,h);
    MatrixXXColMajor _H = H;
    return boost::python::make_tuple(_H, h);
}



BOOST_PYTHON_MODULE(hpp_centroidal_dynamics)
{
    /** BEGIN eigenpy init**/
    eigenpy::enableEigenPy();

    eigenpy::enableEigenPySpecific<MatrixX3ColMajor,MatrixX3ColMajor>();
    eigenpy::enableEigenPySpecific<MatrixXXColMajor,MatrixXXColMajor>();
    eigenpy::enableEigenPySpecific<Vector3,Vector3>();
    eigenpy::enableEigenPySpecific<VectorX,VectorX>();
    /*eigenpy::exposeAngleAxis();
    eigenpy::exposeQuaternion();*/

    /** END eigenpy init**/

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

    enum_<LP_status>("LP_status")
            .value("LP_STATUS_UNKNOWN", LP_STATUS_UNKNOWN)
            .value("LP_STATUS_OPTIMAL", LP_STATUS_OPTIMAL)
            .value("LP_STATUS_INFEASIBLE", LP_STATUS_INFEASIBLE)
            .value("LP_STATUS_UNBOUNDED", LP_STATUS_UNBOUNDED)
            .value("LP_STATUS_MAX_ITER_REACHED", LP_STATUS_MAX_ITER_REACHED)
            .value("LP_STATUS_ERROR", LP_STATUS_ERROR)
            .export_values();

    /** END enum types **/

    //bool (Equilibrium::*setNewContacts)
    //        (const MatrixX3ColMajor&, const MatrixX3ColMajor&, const double, const EquilibriumAlgorithm, const int graspIndex, const double maxGraspForce) = &Equilibrium::setNewContacts;

    class_<Equilibrium>("Equilibrium", init<std::string, double, unsigned int, optional <SolverLP, bool, const unsigned int, const bool> >())
            .def("useWarmStart", &Equilibrium::useWarmStart)
            .def("setUseWarmStart", &Equilibrium::setUseWarmStart)
            .def("getName", &Equilibrium::getName)
            .def("getAlgorithm", &Equilibrium::getAlgorithm)
            .def("setNewContacts", wrapSetNewContacts)
            .def("setNewContacts", wrapSetNewContactsFull)
            .def("computeEquilibriumRobustness", wrapComputeQuasiEquilibriumRobustness)
            .def("computeEquilibriumRobustness", wrapComputeEquilibriumRobustness)
            .def("checkRobustEquilibrium", wrapCheckRobustEquilibrium)
            .def("getPolytopeInequalities", wrapGetPolytopeInequalities)
    ;
}

} // namespace centroidal_dynamics

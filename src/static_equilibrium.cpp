/*
 * Copyright 2015, LAAS-CNRS
 * Author: Andrea Del Prete
 */

#include <robust-equilibrium-lib/static_equilibrium.hh>
#include <robust-equilibrium-lib/logger.hh>
#include <robust-equilibrium-lib/stop-watch.hh>
#include <iostream>
#include <vector>

using namespace std;

namespace robust_equilibrium
{

bool StaticEquilibrium::m_is_cdd_initialized = false;

StaticEquilibrium::StaticEquilibrium(string name, double mass, unsigned int generatorsPerContact, SolverLP solver_type)
{
  if(!m_is_cdd_initialized)
  {
    init_cdd_library();
    m_is_cdd_initialized = true;
  }

  if(generatorsPerContact!=4)
  {
    SEND_WARNING_MSG("Algorithm currently supports only 4 generators per contact!");
    generatorsPerContact = 4;
  }

  m_name = name;
  m_solver = Solver_LP_abstract::getNewSolver(solver_type);

  m_generatorsPerContact = generatorsPerContact;
  m_mass = mass;
  m_gravity.setZero();
  m_gravity(2) = -9.81;

  m_d.setZero();
  m_d.head<3>() = m_mass*m_gravity;
  m_D.setZero();
  m_D.block<3,2>(3,0) = crossMatrix(-m_mass*m_gravity).leftCols<2>();
}

bool StaticEquilibrium::setNewContacts(Cref_matrixX3 contactPoints, Cref_matrixX3 contactNormals,
                                       Cref_vectorX frictionCoefficients, StaticEquilibriumAlgorithm alg)
{
  assert(contactPoints.rows()==contactNormals.rows());
  assert(contactPoints.rows()==frictionCoefficients.rows());

  if(alg==STATIC_EQUILIBRIUM_ALGORITHM_IP)
  {
    SEND_ERROR_MSG("Algorithm IP not implemented yet");
    return false;
  }
  if(alg==STATIC_EQUILIBRIUM_ALGORITHM_DIP)
  {
    SEND_ERROR_MSG("Algorithm DIP not implemented yet");
    return false;
  }

  m_algorithm = alg;

  // compute tangent directions
  long int c = contactPoints.rows();
  int cg = m_generatorsPerContact;
  long int m = c*cg;           // number of generators
  double muu;
  m_T1.resize(c,3);
  m_T2.resize(c,3);
  m_A.resize(6,3*c);
  m_G.resize(3,m);
  m_G_centr.resize(6,m);

  for(long int i=0; i<c; i++)
  {
    // check that contact normals have norm 1
    if(fabs(contactNormals.row(i).norm()-1.0)>1e-6)
    {
      printf("ERROR Contact normals should have norm 1, this has norm %f", contactNormals.row(i).norm());
      return false;
    }
    // compute tangent directions
    m_T1.row(i) = contactNormals.row(i).cross(Vector3::UnitY());
    if(m_T1.row(i).norm()<1e-5)
      m_T1.row(i) = contactNormals.row(i).cross(Vector3::UnitX());
    m_T2.row(i) = contactNormals.row(i).transpose().cross(m_T1.row(i));
    m_T1.row(i).normalize();
    m_T2.row(i).normalize();

//    cout<<"Contact point "<<i<<"\nT1="<<m_T1.row(i)<<"\nT2="<<m_T2.row(i)<<"\n";

    // compute matrix mapping contact forces to gravito-inertial wrench
    m_A.block<3,3>(0, 3*i) = -Matrix3::Identity();
    m_A.block<3,3>(3, 3*i) = crossMatrix(-1.0*contactPoints.row(i).transpose());
//    cout<<"A:\n"<<m_A.block<6,3>(0,3*i)<<"\n";

    // compute generators
    muu = frictionCoefficients(i)/sqrt(2.0);
    m_G.col(cg*i+0) =  muu*m_T1.row(i) + muu*m_T2.row(i) + contactNormals.row(i);
    m_G.col(cg*i+1) =  muu*m_T1.row(i) - muu*m_T2.row(i) + contactNormals.row(i);
    m_G.col(cg*i+2) = -muu*m_T1.row(i) + muu*m_T2.row(i) + contactNormals.row(i);
    m_G.col(cg*i+3) = -muu*m_T1.row(i) - muu*m_T2.row(i) + contactNormals.row(i);

    // normalize generators
    m_G.col(cg*i+0).normalize();
    m_G.col(cg*i+1).normalize();
    m_G.col(cg*i+2).normalize();
    m_G.col(cg*i+3).normalize();

//    cout<<"Generators contact "<<i<<"\n"<<m_G.middleCols<4>(cg*i).transpose()<<"\n";
  }

  // project generators in 6d centroidal space
  for(unsigned int i=0; i<c; i++)
    m_G_centr.block(0,cg*i,6,cg) = m_A.block<6,3>(0,3*i) * m_G.block(0,cg*i,3,cg);
//  cout<<"G_centr:\n"<<m_G_centr.transpose()<<"\n";

  if(m_algorithm==STATIC_EQUILIBRIUM_ALGORITHM_PP)
  {
    if(!computePolytopeProjection(m_G_centr))
      return false;
    m_HD = m_H * m_D;
    m_Hd = m_H * m_d;
  }

  return true;
}


double StaticEquilibrium::computeEquilibriumRobustness(Cref_vector2 com)
{
  const long m = m_G_centr.cols(); // number of gravito-inertial wrench generators

  if(m_algorithm==STATIC_EQUILIBRIUM_ALGORITHM_LP)
  {
    /* Compute the robustness measure of the equilibrium of a specified CoM position
     * by solving the following LP:
          find          b, b0
          minimize      -b0
          subject to    D c + d <= G b    <= D c + d
                        0       <= b - b0 <= Inf
        where
          b         are the coefficient of the contact force generators (f = V b)
          b0        is the robustness measure
          c         is the CoM position
          G         is the matrix whose columns are the gravito-inertial wrench generators
    */
    VectorX b_b0(m+1);
    VectorX c = VectorX::Zero(m+1);
    c(m) = -1.0;
    VectorX lb = -VectorX::Ones(m+1)*1e5;
    VectorX ub = VectorX::Ones(m+1)*1e10;
    VectorX Alb = VectorX::Zero(6+m);
    VectorX Aub = VectorX::Ones(6+m)*1e100;
    MatrixXX A = MatrixXX::Zero(6+m, m+1);
    Alb.head<6>() = m_D * com + m_d;
    Aub.head<6>() = Alb.head<6>();
    A.topLeftCorner(6,m)      = m_G_centr;
    A.bottomLeftCorner(m,m)   = MatrixXX::Identity(m,m);
    A.bottomRightCorner(m,1)  = -VectorX::Ones(m);

    LP_status lpStatus_primal = m_solver->solve(c, lb, ub, A, Alb, Aub, b_b0);
    if(lpStatus_primal==LP_STATUS_OPTIMAL)
      return -1.0*m_solver->getObjectiveValue();

    SEND_ERROR_MSG("Primal LP problem could not be solved: "+toString(lpStatus_primal));
    return -1.0;
  }

  if(m_algorithm==STATIC_EQUILIBRIUM_ALGORITHM_LP2)
  {
    /* Compute the robustness measure of the equilibrium of a specified CoM position
     * by solving the following LP:
            find          b, b0
            minimize      -b0
            subject to    D c + d <= G (b + 1*b0) <= D c + d
                          0       <= b            <= Inf
          where
            b         are the (delta) coefficient of the contact force generators (f = G (b+b0))
            b0        is the robustness measure
            c         is the CoM position
            G         is the matrix whose columns are the gravito-inertial wrench generators
      */
    VectorX b_b0(m+1);
    VectorX c = VectorX::Zero(m+1);
    c(m) = -1.0;
    VectorX lb = VectorX::Zero(m+1);
    lb(m) = -1e10;
    VectorX ub = VectorX::Ones(m+1)*1e10;
    MatrixXX A = MatrixXX::Zero(6, m+1);
    Vector6 Alb = m_D * com + m_d;
    Vector6 Aub = Alb;
    A.leftCols(m)  = m_G_centr;
    A.rightCols(1) = m_G_centr * VectorX::Ones(m);

    LP_status lpStatus_primal = m_solver->solve(c, lb, ub, A, Alb, Aub, b_b0);
    if(lpStatus_primal==LP_STATUS_OPTIMAL)
    {
      return -1.0*m_solver->getObjectiveValue();
    }
    SEND_ERROR_MSG("Primal LP problem could not be solved: "+toString(lpStatus_primal));
    return -1.0;
  }

  if(m_algorithm==STATIC_EQUILIBRIUM_ALGORITHM_DLP)
  {
    /*Compute the robustness measure of the equilibrium of a specified CoM position
      by solving the following dual LP:
        find          v
        minimize      (d+D*com)' v
        subject to    G' v >= 0
                      1' G' v = 1
      where
        -(d+D c)' v   is the robustness measure
        c             is the CoM position
        G             is the matrix whose columns are the gravito-inertial wrench generators
     */
    Vector6 v;
    Vector6 c = m_D*com + m_d;
    Vector6 lb = Vector6::Ones()*-1e100;
    Vector6 ub = Vector6::Ones()*1e100;
    VectorX Alb = VectorX::Zero(m+1);
    Alb(m) = 1.0;
    VectorX Aub = VectorX::Ones(m+1)*1e100;
    Aub(m) = 1.0;
    MatrixX6 A(m+1,6);
    A.topRows(m) = m_G_centr.transpose();
    A.bottomRows<1>() = (m_G_centr*VectorX::Ones(m)).transpose();

    LP_status lpStatus_dual = m_solver->solve(c, lb, ub, A, Alb, Aub, v);
    if(lpStatus_dual==LP_STATUS_OPTIMAL)
      return m_solver->getObjectiveValue();

    SEND_ERROR_MSG("Dual LP problem could not be solved: "+toString(lpStatus_dual));
    return -1.0;
  }

  SEND_ERROR_MSG("checkRobustEquilibrium is not implemented for the specified algorithm");
  return -1.0;
}

bool StaticEquilibrium::checkRobustEquilibrium(Cref_vector2 com, double e_max)
{
  if(e_max!=0.0)
  {
    SEND_ERROR_MSG("checkRobustEquilibrium with e_max!=0 not implemented yet");
    return false;
  }
  if(m_algorithm!=STATIC_EQUILIBRIUM_ALGORITHM_PP)
  {
    SEND_ERROR_MSG("checkRobustEquilibrium is only implemented for the PP algorithm");
    return false;
  }

  VectorX res = m_HD * com + m_Hd;
  for(long i=0; i<res.size(); i++)
    if(res(i)>0.0)
      return false;
  return true;
}

bool StaticEquilibrium::findExtremumOverLine(Cref_vector2 a, Cref_vector2 a0, double e_max, Ref_vector2 com)
{
  const long m = m_G_centr.cols(); // number of gravito-inertial wrench generators

  if(m_algorithm==STATIC_EQUILIBRIUM_ALGORITHM_LP)
  {
    /* Compute the extremum CoM position over the line a*p + a0 that is in robust equilibrium
     * by solving the following LP:
          find          b, p
          minimize      -p
          subject to    D (a p + a0) + d <= G (b + b0) <= D (a p + a0) + d
                        0       <= b <= Inf
        where
          b0+b      are the coefficient of the contact force generators (f = G (b0+b))
          b0        is the robustness measure
          p         is the line parameter
          G         is the matrix whose columns are the gravito-inertial wrench generators
    */
    VectorX b_p(m+1);
    VectorX c = VectorX::Zero(m+1);
    c(m) = -1.0;
    VectorX lb = -VectorX::Zero(m+1);
    lb(m) = -1e5;
    VectorX ub = VectorX::Ones(m+1)*1e10;
    Vector6 Alb = m_D*a0 + m_d - m_G_centr*VectorX::Ones(m)*e_max;
    Vector6 Aub = Alb;
    Matrix6X A = Matrix6X::Zero(6, m+1);
    A.leftCols(m)     = m_G_centr;
    A.rightCols(1)    = -m_D*a;

    LP_status lpStatus_primal = m_solver->solve(c, lb, ub, A, Alb, Aub, b_p);
    if(lpStatus_primal==LP_STATUS_OPTIMAL)
    {
      com = a0 + a*b_p(m);
      return true;
    }

    SEND_ERROR_MSG("Primal LP problem could not be solved: "+toString(lpStatus_primal));
    return false;
  }

  if(m_algorithm==STATIC_EQUILIBRIUM_ALGORITHM_DLP)
  {
    /* Compute the extremum CoM position over the line a*x + a0 that is in robust equilibrium
     * by solving the following dual LP:
          find          v
          minimize      (D a0 + d -G b0)' v
          subject to    0  <= G' v    <= Inf
                        -1 <= a' D' v <= -1
        where
          G         is the matrix whose columns are the gravito-inertial wrench generators
    */
    Vector6 v;
    Vector6 c = m_D*a0 + m_d - m_G_centr*VectorX::Ones(m)*e_max;
    Vector6 lb = Vector6::Ones()*-1e100;
    Vector6 ub = Vector6::Ones()*1e100;
    VectorX Alb = VectorX::Zero(m+1);
    Alb(m) = -1.0;
    VectorX Aub = VectorX::Ones(m+1)*1e100;
    Aub(m) = -1.0;
    MatrixX6 A(m+1,6);
    A.topRows(m) = m_G_centr.transpose();
    A.bottomRows<1>() = (m_D*a).transpose();

    LP_status lpStatus_dual = m_solver->solve(c, lb, ub, A, Alb, Aub, v);
    if(lpStatus_dual==LP_STATUS_OPTIMAL)
    {
      double p = m_solver->getObjectiveValue();
      com = a0 + a*p;
      return true;
    }

    SEND_ERROR_MSG("Dual LP problem could not be solved: "+toString(lpStatus_dual));
    return false;
  }

  SEND_ERROR_MSG("checkRobustEquilibrium is not implemented for the specified algorithm");
  return false;
}

double StaticEquilibrium::findExtremumInDirection(Cref_vector2 direction, double e_max)
{
  SEND_ERROR_MSG("findExtremumInDirection not implemented yet");
  return 0.0;
}

bool StaticEquilibrium::computePolytopeProjection(Cref_matrix6X v)
{
//  getProfiler().start("eigen_to_cdd");
  dd_MatrixPtr V = cone_span_eigen_to_cdd(v.transpose());
//  getProfiler().stop("eigen_to_cdd");

  dd_ErrorType error = dd_NoError;

//  getProfiler().start("dd_DDMatrix2Poly");
  dd_PolyhedraPtr H_= dd_DDMatrix2Poly(V, &error);
//  getProfiler().stop("dd_DDMatrix2Poly");

  if(error != dd_NoError)
  {
    SEND_ERROR_MSG("numerical instability in cddlib. ill formed polytope");
    return false;
  }

//  getProfiler().start("cdd to eigen");
  dd_MatrixPtr b_A = dd_CopyInequalities(H_);
  // get equalities and add them as complementary inequality constraints
  std::vector<long> eq_rows;
  for(long elem=1;elem<=(long)(b_A->linset[0]);++elem)
  {
    if (set_member(elem,b_A->linset))
      eq_rows.push_back(elem);
  }
  int rowsize = (int)b_A->rowsize;
  m_H.resize(rowsize + eq_rows.size(), (int)b_A->colsize-1);
  m_h.resize(rowsize + eq_rows.size());
  for(int i=0; i < rowsize; ++i)
  {
    m_h(i) = (value_type)(*(b_A->matrix[i][0]));
    for(int j=1; j < b_A->colsize; ++j)
      m_H(i, j-1) = -(value_type)(*(b_A->matrix[i][j]));
  }
  int i = 0;
  for(std::vector<long int>::const_iterator cit = eq_rows.begin(); cit != eq_rows.end(); ++cit, ++i)
  {
    m_h(rowsize + i) = -m_h((int)(*cit));
    m_H(rowsize + i) = -m_H((int)(*cit));
  }
//  getProfiler().stop("cdd to eigen");

  return true;
}

} // end namespace robust_equilibrium

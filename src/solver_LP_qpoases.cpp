/*
 * Copyright 2015, LAAS-CNRS
 * Author: Andrea Del Prete
 */

#include <centroidal-dynamics-lib/solver_LP_qpoases.hh>
#include <centroidal-dynamics-lib/logger.hh>

USING_NAMESPACE_QPOASES

namespace centroidal_dynamics
{

  Solver_LP_qpoases::Solver_LP_qpoases(): Solver_LP_abstract()
  {
//    m_options.initialStatusBounds = ST_INACTIVE;
//    m_options.setToReliable();
    m_options.setToDefault();
    m_options.printLevel          = PL_NONE; //PL_LOW
    m_options.enableRegularisation = BT_TRUE;
    m_options.enableEqualities = BT_TRUE;
  }

  LP_status Solver_LP_qpoases::solve(Cref_vectorX c, Cref_vectorX lb, Cref_vectorX ub,
                                     Cref_matrixXX A, Cref_vectorX Alb, Cref_vectorX Aub,
                                     Ref_vectorX sol)
  {
    int n = (int)c.size();  // number of variables
    int m = (int)A.rows();  // number of constraints
    assert(lb.size()==n);
    assert(ub.size()==n);
    assert(A.cols()==n);
    assert(Alb.size()==m);
    assert(Aub.size()==m);

    int iters = m_maxIter;
    double solutionTime = m_maxTime;
    if(n!=m_solver.getNV() || m!=m_solver.getNC())
    {
      m_solver = SQProblem(n, m, HST_ZERO);
      m_solver.setOptions(m_options);
//      m_solver.printOptions();
      m_init_succeeded = false;
      m_H = MatrixXX::Zero(n,n);
    }

    if(!m_useWarmStart || !m_init_succeeded)
    {
      m_status = m_solver.init(NULL, c.data(), A.data(), lb.data(), ub.data(),
                               Alb.data(), Aub.data(), iters, &solutionTime);
      if(m_status==SUCCESSFUL_RETURN)
        m_init_succeeded = true;
    }
    else
    {
      // this doesn't work if I pass NULL instead of m_H.data()
      m_status = m_solver.hotstart(m_H.data(), c.data(), A.data(), lb.data(), ub.data(),
                                  Alb.data(), Aub.data(), iters, &solutionTime);
      if(m_status!=SUCCESSFUL_RETURN)
        m_init_succeeded = false;
    }

    if(m_status==SUCCESSFUL_RETURN)
    {
      m_solver.getPrimalSolution(sol.data());
    }

    return getStatus();
  }

  LP_status Solver_LP_qpoases::getStatus()
  {
    int ss = getSimpleStatus(m_status);
    if(ss==0)
      return LP_STATUS_OPTIMAL;
    if(ss==1)
      return LP_STATUS_MAX_ITER_REACHED;
    if(ss==-2)
       return LP_STATUS_INFEASIBLE;
    if(ss==-3)
      return LP_STATUS_UNBOUNDED;
    return LP_STATUS_ERROR;
  }

} // end namespace centroidal_dynamics

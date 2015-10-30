/*
 * Copyright 2015, LAAS-CNRS
 * Author: Andrea Del Prete
 */

#include <vector>
#include <iostream>
#include <robust-equilibrium-lib/static_equilibrium.hh>
#include <robust-equilibrium-lib/logger.hh>
#include <robust-equilibrium-lib/stop-watch.hh>

using namespace robust_equilibrium;
using namespace Eigen;
using namespace std;

#define PERF_PP "PP"
#define PERF_LP_PREPARATION "LP preparation"
#define PERF_LP_COIN "LP coin"
#define PERF_LP_OASES "LP oases"
#define PERF_LP2_COIN "LP2 coin"
#define PERF_LP2_OASES "LP2 oases"
#define PERF_DLP_COIN "DLP coin"
#define PERF_DLP_OASES "DLP oases"

int main()
{
//  init_cdd_library();
  srand ((unsigned int)(time(NULL)));

  int ret = 0;
  double mass = 70.0;
  double mu = 0.3;  // friction coefficient
  unsigned int generatorsPerContact = 4;
  unsigned int N_CONTACTS = 2;
  unsigned int N_POINTS = 50;   // number of com positions to test
  double MIN_FEET_DISTANCE = 0.3;
  double LX = 0.5*0.2172;        // half foot size in x direction
  double LY = 0.5*0.138;         // half foot size in y direction
  RVector3 CONTACT_POINT_LOWER_BOUNDS, CONTACT_POINT_UPPER_BOUNDS;
  CONTACT_POINT_LOWER_BOUNDS << 0.0,  0.0,  0.0;
  CONTACT_POINT_UPPER_BOUNDS << 0.5,  0.5,  0.5;
  double gamma = atan(mu);   // half friction cone angle
  RVector3 RPY_LOWER_BOUNDS, RPY_UPPER_BOUNDS;
  RPY_LOWER_BOUNDS << -0*gamma, -0*gamma, -M_PI;
  RPY_UPPER_BOUNDS << +0*gamma, +0*gamma, +M_PI;
  double X_MARG = 0.07;
  double Y_MARG = 0.07;
  const int GRID_SIZE = 15;

  StaticEquilibrium solver_PP(mass, generatorsPerContact, SOLVER_LP_CLP);
  StaticEquilibrium solver_LP_coin(mass, generatorsPerContact, SOLVER_LP_CLP);
  StaticEquilibrium solver_LP_oases(mass, generatorsPerContact, SOLVER_LP_QPOASES);
  StaticEquilibrium solver_LP2_coin(mass, generatorsPerContact, SOLVER_LP_CLP);
  StaticEquilibrium solver_LP2_oases(mass, generatorsPerContact, SOLVER_LP_QPOASES);
  StaticEquilibrium solver_DLP_coin(mass, generatorsPerContact, SOLVER_LP_CLP);
  StaticEquilibrium solver_DLP_oases(mass, generatorsPerContact, SOLVER_LP_QPOASES);
  MatrixXX contact_pos = MatrixXX::Zero(N_CONTACTS, 3);
  MatrixXX contact_rpy = MatrixXX::Zero(N_CONTACTS, 3);
  MatrixXX p = MatrixXX::Zero(4*N_CONTACTS,3); // contact points
  MatrixXX N = MatrixXX::Zero(4*N_CONTACTS,3); // contact normals
  VectorX frictionCoefficients(4*N_CONTACTS);
  frictionCoefficients.fill(mu);

  contact_pos << 0.122,  0.361,  0.071,
                 0.243,  0.029,  0.112;
  contact_rpy << 0.205, -0.005, -1.335,
                 -0.02 ,  0.206,  0.506;

  // Generate contact positions and orientations
  bool collision;
  for(unsigned int i=0; i<N_CONTACTS; i++)
  {
    while(true) // generate contact position
    {
      uniform(CONTACT_POINT_LOWER_BOUNDS, CONTACT_POINT_UPPER_BOUNDS, contact_pos.row(i));
      if(i==0)
        break;
      collision = false;
      for(unsigned int j=0; j<i-1; j++)
        if((contact_pos.row(i)-contact_pos.row(j)).norm() < MIN_FEET_DISTANCE)
          collision = true;
      if(collision==false)
        break;
    }

//     generate contact orientation
    uniform(RPY_LOWER_BOUNDS, RPY_UPPER_BOUNDS, contact_rpy.row(i));

    generate_rectangle_contacts(LX, LY, contact_pos.row(i), contact_rpy.row(i),
                                p.middleRows<4>(i*4), N.middleRows<4>(i*4));

    printf("Contact surface %d position (%.3f,%.3f,%.3f) ", i, contact_pos(i,0), contact_pos(i,1), contact_pos(i,2));
    printf("Orientation (%.3f,%.3f,%.3f)\n", contact_rpy(i,0), contact_rpy(i,1), contact_rpy(i,2));
  }

  for(int i=0; i<p.rows(); i++)
  {
    printf("Contact point %d position (%.3f,%.3f,%.3f) ", i, p(i,0), p(i,1), p(i,2));
    printf("Normal (%.3f,%.3f,%.3f)\n", N(i,0), N(i,1), N(i,2));
  }

  getProfiler().start(PERF_LP_PREPARATION);
  if(!solver_LP_coin.setNewContacts(p, N, frictionCoefficients, STATIC_EQUILIBRIUM_ALGORITHM_LP))
  {
    printf("Error while setting new contacts");
    return -1;
  }
  getProfiler().stop(PERF_LP_PREPARATION);
  getProfiler().start(PERF_LP_PREPARATION);
  if(!solver_LP_oases.setNewContacts(p, N, frictionCoefficients, STATIC_EQUILIBRIUM_ALGORITHM_LP))
  {
    printf("Error while setting new contacts");
    return -1;
  }
  getProfiler().stop(PERF_LP_PREPARATION);
  getProfiler().start(PERF_LP_PREPARATION);
  if(!solver_LP2_coin.setNewContacts(p, N, frictionCoefficients, STATIC_EQUILIBRIUM_ALGORITHM_LP2))
  {
    printf("Error while setting new contacts");
    return -1;
  }
  getProfiler().stop(PERF_LP_PREPARATION);
  getProfiler().start(PERF_LP_PREPARATION);
  if(!solver_LP2_oases.setNewContacts(p, N, frictionCoefficients, STATIC_EQUILIBRIUM_ALGORITHM_LP2))
  {
    printf("Error while setting new contacts");
    return -1;
  }
  getProfiler().stop(PERF_LP_PREPARATION);
  getProfiler().start(PERF_LP_PREPARATION);
  if(!solver_DLP_coin.setNewContacts(p, N, frictionCoefficients, STATIC_EQUILIBRIUM_ALGORITHM_DLP))
  {
    printf("Error while setting new contacts");
    return -1;
  }
  getProfiler().stop(PERF_LP_PREPARATION);
  getProfiler().start(PERF_LP_PREPARATION);
  if(!solver_DLP_oases.setNewContacts(p, N, frictionCoefficients, STATIC_EQUILIBRIUM_ALGORITHM_DLP))
  {
    printf("Error while setting new contacts");
    return -1;
  }
  getProfiler().stop(PERF_LP_PREPARATION);

  getProfiler().start(PERF_PP);
  bool res = solver_PP.setNewContacts(p, N, frictionCoefficients, STATIC_EQUILIBRIUM_ALGORITHM_PP);
  getProfiler().stop(PERF_PP);
  if(!res)
  {
    printf("Error while setting new contacts");
    return -1;
  }

  RVector2 com_LB, com_UB;
  com_LB(0) = p.col(0).minCoeff()-X_MARG;
  com_UB(0) = p.col(0).maxCoeff()+X_MARG;
  com_LB(1) = p.col(1).minCoeff()-Y_MARG;
  com_UB(1) = p.col(1).maxCoeff()+Y_MARG;

  MatrixXi contactPointCoord(4*N_CONTACTS,2);
  VectorX minDistContactPoint = 1e10*VectorX::Ones(4*N_CONTACTS);

  VectorX x_range(GRID_SIZE), y_range(GRID_SIZE);
  x_range.setLinSpaced(GRID_SIZE,com_LB(0),com_UB(0));
  y_range.setLinSpaced(GRID_SIZE,com_LB(1),com_UB(1));
  MatrixXX equilibrium(GRID_SIZE, GRID_SIZE);
  Vector2 com;
  cout<<"Gonna test equilibrium on a 2d grid of "<<GRID_SIZE<<"X"<<GRID_SIZE<<" points ";
  cout<<"ranging from "<<com_LB<<" to "<<com_UB<<endl;
  for(unsigned int i=0; i<GRID_SIZE; i++)
  {
    com(1) = y_range(GRID_SIZE-1-i);
    for(unsigned int j=0; j<GRID_SIZE; j++)
    {
//      uniform(com_LB, com_UB, com);
      com(0) = x_range(j);

      getProfiler().start(PERF_LP_COIN);
      double rob_coin  = solver_LP_coin.computeEquilibriumRobustness(com);
      getProfiler().stop(PERF_LP_COIN);

      getProfiler().start(PERF_LP_OASES);
      double rob_oases = solver_LP_oases.computeEquilibriumRobustness(com);
      getProfiler().stop(PERF_LP_OASES);

      getProfiler().start(PERF_LP2_COIN);
      double rob_coin2  = solver_LP2_coin.computeEquilibriumRobustness(com);
      getProfiler().stop(PERF_LP2_COIN);

      getProfiler().start(PERF_LP2_OASES);
      double rob_oases2 = solver_LP2_oases.computeEquilibriumRobustness(com);
      getProfiler().stop(PERF_LP2_OASES);

      getProfiler().start(PERF_DLP_COIN);
      double rob_DLP_coin  = solver_DLP_coin.computeEquilibriumRobustness(com);
      getProfiler().stop(PERF_DLP_COIN);

      getProfiler().start(PERF_DLP_OASES);
      double rob_DLP_oases = solver_DLP_oases.computeEquilibriumRobustness(com);
      getProfiler().stop(PERF_DLP_OASES);

      if(fabs(rob_coin-rob_oases)>1e-5)
        SEND_ERROR_MSG("LP_coin and LP_oases returned different results: "+toString(rob_coin)+" VS "+toString(rob_oases));

      if(fabs(rob_coin-rob_oases2)>1e-5)
        SEND_ERROR_MSG("LP_coin and LP2_oases returned different results: "+toString(rob_coin)+" VS "+toString(rob_oases2));

//      if(fabs(rob_coin-rob_coin2)>1e-5)
//        SEND_ERROR_MSG("LP_coin and LP2_coin returned different results: "+toString(rob_coin)+" VS "+toString(rob_coin2));

      if(fabs(rob_coin-rob_DLP_oases)>1e-5)
        SEND_ERROR_MSG("LP_coin and DLP_oases returned different results: "+toString(rob_coin)+" VS "+toString(rob_DLP_oases));

      if(fabs(rob_coin-rob_DLP_coin)>1e-5)
        SEND_ERROR_MSG("LP_coin and DLP_coin returned different results: "+toString(rob_coin)+" VS "+toString(rob_DLP_coin));

      if(solver_PP.checkRobustEquilibrium(com, 0.0))
      {
        equilibrium(i,j) = 1.0;
        if(rob_coin<=0)
          SEND_ERROR_MSG("PP says com is in equilibrium but LP_coin find negative robustness "+toString(rob_coin));
        if(rob_oases<=0)
          SEND_ERROR_MSG("PP says com is in equilibrium but LP_oases find negative robustness "+toString(rob_oases));
        if(rob_coin2<=0)
          SEND_ERROR_MSG("PP says com is in equilibrium but LP_coin2 find negative robustness "+toString(rob_coin2));
        if(rob_oases2<=0)
          SEND_ERROR_MSG("PP says com is in equilibrium but LP_oases2 find negative robustness "+toString(rob_oases2));
        if(rob_DLP_oases<=0)
          SEND_ERROR_MSG("PP says com is in equilibrium but DLP_oases find negative robustness "+toString(rob_DLP_oases));
        if(rob_DLP_coin<=0)
          SEND_ERROR_MSG("PP says com is in equilibrium but DLP_coin negative robustness "+toString(rob_DLP_coin));

        if(rob_coin>9.0)
          rob_coin = 9.0;
        printf("%d ", (int)rob_coin);
      }
      else
      {
        equilibrium(i,j) = 0.0;
        if(rob_coin>0)
          SEND_ERROR_MSG("PP says com is NOT in equilibrium but LP_coin find positive robustness "+toString(rob_coin));
        if(rob_oases>0)
          SEND_ERROR_MSG("PP says com is NOT in equilibrium but LP_oases find positive robustness "+toString(rob_oases));
        if(rob_coin2>0)
          SEND_ERROR_MSG("PP says com is NOT in equilibrium but LP_coin2 find positive robustness "+toString(rob_coin2));
        if(rob_oases2>0)
          SEND_ERROR_MSG("PP says com is NOT in equilibrium but LP_oases2 find positive robustness "+toString(rob_oases2));
        if(rob_DLP_coin>0)
          SEND_ERROR_MSG("PP says com is NOT in equilibrium but DLP_coin find positive robustness "+toString(rob_DLP_coin));
        if(rob_DLP_oases>0)
          SEND_ERROR_MSG("PP says com is NOT in equilibrium but DLP_oases find positive robustness "+toString(rob_DLP_oases));
        printf("- ");
      }

      // look for contact point positions on grid
      for(long k=0; k<4*N_CONTACTS; k++)
      {
        double dist = (p.block<1,2>(k,0) - com.transpose()).norm();
        if(dist < minDistContactPoint(k))
        {
          minDistContactPoint(k) = dist;
          contactPointCoord(k,0) = i;
          contactPointCoord(k,1) = j;
        }
      }
    }
    printf("\n");
  }

//  cout<<"Max dist between contact points and grid points "<<minDistContactPoint.maxCoeff()<<"\n";

  cout<<"\nContact point position on the same grid\n";
  bool contactPointDrawn;
  for(unsigned int i=0; i<GRID_SIZE; i++)
  {
    for(unsigned int j=0; j<GRID_SIZE; j++)
    {
      contactPointDrawn = false;
      for(long k=0; k<4*N_CONTACTS; k++)
      {
        if(contactPointCoord(k,0)==i && contactPointCoord(k,1)==j)
        {
          cout<<"X ";
          contactPointDrawn = true;
          break;
        }
      }
      if(contactPointDrawn==false)
        cout<<"- ";
    }
    printf("\n");
  }

  getProfiler().report_all();

  return ret;
}

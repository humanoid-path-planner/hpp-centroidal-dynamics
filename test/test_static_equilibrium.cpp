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

#define PERF_PP "Polytope Projection"
#define PERF_LP_PREPARATION "Computation of GIWC generators"
#define PERF_LP_COIN "Compute Equilibrium Robustness with LP coin"
#define PERF_LP_OASES "Compute Equilibrium Robustness with LP oases"
#define PERF_LP2_COIN "Compute Equilibrium Robustness with LP2 coin"
#define PERF_LP2_OASES "Compute Equilibrium Robustness with LP2 oases"
#define PERF_DLP_COIN "Compute Equilibrium Robustness with DLP coin"
#define PERF_DLP_OASES "Compute Equilibrium Robustness with DLP oases"

#define EPS 1e-5  // required precision

/** Test two different solvers on the method StaticEquilibrium::computeEquilibriumRobustness.
 */
int test_computeEquilibriumRobustness(StaticEquilibrium solver_1, StaticEquilibrium solver_2, Cref_matrixXX comPositions,
                                      const char* PERF_STRING_1, const char* PERF_STRING_2, int verb=0)
{
  int error_counter = 0;
  for(unsigned int i=0; i<comPositions.rows(); i++)
  {
    getProfiler().start(PERF_STRING_1);
    double rob_1  = solver_1.computeEquilibriumRobustness(comPositions.row(i));
    getProfiler().stop(PERF_STRING_1);

    getProfiler().start(PERF_STRING_2);
    double rob_2 = solver_2.computeEquilibriumRobustness(comPositions.row(i));
    getProfiler().stop(PERF_STRING_2);

    if(fabs(rob_1-rob_2)>EPS)
    {
      if(verb>0)
        SEND_ERROR_MSG(solver_1.getName()+" and "+solver_2.getName()+" returned different results: "+toString(rob_1)+" VS "+toString(rob_2));
      error_counter++;
    }
  }

  SEND_INFO_MSG("Test computeEquilibriumRobustness "+solver_1.getName()+" VS "+solver_2.getName()+": "+toString(error_counter)+" error(s).");
  return error_counter;
}

/** Test method StaticEquilibrium::findExtremumOverLine. The test works in this way: first it
 * calls the method findExtremumOverLine of the solver to test to find the extremum over a random
 * line with a specified robustness. Then it checks that the point found really has the specified
 * robustness by using the ground-truth solver.
 */
int test_findExtremumOverLine(StaticEquilibrium solver_to_test, StaticEquilibrium solver_ground_truth,
                              Cref_vector2 a0, int N_TESTS, double e_max,
                              const char* PERF_STRING_TEST, const char* PERF_STRING_GROUND_TRUTH, int verb=0)
{
  int error_counter = 0;
  Vector2 a, com;
  bool status;
  double desired_robustness;
  for(unsigned int i=0; i<N_TESTS; i++)
  {
    uniform(-1.0*Vector2::Ones(), Vector2::Ones(), a);
    desired_robustness = (rand()/ value_type(RAND_MAX))*e_max;

    getProfiler().start(PERF_STRING_TEST);
    status  = solver_to_test.findExtremumOverLine(a, a0, desired_robustness, com);
    getProfiler().stop(PERF_STRING_TEST);

    if(status==false)
    {
      error_counter++;
      if(verb>0)
        SEND_ERROR_MSG(solver_to_test.getName()+" failed to find extremum over line");
      continue;
    }

    getProfiler().start(PERF_STRING_GROUND_TRUTH);
    double robustness = solver_ground_truth.computeEquilibriumRobustness(com);
    getProfiler().stop(PERF_STRING_GROUND_TRUTH);

    if(fabs(robustness-desired_robustness)>EPS)
    {
      if(verb>0)
        SEND_ERROR_MSG(solver_to_test.getName()+" found this com position: "+toString(com)+
                       " which should have robustness "+toString(desired_robustness)+
                       " but actually has robustness "+toString(robustness));
      error_counter++;
    }
  }

  SEND_INFO_MSG("Test findExtremumOverLine "+solver_to_test.getName()+" VS "+solver_ground_truth.getName()+": "+toString(error_counter)+" error(s).");
  return error_counter;
}

/** Draw a grid on the screen using the robustness computed with the method
 *  StaticEquilibrium::computeEquilibriumRobustness.
 */
void drawRobustnessGrid(StaticEquilibrium solver, Cref_matrixXX comPositions)
{
  int grid_size = (int)sqrt(comPositions.rows());
  for(unsigned int i=0; i<comPositions.rows(); i++)
  {
    double rob = solver.computeEquilibriumRobustness(comPositions.row(i));
    if(rob>=0.0)
    {
      if(rob>9.0)
        rob = 9.0;
      printf("%d ", (int)rob);
    }
    else
      printf("- ");
    if((i+1)%grid_size==0)
      printf("\n");
  }
}

int main()
{
  srand ((unsigned int)(time(NULL)));
  RVector3 CONTACT_POINT_LOWER_BOUNDS, CONTACT_POINT_UPPER_BOUNDS;
  RVector3 RPY_LOWER_BOUNDS, RPY_UPPER_BOUNDS;

  /************************************** USER PARAMETERS *******************************/
  double mass = 70.0;
  double mu = 0.3;  // friction coefficient
  unsigned int generatorsPerContact = 4;
  unsigned int N_CONTACTS = 2;
  double MIN_FEET_DISTANCE = 0.3;
  double LX = 0.5*0.2172;        // half foot size in x direction
  double LY = 0.5*0.138;         // half foot size in y direction
  CONTACT_POINT_LOWER_BOUNDS << 0.0,  0.0,  0.0;
  CONTACT_POINT_UPPER_BOUNDS << 0.5,  0.5,  0.5;
  double gamma = atan(mu);   // half friction cone angle
  RPY_LOWER_BOUNDS << -0*gamma, -0*gamma, -M_PI;
  RPY_UPPER_BOUNDS << +0*gamma, +0*gamma, +M_PI;
  double X_MARG = 0.07;
  double Y_MARG = 0.07;
  const int GRID_SIZE = 15;
  /************************************ END USER PARAMETERS *****************************/

  StaticEquilibrium solver_PP("PP", mass, generatorsPerContact, SOLVER_LP_CLP);
  StaticEquilibrium solver_LP_coin("LP coin", mass, generatorsPerContact, SOLVER_LP_CLP);
  StaticEquilibrium solver_LP_oases("LP oases", mass, generatorsPerContact, SOLVER_LP_QPOASES);
  StaticEquilibrium solver_LP2_coin("LP2 coin", mass, generatorsPerContact, SOLVER_LP_CLP);
  StaticEquilibrium solver_LP2_oases("LP2 oases", mass, generatorsPerContact, SOLVER_LP_QPOASES);
  StaticEquilibrium solver_DLP_coin("DLP coin", mass, generatorsPerContact, SOLVER_LP_CLP);
  StaticEquilibrium solver_DLP_oases("DLP oases", mass, generatorsPerContact, SOLVER_LP_QPOASES);

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
  MatrixXX comPositions(GRID_SIZE*GRID_SIZE, 2);
  cout<<"Gonna test equilibrium on a 2d grid of "<<GRID_SIZE<<"X"<<GRID_SIZE<<" points ";
  cout<<"ranging from "<<com_LB<<" to "<<com_UB<<endl;
  for(unsigned int i=0; i<GRID_SIZE; i++)
  {
    for(unsigned int j=0; j<GRID_SIZE; j++)
    {
      comPositions(i*GRID_SIZE+j, 1) = y_range(GRID_SIZE-1-i);
      comPositions(i*GRID_SIZE+j, 0) = x_range(j);

      // look for contact point positions on grid
      for(long k=0; k<4*N_CONTACTS; k++)
      {
        double dist = (p.block<1,2>(k,0) - comPositions.row(i*GRID_SIZE+j)).norm();
        if(dist < minDistContactPoint(k))
        {
          minDistContactPoint(k) = dist;
          contactPointCoord(k,0) = i;
          contactPointCoord(k,1) = j;
        }
      }
    }
  }

  cout<<"\nContact point positions\n";
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

  drawRobustnessGrid(solver_DLP_oases, comPositions);

  test_computeEquilibriumRobustness(solver_LP_coin, solver_LP_oases, comPositions, PERF_LP_COIN, PERF_LP_OASES);
  test_computeEquilibriumRobustness(solver_LP_coin, solver_LP2_coin, comPositions, PERF_LP_COIN, PERF_LP2_COIN);
  test_computeEquilibriumRobustness(solver_LP_coin, solver_LP2_oases, comPositions, PERF_LP_COIN, PERF_LP2_OASES);
  test_computeEquilibriumRobustness(solver_LP_coin, solver_DLP_coin, comPositions, PERF_LP_COIN, PERF_DLP_COIN);
  test_computeEquilibriumRobustness(solver_LP_coin, solver_DLP_oases, comPositions, PERF_LP_COIN, PERF_DLP_OASES);

  Vector2 a0 = 0.5*(com_LB+com_UB);
  const int N_TESTS = 100;
  const double E_MAX = 5.0;
  test_findExtremumOverLine(solver_LP_oases, solver_DLP_oases, a0, N_TESTS, E_MAX, "EXTREMUM OVER LINE LP OASES", PERF_DLP_OASES, 1);
  test_findExtremumOverLine(solver_DLP_oases, solver_DLP_oases, a0, N_TESTS, E_MAX, "EXTREMUM OVER LINE DLP OASES", PERF_DLP_OASES, 1);

  getProfiler().report_all();

  return 0;
}

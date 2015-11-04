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

#define EPS 1e-3  // required precision

/** Check the coherence between the method StaticEquilibrium::computeEquilibriumRobustness
 * and the method StaticEquilibrium::checkRobustEquilibrium.
 * @param solver_1 Solver used to test computeEquilibriumRobustness.
 * @param solver_2 Solver used to test checkRobustEquilibrium.
 * @param comPositions List of 2d com positions on which to perform the tests.
 * @param PERF_STRING_1 String to use for logging the computation times of solver_1
 * @param PERF_STRING_2 String to use for logging the computation times of solver_2
 * @param verb Verbosity level, 0 print nothing, 1 print summary, 2 print everything
 */
int test_computeEquilibriumRobustness_vs_checkEquilibrium(StaticEquilibrium solver_1,
                                                          StaticEquilibrium solver_2,
                                                          Cref_matrixXX comPositions,
                                                          const char* PERF_STRING_1=NULL,
                                                          const char* PERF_STRING_2=NULL,
                                                          int verb=0)
{
  int error_counter = 0;
  double rob;
  LP_status status;
  bool equilibrium;
  for(unsigned int i=0; i<comPositions.rows(); i++)
  {
    if(PERF_STRING_1!=NULL)
      getProfiler().start(PERF_STRING_1);
    status = solver_1.computeEquilibriumRobustness(comPositions.row(i), rob);
    if(PERF_STRING_1!=NULL)
      getProfiler().stop(PERF_STRING_1);

    if(status!=LP_STATUS_OPTIMAL)
    {
      if(verb>1)
        SEND_ERROR_MSG(solver_1.getName()+" failed to compute robustness of com position "+toString(comPositions.row(i)));
      error_counter++;
      continue;
    }

    if(PERF_STRING_2!=NULL)
      getProfiler().start(PERF_STRING_2);
    status= solver_2.checkRobustEquilibrium(comPositions.row(i), equilibrium);
    if(PERF_STRING_2!=NULL)
      getProfiler().stop(PERF_STRING_2);

    if(status!=LP_STATUS_OPTIMAL)
    {
      if(verb>1)
        SEND_ERROR_MSG(solver_2.getName()+" failed to check equilibrium of com position "+toString(comPositions.row(i)));
      error_counter++;
      continue;
    }

    if(equilibrium==true && rob<0.0)
    {
      if(verb>1)
        SEND_ERROR_MSG(solver_2.getName()+" says com is in equilibrium while "+solver_1.getName()+" computed a negative robustness measure "+toString(rob));
      error_counter++;
    }
    else if(equilibrium==false && rob>0.0)
    {
      if(verb>1)
        SEND_ERROR_MSG(solver_2.getName()+" says com is not in equilibrium while "+solver_1.getName()+" computed a positive robustness measure "+toString(rob));
      error_counter++;
    }
  }

  if(verb>0)
    cout<<"Test test_computeEquilibriumRobustness_vs_checkEquilibrium "+solver_1.getName()+" VS "+solver_2.getName()+": "+toString(error_counter)+" error(s).\n";
  return error_counter;
}

/** Test two different solvers on the method StaticEquilibrium::computeEquilibriumRobustness.
 * @param solver_1 First solver to test.
 * @param solver_2 Second solver to test.
 * @param comPositions List of 2d com positions on which to perform the tests.
 * @param PERF_STRING_1 String to use for logging the computation times of solver_1
 * @param PERF_STRING_2 String to use for logging the computation times of solver_2
 * @param verb Verbosity level, 0 print nothing, 1 print summary, 2 print everything
 */
int test_computeEquilibriumRobustness(StaticEquilibrium solver_1, StaticEquilibrium solver_2, Cref_matrixXX comPositions,
                                      const char* PERF_STRING_1, const char* PERF_STRING_2, int verb=0)
{
  int error_counter = 0;
  double rob_1, rob_2;
  LP_status status;
  for(unsigned int i=0; i<comPositions.rows(); i++)
  {
    getProfiler().start(PERF_STRING_1);
    status = solver_1.computeEquilibriumRobustness(comPositions.row(i), rob_1);
    getProfiler().stop(PERF_STRING_1);

    if(status!=LP_STATUS_OPTIMAL)
    {
      if(verb>1)
        SEND_ERROR_MSG(solver_1.getName()+" failed to compute robustness of com position "+toString(comPositions.row(i)));
      error_counter++;
      continue;
    }

    getProfiler().start(PERF_STRING_2);
    status = solver_2.computeEquilibriumRobustness(comPositions.row(i), rob_2);
    getProfiler().stop(PERF_STRING_2);

    if(status!=LP_STATUS_OPTIMAL)
    {
      if(verb>1)
        SEND_ERROR_MSG(solver_2.getName()+" failed to compute robustness of com position "+toString(comPositions.row(i)));
      error_counter++;
      continue;
    }

    if(fabs(rob_1-rob_2)>EPS)
    {
      if(verb>1)
        SEND_ERROR_MSG(solver_1.getName()+" and "+solver_2.getName()+" returned different results: "+toString(rob_1)+" VS "+toString(rob_2));
      error_counter++;
    }
  }

  if(verb>0)
    cout<<"Test computeEquilibriumRobustness "+solver_1.getName()+" VS "+solver_2.getName()+": "+toString(error_counter)+" error(s).\n";
  return error_counter;
}

/** Test method StaticEquilibrium::findExtremumOverLine. The test works in this way: first it
 * calls the method findExtremumOverLine of the solver to test to find the extremum over a random
 * line with a specified robustness. Then it checks that the point found really has the specified
 * robustness by using the ground-truth solver.
 * @param solver_to_test Solver to test.
 * @param solver_ground_truth Second solver to use as ground truth.
 * @param a0 A 2d com position that allows for static equilibrium.
 * @param N_TESTS Number of tests to perform.
 * @param e_max Maximum value for the desired robustness.
 * @param PERF_STRING_TEST String to use for logging the computation times of solver_to_test
 * @param PERF_STRING_GROUND_TRUTH String to use for logging the computation times of solver_ground_truth
 * @param verb Verbosity level, 0 print nothing, 1 print summary, 2 print everything
 */
int test_findExtremumOverLine(StaticEquilibrium &solver_to_test, StaticEquilibrium &solver_ground_truth,
                              Cref_vector3 a0, int N_TESTS, double e_max,
                              const char* PERF_STRING_TEST, const char* PERF_STRING_GROUND_TRUTH, int verb=0)
{
  int error_counter = 0;
  Vector3 a, com;
  LP_status status;
  double desired_robustness, robustness;
  for(unsigned int i=0; i<N_TESTS; i++)
  {
    uniform(-1.0*Vector3::Ones(), Vector3::Ones(), a);
    if(e_max>=0.0)
      desired_robustness = (rand()/ value_type(RAND_MAX))*e_max;
    else
      desired_robustness = e_max - EPS;

    getProfiler().start(PERF_STRING_TEST);
    status  = solver_to_test.findExtremumOverLine(a, a0, desired_robustness, com);
    getProfiler().stop(PERF_STRING_TEST);

    if(status!=LP_STATUS_OPTIMAL)
    {
      status = solver_ground_truth.computeEquilibriumRobustness(a0, robustness);
      if(status!=LP_STATUS_OPTIMAL)
      {
        error_counter++;
        if(verb>1)
          SEND_ERROR_MSG(solver_ground_truth.getName()+" failed to compute equilibrium robustness of com position "+toString(a0.transpose()));
      }
      else if(robustness>desired_robustness)
      {
        error_counter++;
        if(verb>1)
          SEND_ERROR_MSG(solver_to_test.getName()+" failed to find extremum over line starting from "+
                         toString(a0.transpose())+" with robustness "+toString(desired_robustness)+" while "+
                         solver_ground_truth.getName()+" says this position has robustness "+toString(robustness));
      }
      continue;
    }

    getProfiler().start(PERF_STRING_GROUND_TRUTH);
    status = solver_ground_truth.computeEquilibriumRobustness(com, robustness);
    getProfiler().stop(PERF_STRING_GROUND_TRUTH);

    if(status!=LP_STATUS_OPTIMAL)
    {
      error_counter++;
      if(verb>1)
        SEND_ERROR_MSG(solver_ground_truth.getName()+" failed to compute equilibrium robustness of com posiiton "+toString(com.transpose()));
    }
    else if(fabs(robustness-desired_robustness)>EPS)
    {
      if(verb>1)
        SEND_ERROR_MSG(solver_to_test.getName()+" found this extremum: "+toString(com.transpose())+
                       " over the line starting at "+toString(a0.transpose())+" in direction "+toString(a.transpose())+
                       " which should have robustness "+toString(desired_robustness)+
                       " but actually has robustness "+toString(robustness));
      error_counter++;
    }
  }

  if(verb>0)
    cout<<"Test findExtremumOverLine "+solver_to_test.getName()+" VS "+solver_ground_truth.getName()+": "+toString(error_counter)+" error(s).\n";
  return error_counter;
}

/** Draw a grid on the screen using the robustness computed with the method
 *  StaticEquilibrium::computeEquilibriumRobustness.
 * @param solver The solver to use for computing the equilibrium robustness.
 * @param comPositions Grid of CoM positions in the form of an Nx2 matrix.
 */
void drawRobustnessGrid(StaticEquilibrium &solver, Cref_matrixXX comPositions)
{
  int grid_size = (int)sqrt(comPositions.rows());
  double rob ;
  LP_status status;
  for(unsigned int i=0; i<comPositions.rows(); i++)
  {
    status = solver.computeEquilibriumRobustness(comPositions.row(i), rob);
    if(status!=LP_STATUS_OPTIMAL)
    {
      SEND_ERROR_MSG("Faild to compute equilibrium robustness of com position "+toString(comPositions.row(i))+", error code "+toString(status));
      rob = -1.0;
    }

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

void testWithLoadedData()
{
  cout<<"*** TEST WITH LOADED DATA ***\n";

  double mass = 55.8836;
  double mu = 0.5;  // friction coefficient
  unsigned int generatorsPerContact = 4;
  string file_path = "../test_data/";
  double expected_robustness = 17.1222;

  const int N_SOLVERS = 3;
  string solverNames[] = {"LP oases", "LP2 oases", "DLP oases"};
  StaticEquilibriumAlgorithm algorithms[] = {STATIC_EQUILIBRIUM_ALGORITHM_LP,
                                             STATIC_EQUILIBRIUM_ALGORITHM_LP2,
                                             STATIC_EQUILIBRIUM_ALGORITHM_DLP};

  MatrixXX contactPoints, contactNormals;
  Vector3 com;
  if(!readMatrixFromFile(file_path+"positions.dat", contactPoints))
  {
    SEND_ERROR_MSG("Impossible to read positions from file");
    return;
  }
  if(!readMatrixFromFile(file_path+"normals.dat", contactNormals))
  {
    SEND_ERROR_MSG("Impossible to read normals from file");
    return;
  }
  if(!readMatrixFromFile(file_path+"com.dat", com))
  {
    SEND_ERROR_MSG("Impossible to read com from file");
    return;
  }

  StaticEquilibrium* solvers[N_SOLVERS];
  double robustness[N_SOLVERS];
  for(int s=0; s<N_SOLVERS; s++)
  {
    solvers[s] = new StaticEquilibrium(solverNames[s], mass, generatorsPerContact, SOLVER_LP_QPOASES);
    if(!solvers[s]->setNewContacts(contactPoints, contactNormals, mu, algorithms[s]))
    {
      SEND_ERROR_MSG("Error while setting new contacts for solver "+solvers[s]->getName());
      continue;
    }
    LP_status status = solvers[s]->computeEquilibriumRobustness(com, robustness[s]);
    if(status==LP_STATUS_OPTIMAL)
    {
      if(fabs(expected_robustness-robustness[s])>EPS)
        cout<<"[ERROR] Solver "<<solvers[s]->getName()<<" computed robustness "<<robustness[s]<<" rather than "<<expected_robustness<<endl;
    }
    else
      SEND_ERROR_MSG("Solver "+solvers[s]->getName()+" failed to compute robustness, error code "+toString(status));
  }
  cout<<"*** END TEST WITH LOADED DATA ***\n\n";
}

int main()
{
  testWithLoadedData();

  cout<<"*** TEST WITH RANDOMLY GENERATED DATA ***\n";
  unsigned int seed = (unsigned int)(time(NULL));
//  seed = 1446555515;
  srand (seed);
  cout<<"Initialize random number generator with seed "<<seed<<" (in case you wanna repeat the same test later)\n";

  RVector3 CONTACT_POINT_LOWER_BOUNDS, CONTACT_POINT_UPPER_BOUNDS;
  RVector3 RPY_LOWER_BOUNDS, RPY_UPPER_BOUNDS;

  /************************************** USER PARAMETERS *******************************/
  double mass = 55.0;
  double mu = 0.3;  // friction coefficient
  unsigned int generatorsPerContact = 4;
  unsigned int N_CONTACTS = 2;
  double MIN_FEET_DISTANCE = 0.3;
  double LX = 0.5*0.2172;        // half contact surface size in x direction
  double LY = 0.5*0.138;         // half contact surface size in y direction
  CONTACT_POINT_LOWER_BOUNDS << 0.0,  0.0,  0.0;
  CONTACT_POINT_UPPER_BOUNDS << 0.5,  0.5,  0.5;
  double gamma = atan(mu);   // half friction cone angle
  RPY_LOWER_BOUNDS << -2*gamma, -2*gamma, -M_PI;
  RPY_UPPER_BOUNDS << +2*gamma, +2*gamma, +M_PI;
  double X_MARG = 0.07;
  double Y_MARG = 0.07;
  const int GRID_SIZE = 10;
  /************************************ END USER PARAMETERS *****************************/

  cout<<"Number of contacts: "<<N_CONTACTS<<endl;
  cout<<"Number of generators per contact: "<<generatorsPerContact<<endl;

  StaticEquilibrium solver_PP("PP", mass, generatorsPerContact, SOLVER_LP_QPOASES);
  StaticEquilibrium solver_LP_oases("LP oases", mass, generatorsPerContact, SOLVER_LP_QPOASES);
  StaticEquilibrium solver_LP2_oases("LP2 oases", mass, generatorsPerContact, SOLVER_LP_QPOASES);
  StaticEquilibrium solver_DLP_oases("DLP oases", mass, generatorsPerContact, SOLVER_LP_QPOASES);

#ifdef CLP_FOUND
  StaticEquilibrium solver_LP_coin("LP coin", mass, generatorsPerContact, SOLVER_LP_CLP);
  StaticEquilibrium solver_LP2_coin("LP2 coin", mass, generatorsPerContact, SOLVER_LP_CLP);
  StaticEquilibrium solver_DLP_coin("DLP coin", mass, generatorsPerContact, SOLVER_LP_CLP);
#endif

  MatrixXX contact_pos = MatrixXX::Zero(N_CONTACTS, 3);
  MatrixXX contact_rpy = MatrixXX::Zero(N_CONTACTS, 3);
  MatrixXX p = MatrixXX::Zero(4*N_CONTACTS,3); // contact points
  MatrixXX N = MatrixXX::Zero(4*N_CONTACTS,3); // contact normals

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

  // compute upper and lower bounds of com positions to test
  RVector2 com_LB, com_UB;
  com_LB(0) = p.col(0).minCoeff()-X_MARG;
  com_UB(0) = p.col(0).maxCoeff()+X_MARG;
  com_LB(1) = p.col(1).minCoeff()-Y_MARG;
  com_UB(1) = p.col(1).maxCoeff()+Y_MARG;

  MatrixXi contactPointCoord(4*N_CONTACTS,2);
  VectorX minDistContactPoint = 1e10*VectorX::Ones(4*N_CONTACTS);

  // create grid of com positions to test
  VectorX x_range(GRID_SIZE), y_range(GRID_SIZE);
  x_range.setLinSpaced(GRID_SIZE,com_LB(0),com_UB(0));
  y_range.setLinSpaced(GRID_SIZE,com_LB(1),com_UB(1));
  MatrixXX comPositions(GRID_SIZE*GRID_SIZE, 3);
  cout<<"Gonna test equilibrium on a 2d grid of "<<GRID_SIZE<<"X"<<GRID_SIZE<<" points ";
  cout<<"ranging from "<<com_LB<<" to "<<com_UB<<endl;
  for(unsigned int i=0; i<GRID_SIZE; i++)
  {
    for(unsigned int j=0; j<GRID_SIZE; j++)
    {
      comPositions(i*GRID_SIZE+j, 1) = y_range(GRID_SIZE-1-i);
      comPositions(i*GRID_SIZE+j, 0) = x_range(j);
      comPositions(i*GRID_SIZE+j, 2) = 0.0;

      // look for contact point positions on grid
      for(long k=0; k<4*N_CONTACTS; k++)
      {
        double dist = (p.block<1,2>(k,0) - comPositions.block<1,2>(i*GRID_SIZE+j,0)).norm();
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
  if(!solver_LP_oases.setNewContacts(p, N, mu, STATIC_EQUILIBRIUM_ALGORITHM_LP))
  {
    printf("Error while setting new contacts");
    return -1;
  }
  getProfiler().stop(PERF_LP_PREPARATION);

  getProfiler().start(PERF_LP_PREPARATION);
  if(!solver_LP2_oases.setNewContacts(p, N, mu, STATIC_EQUILIBRIUM_ALGORITHM_LP2))
  {
    printf("Error while setting new contacts");
    return -1;
  }
  getProfiler().stop(PERF_LP_PREPARATION);

  getProfiler().start(PERF_LP_PREPARATION);
  if(!solver_DLP_oases.setNewContacts(p, N, mu, STATIC_EQUILIBRIUM_ALGORITHM_DLP))
  {
    printf("Error while setting new contacts");
    return -1;
  }
  getProfiler().stop(PERF_LP_PREPARATION);

  getProfiler().start(PERF_PP);
  bool res = solver_PP.setNewContacts(p, N, mu, STATIC_EQUILIBRIUM_ALGORITHM_PP);
  getProfiler().stop(PERF_PP);
  if(!res)
  {
    printf("Error while setting new contacts");
    return -1;
  }

#ifdef CLP_FOUND
  getProfiler().start(PERF_LP_PREPARATION);
  if(!solver_LP_coin.setNewContacts(p, N, mu, STATIC_EQUILIBRIUM_ALGORITHM_LP))
  {
    printf("Error while setting new contacts");
    return -1;
  }
  getProfiler().stop(PERF_LP_PREPARATION);

  getProfiler().start(PERF_LP_PREPARATION);
  if(!solver_LP2_coin.setNewContacts(p, N, mu, STATIC_EQUILIBRIUM_ALGORITHM_LP2))
  {
    printf("Error while setting new contacts");
    return -1;
  }
  getProfiler().stop(PERF_LP_PREPARATION);

  getProfiler().start(PERF_LP_PREPARATION);
  if(!solver_DLP_coin.setNewContacts(p, N, mu, STATIC_EQUILIBRIUM_ALGORITHM_DLP))
  {
    printf("Error while setting new contacts");
    return -1;
  }
  getProfiler().stop(PERF_LP_PREPARATION);
#endif


  cout<<"\nRobustness grid computed with DLP oases\n";
  drawRobustnessGrid(solver_DLP_oases, comPositions);

  test_computeEquilibriumRobustness(solver_DLP_oases, solver_LP_oases, comPositions, PERF_DLP_OASES, PERF_LP_OASES, 1);
  test_computeEquilibriumRobustness(solver_DLP_oases, solver_LP2_oases, comPositions, PERF_DLP_OASES, PERF_LP2_OASES, 1);

  test_computeEquilibriumRobustness_vs_checkEquilibrium(solver_LP_oases, solver_PP, comPositions, PERF_LP_OASES, NULL, 1);
  test_computeEquilibriumRobustness_vs_checkEquilibrium(solver_LP2_oases, solver_PP, comPositions, PERF_LP2_OASES, NULL, 1);
  test_computeEquilibriumRobustness_vs_checkEquilibrium(solver_DLP_oases, solver_PP, comPositions, PERF_DLP_OASES, NULL, 1);

#ifdef CLP_FOUND
  test_computeEquilibriumRobustness(solver_DLP_oases, solver_LP2_coin, comPositions, PERF_DLP_OASES, PERF_LP2_COIN, 1);
  test_computeEquilibriumRobustness(solver_DLP_oases, solver_DLP_coin, comPositions, PERF_DLP_OASES, PERF_DLP_COIN, 1);
  test_computeEquilibriumRobustness(solver_DLP_oases, solver_LP_coin, comPositions, PERF_DLP_OASES, PERF_LP_COIN, 1);

  test_computeEquilibriumRobustness_vs_checkEquilibrium(solver_LP_coin, solver_PP, comPositions, PERF_LP_COIN, NULL, 1);
  test_computeEquilibriumRobustness_vs_checkEquilibrium(solver_LP2_coin, solver_PP, comPositions, PERF_LP2_COIN, NULL, 1);
  test_computeEquilibriumRobustness_vs_checkEquilibrium(solver_DLP_coin, solver_PP, comPositions, PERF_DLP_COIN, NULL, 1);
#endif


  const int N_TESTS = 1000;
  Vector3 a0 = Vector3::Zero();
  a0.head<2>() = 0.5*(com_LB+com_UB);
  double e_max;
  LP_status status = solver_LP_oases.computeEquilibriumRobustness(a0, e_max);
  if(status!=LP_STATUS_OPTIMAL)
  {
    SEND_ERROR_MSG(solver_LP_oases.getName()+" failed to compute robustness of com position "+toString(a0.transpose())+", error code: "+toString(status));
  }
  else
  {
    test_findExtremumOverLine(solver_LP_oases, solver_DLP_oases, a0, N_TESTS, e_max, "EXTREMUM OVER LINE LP OASES", PERF_DLP_OASES, 2);
    test_findExtremumOverLine(solver_DLP_oases, solver_DLP_oases, a0, N_TESTS, e_max, "EXTREMUM OVER LINE DLP OASES", PERF_DLP_OASES, 2);
#ifdef CLP_FOUND
    test_findExtremumOverLine(solver_LP_coin, solver_LP_coin, a0, N_TESTS, e_max, "EXTREMUM OVER LINE LP COIN", PERF_LP_COIN, 2);
    test_findExtremumOverLine(solver_DLP_coin, solver_LP_coin, a0, N_TESTS, e_max, "EXTREMUM OVER LINE DLP COIN", PERF_LP_COIN, 2);
#endif
  }

  getProfiler().report_all();

  cout<<"*** END TEST WITH RANDOMLY GENERATED DATA ***\n";

  return 0;
}

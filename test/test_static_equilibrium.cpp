/*
 * Copyright 2015, LAAS-CNRS
 * Author: Andrea Del Prete
 */

#include <vector>
#include <iostream>
#include <hpp/centroidal-dynamics/centroidal_dynamics.hh>
#include <hpp/centroidal-dynamics/logger.hh>
#include <hpp/centroidal-dynamics/stop-watch.hh>

using namespace centroidal_dynamics;
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

/** Check the coherence between the method Equilibrium::computeEquilibriumRobustness
 * and the method Equilibrium::checkRobustEquilibrium.
 * @param solver_1 Solver used to test computeEquilibriumRobustness.
 * @param solver_2 Solver used to test checkRobustEquilibrium.
 * @param comPositions List of 2d com positions on which to perform the tests.
 * @param PERF_STRING_1 String to use for logging the computation times of solver_1
 * @param PERF_STRING_2 String to use for logging the computation times of solver_2
 * @param verb Verbosity level, 0 print nothing, 1 print summary, 2 print everything
 */
int test_computeEquilibriumRobustness_vs_checkEquilibrium(Equilibrium* solver_1, Equilibrium* solver_2,
                                                          Cref_matrixXX comPositions, const string& PERF_STRING_1 = "",
                                                          const string& PERF_STRING_2 = "", int verb = 0) {
  int error_counter = 0;
  double rob;
  LP_status status;
  bool equilibrium;
  for (unsigned int i = 0; i < comPositions.rows(); i++) {
    if (!PERF_STRING_1.empty()) getProfiler().start(PERF_STRING_1);
    status = solver_1->computeEquilibriumRobustness(comPositions.row(i).transpose(), rob);
    if (!PERF_STRING_1.empty()) getProfiler().stop(PERF_STRING_1);

    if (status != LP_STATUS_OPTIMAL) {
      if (verb > 1)
        SEND_ERROR_MSG(solver_1->getName() + " failed to compute robustness of com position " +
                       toString(comPositions.row(i)));
      error_counter++;
      continue;
    }

    if (!PERF_STRING_2.empty()) getProfiler().start(PERF_STRING_2);
    status = solver_2->checkRobustEquilibrium(comPositions.row(i).transpose(), equilibrium);
    if (!PERF_STRING_2.empty()) getProfiler().stop(PERF_STRING_2);

    if (status != LP_STATUS_OPTIMAL) {
      if (verb > 1)
        SEND_ERROR_MSG(solver_2->getName() + " failed to check equilibrium of com position " +
                       toString(comPositions.row(i)));
      error_counter++;
      continue;
    }

    if (equilibrium == true && rob < 0.0) {
      if (verb > 1)
        SEND_ERROR_MSG(solver_2->getName() + " says com is in equilibrium while " + solver_1->getName() +
                       " computed a negative robustness measure " + toString(rob));
      error_counter++;
    } else if (equilibrium == false && rob > 0.0) {
      if (verb > 1)
        SEND_ERROR_MSG(solver_2->getName() + " says com is not in equilibrium while " + solver_1->getName() +
                       " computed a positive robustness measure " + toString(rob));
      error_counter++;
    }
  }

  if (verb > 0)
    cout << "Test test_computeEquilibriumRobustness_vs_checkEquilibrium " + solver_1->getName() + " VS " +
                solver_2->getName() + ": " + toString(error_counter) + " error(s).\n";
  return error_counter;
}

/** Test two different solvers on the method Equilibrium::computeEquilibriumRobustness.
 * @param solver_1 First solver to test.
 * @param solver_2 Second solver to test.
 * @param comPositions List of 2d com positions on which to perform the tests.
 * @param PERF_STRING_1 String to use for logging the computation times of solver_1
 * @param PERF_STRING_2 String to use for logging the computation times of solver_2
 * @param verb Verbosity level, 0 print nothing, 1 print summary, 2 print everything
 */
int test_computeEquilibriumRobustness(Equilibrium* solver_1, Equilibrium* solver_2, Cref_matrixXX comPositions,
                                      const string& PERF_STRING_1, const string& PERF_STRING_2, int verb = 0) {
  int error_counter = 0;
  double rob_1, rob_2;
  LP_status status;
  for (unsigned int i = 0; i < comPositions.rows(); i++) {
    getProfiler().start(PERF_STRING_1);
    status = solver_1->computeEquilibriumRobustness(comPositions.row(i).transpose(), rob_1);
    getProfiler().stop(PERF_STRING_1);

    if (status != LP_STATUS_OPTIMAL) {
      if (verb > 1)
        SEND_ERROR_MSG(solver_1->getName() + " failed to compute robustness of com position " +
                       toString(comPositions.row(i)));
      error_counter++;
      continue;
    }

    getProfiler().start(PERF_STRING_2);
    status = solver_2->computeEquilibriumRobustness(comPositions.row(i).transpose(), rob_2);
    getProfiler().stop(PERF_STRING_2);

    if (status != LP_STATUS_OPTIMAL) {
      if (verb > 1)
        SEND_ERROR_MSG(solver_2->getName() + " failed to compute robustness of com position " +
                       toString(comPositions.row(i)));
      error_counter++;
      continue;
    }

    if (fabs(rob_1 - rob_2) > EPS) {
      if (verb > 1)
        SEND_ERROR_MSG(solver_1->getName() + " and " + solver_2->getName() +
                       " returned different results: " + toString(rob_1) + " VS " + toString(rob_2));
      error_counter++;
    }
  }

  if (verb > 0)
    cout << "Test computeEquilibriumRobustness " + solver_1->getName() + " VS " + solver_2->getName() + ": " +
                toString(error_counter) + " error(s).\n";
  return error_counter;
}

/** Test method Equilibrium::findExtremumOverLine. The test works in this way: first it
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
int test_findExtremumOverLine(Equilibrium* solver_to_test, Equilibrium* solver_ground_truth, Cref_vector3 a0,
                              int N_TESTS, double e_max, const string& PERF_STRING_TEST,
                              const string& PERF_STRING_GROUND_TRUTH, int verb = 0) {
  int error_counter = 0;
  centroidal_dynamics::Vector3 a, com;
  LP_status status;
  double desired_robustness, robustness;
  for (unsigned int i = 0; i < N_TESTS; i++) {
    uniform3(-1.0 * centroidal_dynamics::Vector3::Ones(), centroidal_dynamics::Vector3::Ones(), a);
    if (e_max >= 0.0)
      desired_robustness = (rand() / value_type(RAND_MAX)) * e_max;
    else
      desired_robustness = e_max - EPS;

    getProfiler().start(PERF_STRING_TEST);
    status = solver_to_test->findExtremumOverLine(a, a0, desired_robustness, com);
    getProfiler().stop(PERF_STRING_TEST);

    if (status != LP_STATUS_OPTIMAL) {
      status = solver_ground_truth->computeEquilibriumRobustness(a0, robustness);
      if (status != LP_STATUS_OPTIMAL) {
        error_counter++;
        if (verb > 1)
          SEND_ERROR_MSG(solver_ground_truth->getName() +
                         " failed to compute equilibrium robustness of com position " + toString(a0.transpose()));
      } else if (robustness > desired_robustness) {
        error_counter++;
        if (verb > 1)
          SEND_ERROR_MSG(solver_to_test->getName() + " failed to find extremum over line starting from " +
                         toString(a0.transpose()) + " with robustness " + toString(desired_robustness) + " while " +
                         solver_ground_truth->getName() + " says this position has robustness " +
                         toString(robustness));
      }
      continue;
    }

    getProfiler().start(PERF_STRING_GROUND_TRUTH);
    status = solver_ground_truth->computeEquilibriumRobustness(com, robustness);
    getProfiler().stop(PERF_STRING_GROUND_TRUTH);

    if (status != LP_STATUS_OPTIMAL) {
      error_counter++;
      if (verb > 1)
        SEND_ERROR_MSG(solver_ground_truth->getName() + " failed to compute equilibrium robustness of com posiiton " +
                       toString(com.transpose()));
    } else if (fabs(robustness - desired_robustness) > EPS) {
      if (verb > 1)
        SEND_ERROR_MSG(solver_to_test->getName() + " found this extremum: " + toString(com.transpose()) +
                       " over the line starting at " + toString(a0.transpose()) + " in direction " +
                       toString(a.transpose()) + " which should have robustness " + toString(desired_robustness) +
                       " but actually has robustness " + toString(robustness));
      error_counter++;
    }
  }

  if (verb > 0)
    cout << "Test findExtremumOverLine " + solver_to_test->getName() + " VS " + solver_ground_truth->getName() + ": " +
                toString(error_counter) + " error(s).\n";
  return error_counter;
}

/** Draw a grid on the screen using the robustness computed with the method
 *  Equilibrium::computeEquilibriumRobustness.
 * @param solver The solver to use for computing the equilibrium robustness.
 * @param comPositions Grid of CoM positions in the form of an Nx2 matrix.
 */
void drawRobustnessGrid(int N_CONTACTS, int GRID_SIZE, Equilibrium* solver, Cref_matrixXX comPositions,
                        Cref_matrixXX p) {
  MatrixXi contactPointCoord(4 * N_CONTACTS, 2);
  centroidal_dynamics::VectorX minDistContactPoint = 1e10 * centroidal_dynamics::VectorX::Ones(4 * N_CONTACTS);

  // create grid of com positions to test
  for (unsigned int i = 0; i < GRID_SIZE; i++) {
    for (unsigned int j = 0; j < GRID_SIZE; j++) {
      // look for contact point positions on grid
      for (long k = 0; k < 4 * N_CONTACTS; k++) {
        double dist = (p.block<1, 2>(k, 0) - comPositions.block<1, 2>(i * GRID_SIZE + j, 0)).norm();
        if (dist < minDistContactPoint(k)) {
          minDistContactPoint(k) = dist;
          contactPointCoord(k, 0) = i;
          contactPointCoord(k, 1) = j;
        }
      }
    }
  }

  cout << "\nContact point positions\n";
  bool contactPointDrawn;
  for (unsigned int i = 0; i < GRID_SIZE; i++) {
    for (unsigned int j = 0; j < GRID_SIZE; j++) {
      contactPointDrawn = false;
      for (long k = 0; k < 4 * N_CONTACTS; k++) {
        if (contactPointCoord(k, 0) == i && contactPointCoord(k, 1) == j) {
          cout << "X ";
          contactPointDrawn = true;
          break;
        }
      }
      if (contactPointDrawn == false) cout << "- ";
    }
    printf("\n");
  }

  cout << "\nRobustness grid computed with solver " << solver->getName() << endl;
  int grid_size = (int)sqrt(comPositions.rows());
  double rob;
  LP_status status;
  for (unsigned int i = 0; i < comPositions.rows(); i++) {
    status = solver->computeEquilibriumRobustness(comPositions.row(i).transpose(), rob);
    if (status != LP_STATUS_OPTIMAL) {
      SEND_ERROR_MSG("Faild to compute equilibrium robustness of com position " + toString(comPositions.row(i)) +
                     ", error code " + toString(status));
      rob = -1.0;
    }

    if (rob >= 0.0) {
      if (rob > 9.0) rob = 9.0;
      printf("%d ", (int)rob);
    } else
      printf("- ");
    if ((i + 1) % grid_size == 0) printf("\n");
  }
}

void generateContacts(unsigned int N_CONTACTS, double MIN_CONTACT_DISTANCE, double LX, double LY,
                      RVector3& CONTACT_POINT_LOWER_BOUNDS, RVector3& CONTACT_POINT_UPPER_BOUNDS,
                      RVector3& RPY_LOWER_BOUNDS, RVector3& RPY_UPPER_BOUNDS, centroidal_dynamics::MatrixX3& p,
                      centroidal_dynamics::MatrixX3& N) {
  MatrixXX contact_pos = MatrixXX::Zero(N_CONTACTS, 3);
  MatrixXX contact_rpy = MatrixXX::Zero(N_CONTACTS, 3);
  p.setZero(4 * N_CONTACTS, 3);  // contact points
  N.setZero(4 * N_CONTACTS, 3);  // contact normals

  // Generate contact positions and orientations
  bool collision;
  for (unsigned int i = 0; i < N_CONTACTS; i++) {
    while (true)  // generate contact position
    {
      uniform(CONTACT_POINT_LOWER_BOUNDS, CONTACT_POINT_UPPER_BOUNDS, contact_pos.row(i));
      if (i == 0) break;
      collision = false;
      for (unsigned int j = 0; j < i - 1; j++)
        if ((contact_pos.row(i) - contact_pos.row(j)).norm() < MIN_CONTACT_DISTANCE) collision = true;
      if (collision == false) break;
    }

    //     generate contact orientation
    uniform(RPY_LOWER_BOUNDS, RPY_UPPER_BOUNDS, contact_rpy.row(i));
    generate_rectangle_contacts(LX, LY, contact_pos.row(i).transpose(), contact_rpy.row(i).transpose(),
                                p.middleRows<4>(i * 4), N.middleRows<4>(i * 4));
    //    printf("Contact surface %d position (%.3f,%.3f,%.3f) ", i, contact_pos(i,0), contact_pos(i,1),
    //    contact_pos(i,2)); printf("Orientation (%.3f,%.3f,%.3f)\n", contact_rpy(i,0), contact_rpy(i,1),
    //    contact_rpy(i,2));
  }

  //  for(int i=0; i<p.rows(); i++)
  //  {
  //    printf("Contact point %d position (%.3f,%.3f,%.3f) ", i, p(i,0), p(i,1), p(i,2));
  //    printf("Normal (%.3f,%.3f,%.3f)\n", N(i,0), N(i,1), N(i,2));
  //  }
}

void testWithLoadedData() {
  cout << "*** TEST WITH LOADED DATA ***\n";

  double mass = 55.8836;
  double mu = 0.5;  // friction coefficient
  unsigned int generatorsPerContact = 4;
  string file_path = "../test_data/";
  double expected_robustness = 17.1222;

  const int N_SOLVERS = 3;
  string solverNames[] = {"LP oases", "LP2 oases", "DLP oases"};
  EquilibriumAlgorithm algorithms[] = {EQUILIBRIUM_ALGORITHM_LP, EQUILIBRIUM_ALGORITHM_LP2, EQUILIBRIUM_ALGORITHM_DLP};

  MatrixXX contactPoints, contactNormals;
  centroidal_dynamics::Vector3 com;
  if (!readMatrixFromFile(file_path + "positions.dat", contactPoints)) {
    SEND_ERROR_MSG("Impossible to read positions from file");
    return;
  }
  if (!readMatrixFromFile(file_path + "normals.dat", contactNormals)) {
    SEND_ERROR_MSG("Impossible to read normals from file");
    return;
  }
  if (!readMatrixFromFile(file_path + "com.dat", com)) {
    SEND_ERROR_MSG("Impossible to read com from file");
    return;
  }

  // this is a test to ensure that a matrixXX can be cast into a MatrixX3
  const centroidal_dynamics::MatrixX3& cp = contactPoints;
  const centroidal_dynamics::MatrixX3& cn = contactNormals;
  Equilibrium* solvers[N_SOLVERS];
  double robustness[N_SOLVERS];
  for (int s = 0; s < N_SOLVERS; s++) {
    solvers[s] = new Equilibrium(solverNames[s], mass, generatorsPerContact, SOLVER_LP_QPOASES);

    if (!solvers[s]->setNewContacts(cp, cn, mu, algorithms[s])) {
      SEND_ERROR_MSG("Error while setting new contacts for solver " + solvers[s]->getName());
      continue;
    }
    LP_status status = solvers[s]->computeEquilibriumRobustness(com, robustness[s]);
    if (status == LP_STATUS_OPTIMAL) {
      if (fabs(expected_robustness - robustness[s]) > EPS)
        cout << "[ERROR] Solver " << solvers[s]->getName() << " computed robustness " << robustness[s]
             << " rather than " << expected_robustness << endl;
    } else
      SEND_ERROR_MSG("Solver " + solvers[s]->getName() + " failed to compute robustness, error code " +
                     toString(status));
  }
  cout << "*** END TEST WITH LOADED DATA ***\n\n";
}

int main() {
  testWithLoadedData();

  cout << "*** TEST WITH RANDOMLY GENERATED DATA ***\n";
  unsigned int seed = (unsigned int)(time(NULL));
  //  seed = 1446555515;
  srand(seed);
  cout << "Initialize random number generator with seed " << seed
       << " (in case you wanna repeat the same test later)\n";

  RVector3 CONTACT_POINT_LOWER_BOUNDS, CONTACT_POINT_UPPER_BOUNDS;
  RVector3 RPY_LOWER_BOUNDS, RPY_UPPER_BOUNDS;

  /************************************** USER PARAMETERS *******************************/
  unsigned int N_TESTS = 10;
  double mass = 55.0;
  double mu = 0.3;  // friction coefficient
  unsigned int generatorsPerContact = 4;
  unsigned int N_CONTACTS = 2;
  double MIN_CONTACT_DISTANCE = 0.3;
  double LX = 0.5 * 0.2172;  // half contact surface size in x direction
  double LY = 0.5 * 0.138;   // half contact surface size in y direction
  CONTACT_POINT_LOWER_BOUNDS << 0.0, 0.0, 0.0;
  CONTACT_POINT_UPPER_BOUNDS << 0.5, 0.5, 0.5;
  double gamma = atan(mu);  // half friction cone angle
  RPY_LOWER_BOUNDS << -2 * gamma, -2 * gamma, -M_PI;
  RPY_UPPER_BOUNDS << +2 * gamma, +2 * gamma, +M_PI;
  double X_MARG = 0.07;
  double Y_MARG = 0.07;
  const int GRID_SIZE = 10;
  bool DRAW_CONTACT_POINTS = false;
  /************************************ END USER PARAMETERS *****************************/

#ifdef CLP_FOUND
  const int N_SOLVERS = 6;
  string solverNames[] = {"LP oases", "LP2 oases", "DLP oases", "LP coin", "LP2 coin", "DLP coin"};
  EquilibriumAlgorithm algorithms[] = {EQUILIBRIUM_ALGORITHM_LP, EQUILIBRIUM_ALGORITHM_LP2, EQUILIBRIUM_ALGORITHM_DLP,
                                       EQUILIBRIUM_ALGORITHM_LP, EQUILIBRIUM_ALGORITHM_LP2, EQUILIBRIUM_ALGORITHM_DLP};
  SolverLP lp_solver_types[] = {SOLVER_LP_QPOASES, SOLVER_LP_QPOASES, SOLVER_LP_QPOASES,
                                SOLVER_LP_CLP,     SOLVER_LP_CLP,     SOLVER_LP_CLP};
#else
  const int N_SOLVERS = 3;
  string solverNames[] = {"LP oases", "LP2 oases", "DLP oases"};
  EquilibriumAlgorithm algorithms[] = {EQUILIBRIUM_ALGORITHM_LP, EQUILIBRIUM_ALGORITHM_LP2, EQUILIBRIUM_ALGORITHM_DLP};
  SolverLP lp_solver_types[] = {SOLVER_LP_QPOASES, SOLVER_LP_QPOASES, SOLVER_LP_QPOASES};
#endif

  cout << "Number of contacts: " << N_CONTACTS << endl;
  cout << "Number of generators per contact: " << generatorsPerContact << endl;
  cout << "Gonna test equilibrium on a 2d grid of " << GRID_SIZE << "X" << GRID_SIZE << " points " << endl;

  Equilibrium* solver_PP = new Equilibrium("PP", mass, generatorsPerContact, SOLVER_LP_QPOASES);
  Equilibrium* solvers[N_SOLVERS];
  for (int s = 0; s < N_SOLVERS; s++)
    solvers[s] = new Equilibrium(solverNames[s], mass, generatorsPerContact, lp_solver_types[s]);

  centroidal_dynamics::MatrixX3 p, N;
  RVector2 com_LB, com_UB;
  centroidal_dynamics::VectorX x_range(GRID_SIZE), y_range(GRID_SIZE);
  MatrixXX comPositions(GRID_SIZE * GRID_SIZE, 3);
  for (unsigned n_test = 0; n_test < N_TESTS; n_test++) {
    generateContacts(N_CONTACTS, MIN_CONTACT_DISTANCE, LX, LY, CONTACT_POINT_LOWER_BOUNDS, CONTACT_POINT_UPPER_BOUNDS,
                     RPY_LOWER_BOUNDS, RPY_UPPER_BOUNDS, p, N);

    for (int s = 0; s < N_SOLVERS; s++) {
      getProfiler().start(PERF_LP_PREPARATION);
      if (!solvers[s]->setNewContacts(p, N, mu, algorithms[s])) {
        SEND_ERROR_MSG("Error while setting new contacts for solver " + solvers[s]->getName());
        return -1;
      }
      getProfiler().stop(PERF_LP_PREPARATION);
    }
    getProfiler().start(PERF_PP);
    if (!solver_PP->setNewContacts(p, N, mu, EQUILIBRIUM_ALGORITHM_PP)) {
      SEND_ERROR_MSG("Error while setting new contacts for solver " + solver_PP->getName());
      return -1;
    }
    getProfiler().stop(PERF_PP);

    // compute upper and lower bounds of com positions to test
    com_LB(0) = p.col(0).minCoeff() - X_MARG;
    com_UB(0) = p.col(0).maxCoeff() + X_MARG;
    com_LB(1) = p.col(1).minCoeff() - Y_MARG;
    com_UB(1) = p.col(1).maxCoeff() + Y_MARG;

    // create grid of com positions to test

    x_range.setLinSpaced(GRID_SIZE, com_LB(0), com_UB(0));
    y_range.setLinSpaced(GRID_SIZE, com_LB(1), com_UB(1));
    //    cout<<"ranging from "<<com_LB<<" to "<<com_UB<<endl;
    for (unsigned int i = 0; i < GRID_SIZE; i++) {
      for (unsigned int j = 0; j < GRID_SIZE; j++) {
        comPositions(i * GRID_SIZE + j, 1) = y_range(GRID_SIZE - 1 - i);
        comPositions(i * GRID_SIZE + j, 0) = x_range(j);
        comPositions(i * GRID_SIZE + j, 2) = 0.0;
      }
    }

    if (DRAW_CONTACT_POINTS) drawRobustnessGrid(N_CONTACTS, GRID_SIZE, solvers[0], comPositions, p);

    string test_name = "Compute equilibrium robustness ";
    for (int s = 1; s < N_SOLVERS; s++) {
      test_computeEquilibriumRobustness(solvers[0], solvers[s], comPositions, test_name + solvers[0]->getName(),
                                        test_name + solvers[s]->getName(), 1);
    }

    for (int s = 0; s < N_SOLVERS; s++) {
      test_computeEquilibriumRobustness_vs_checkEquilibrium(solvers[s], solver_PP, comPositions,
                                                            test_name + solvers[s]->getName(), "", 1);
    }

    const int N_TESTS_EXTREMUM = 100;
    centroidal_dynamics::Vector3 a0 = centroidal_dynamics::Vector3::Zero();
    a0.head<2>() = 0.5 * (com_LB + com_UB);
    double e_max;
    LP_status status = solvers[0]->computeEquilibriumRobustness(a0, e_max);
    if (status != LP_STATUS_OPTIMAL)
      SEND_ERROR_MSG(solvers[0]->getName() + " failed to compute robustness of com position " +
                     toString(a0.transpose()) + ", error code: " + toString(status));
    else {
      test_name = "EXTREMUM OVER LINE ";
      string test_name2 = "Compute equilibrium robustness ";
      for (int s = 1; s < N_SOLVERS; s++) {
        if (solvers[s]->getAlgorithm() != EQUILIBRIUM_ALGORITHM_LP2)
          test_findExtremumOverLine(solvers[s], solvers[0], a0, N_TESTS_EXTREMUM, e_max,
                                    test_name + solvers[s]->getName(), test_name2 + solvers[0]->getName(), 1);
      }
    }
  }

  getProfiler().report_all();

  cout << "*** END TEST WITH RANDOMLY GENERATED DATA ***\n";

  return 0;
}

/*
 * Copyright 2015, LAAS-CNRS
 * Author: Andrea Del Prete
 */

#ifdef CLP_FOUND
#include "ClpSimplex.hpp"
#include "CoinTime.hpp"
#include "CoinBuild.hpp"
#include "CoinModel.hpp"
#include <hpp/centroidal-dynamics/solver_LP_clp.hh>
#endif

#include <qpOASES.hpp>
#include <hpp/centroidal-dynamics/solver_LP_qpoases.hh>
#include <hpp/centroidal-dynamics/logger.hh>

#include <iostream>
#include <iomanip>

using namespace std;
using namespace centroidal_dynamics;
USING_NAMESPACE_QPOASES

#define EPS 1e-6

#ifdef CLP_FOUND
/** Example addRows.cpp */
void test_addRows() {
  const int N_ROWS_TO_ADD = 30000;
  try {
    // Empty model
    ClpSimplex model;

    // Objective - just nonzeros
    int objIndex[] = {0, 2};
    double objValue[] = {1.0, 4.0};
    // Upper bounds - as dense vector
    double upper[] = {2.0, COIN_DBL_MAX, 4.0};

    // Create space for 3 columns
    model.resize(0, 3);
    // Fill in
    int i;
    // Virtuous way
    // First objective
    for (i = 0; i < 2; i++) model.setObjectiveCoefficient(objIndex[i], objValue[i]);
    // Now bounds (lower will be zero by default but do again)
    for (i = 0; i < 3; i++) {
      model.setColumnLower(i, 0.0);
      model.setColumnUpper(i, upper[i]);
    }
    /*
         We could also have done in non-virtuous way e.g.
         double * objective = model.objective();
         and then set directly
       */
    // Faster to add rows all at once - but this is easier to show
    // Now add row 1 as >= 2.0
    int row1Index[] = {0, 2};
    double row1Value[] = {1.0, 1.0};
    model.addRow(2, row1Index, row1Value, 2.0, COIN_DBL_MAX);
    // Now add row 2 as == 1.0
    int row2Index[] = {0, 1, 2};
    double row2Value[] = {1.0, -5.0, 1.0};
    model.addRow(3, row2Index, row2Value, 1.0, 1.0);
    // solve
    model.dual();

    /*
         Adding one row at a time has a significant overhead so let's
         try a more complicated but faster way

         First time adding in 10000 rows one by one
       */
    model.allSlackBasis();
    ClpSimplex modelSave = model;
    double time1 = CoinCpuTime();
    int k;
    for (k = 0; k < N_ROWS_TO_ADD; k++) {
      int row2Index[] = {0, 1, 2};
      double row2Value[] = {1.0, -5.0, 1.0};
      model.addRow(3, row2Index, row2Value, 1.0, 1.0);
    }
    printf("Time for 10000 addRow is %g\n", CoinCpuTime() - time1);
    model.dual();
    model = modelSave;
    // Now use build
    CoinBuild buildObject;
    time1 = CoinCpuTime();
    for (k = 0; k < N_ROWS_TO_ADD; k++) {
      int row2Index[] = {0, 1, 2};
      double row2Value[] = {1.0, -5.0, 1.0};
      buildObject.addRow(3, row2Index, row2Value, 1.0, 1.0);
    }
    model.addRows(buildObject);
    printf("Time for 10000 addRow using CoinBuild is %g\n", CoinCpuTime() - time1);
    model.dual();
    model = modelSave;
    int del[] = {0, 1, 2};
    model.deleteRows(2, del);
    // Now use build +-1
    CoinBuild buildObject2;
    time1 = CoinCpuTime();
    for (k = 0; k < N_ROWS_TO_ADD; k++) {
      int row2Index[] = {0, 1, 2};
      double row2Value[] = {1.0, -1.0, 1.0};
      buildObject2.addRow(3, row2Index, row2Value, 1.0, 1.0);
    }
    model.addRows(buildObject2, true);
    printf("Time for 10000 addRow using CoinBuild+-1 is %g\n", CoinCpuTime() - time1);
    model.dual();
    model = modelSave;
    model.deleteRows(2, del);
    // Now use build +-1
    CoinModel modelObject2;
    time1 = CoinCpuTime();
    for (k = 0; k < N_ROWS_TO_ADD; k++) {
      int row2Index[] = {0, 1, 2};
      double row2Value[] = {1.0, -1.0, 1.0};
      modelObject2.addRow(3, row2Index, row2Value, 1.0, 1.0);
    }
    model.addRows(modelObject2, true);
    printf("Time for 10000 addRow using CoinModel+-1 is %g\n", CoinCpuTime() - time1);
    model.dual();
    model = ClpSimplex();
    // Now use build +-1
    CoinModel modelObject3;
    time1 = CoinCpuTime();
    for (k = 0; k < N_ROWS_TO_ADD; k++) {
      int row2Index[] = {0, 1, 2};
      double row2Value[] = {1.0, -1.0, 1.0};
      modelObject3.addRow(3, row2Index, row2Value, 1.0, 1.0);
    }
    model.loadProblem(modelObject3, true);
    printf("Time for 10000 addRow using CoinModel load +-1 is %g\n", CoinCpuTime() - time1);
    model.writeMps("xx.mps");
    model.dual();
    model = modelSave;
    // Now use model
    CoinModel modelObject;
    time1 = CoinCpuTime();
    for (k = 0; k < N_ROWS_TO_ADD; k++) {
      int row2Index[] = {0, 1, 2};
      double row2Value[] = {1.0, -5.0, 1.0};
      modelObject.addRow(3, row2Index, row2Value, 1.0, 1.0);
    }
    model.addRows(modelObject);
    printf("Time for 10000 addRow using CoinModel is %g\n", CoinCpuTime() - time1);
    model.dual();
    model.writeMps("b.mps");
    // Method using least memory - but most complicated
    time1 = CoinCpuTime();
    // Assumes we know exact size of model and matrix
    // Empty model
    ClpSimplex model2;
    {
      // Create space for 3 columns and 10000 rows
      int numberRows = N_ROWS_TO_ADD;
      int numberColumns = 3;
      // This is fully dense - but would not normally be so
      int numberElements = numberRows * numberColumns;
      // Arrays will be set to default values
      model2.resize(numberRows, numberColumns);
      double *elements = new double[numberElements];
      CoinBigIndex *starts = new CoinBigIndex[numberColumns + 1];
      int *rows = new int[numberElements];
      ;
      int *lengths = new int[numberColumns];
      // Now fill in - totally unsafe but ....
      // no need as defaults to 0.0 double * columnLower = model2.columnLower();
      double *columnUpper = model2.columnUpper();
      double *objective = model2.objective();
      double *rowLower = model2.rowLower();
      double *rowUpper = model2.rowUpper();
      // Columns - objective was packed
      for (k = 0; k < 2; k++) {
        int iColumn = objIndex[k];
        objective[iColumn] = objValue[k];
      }
      for (k = 0; k < numberColumns; k++) columnUpper[k] = upper[k];
      // Rows
      for (k = 0; k < numberRows; k++) {
        rowLower[k] = 1.0;
        rowUpper[k] = 1.0;
      }
      // Now elements
      double row2Value[] = {1.0, -5.0, 1.0};
      CoinBigIndex put = 0;
      for (k = 0; k < numberColumns; k++) {
        starts[k] = put;
        lengths[k] = numberRows;
        double value = row2Value[k];
        for (int i = 0; i < numberRows; i++) {
          rows[put] = i;
          elements[put] = value;
          put++;
        }
      }
      starts[numberColumns] = put;
      // assign to matrix
      CoinPackedMatrix *matrix = new CoinPackedMatrix(true, 0.0, 0.0);
      matrix->assignMatrix(true, numberRows, numberColumns, numberElements, elements, rows, starts, lengths);
      ClpPackedMatrix *clpMatrix = new ClpPackedMatrix(matrix);
      model2.replaceMatrix(clpMatrix, true);
      printf("Time for 10000 addRow using hand written code is %g\n", CoinCpuTime() - time1);
      // If matrix is really big could switch off creation of row copy
      // model2.setSpecialOptions(256);
    }
    model2.dual();
    model2.writeMps("a.mps");
    // Print column solution
    int numberColumns = model.numberColumns();

    // Alternatively getColSolution()
    double *columnPrimal = model.primalColumnSolution();
    // Alternatively getReducedCost()
    double *columnDual = model.dualColumnSolution();
    // Alternatively getColLower()
    double *columnLower = model.columnLower();
    // Alternatively getColUpper()
    double *columnUpper = model.columnUpper();
    // Alternatively getObjCoefficients()
    double *columnObjective = model.objective();

    int iColumn;

    std::cout << "               Primal          Dual         Lower         Upper          Cost" << std::endl;

    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      double value;
      std::cout << std::setw(6) << iColumn << " ";
      value = columnPrimal[iColumn];
      if (fabs(value) < 1.0e5)
        std::cout << setiosflags(std::ios::fixed | std::ios::showpoint) << std::setw(14) << value;
      else
        std::cout << setiosflags(std::ios::scientific) << std::setw(14) << value;
      value = columnDual[iColumn];
      if (fabs(value) < 1.0e5)
        std::cout << setiosflags(std::ios::fixed | std::ios::showpoint) << std::setw(14) << value;
      else
        std::cout << setiosflags(std::ios::scientific) << std::setw(14) << value;
      value = columnLower[iColumn];
      if (fabs(value) < 1.0e5)
        std::cout << setiosflags(std::ios::fixed | std::ios::showpoint) << std::setw(14) << value;
      else
        std::cout << setiosflags(std::ios::scientific) << std::setw(14) << value;
      value = columnUpper[iColumn];
      if (fabs(value) < 1.0e5)
        std::cout << setiosflags(std::ios::fixed | std::ios::showpoint) << std::setw(14) << value;
      else
        std::cout << setiosflags(std::ios::scientific) << std::setw(14) << value;
      value = columnObjective[iColumn];
      if (fabs(value) < 1.0e5)
        std::cout << setiosflags(std::ios::fixed | std::ios::showpoint) << std::setw(14) << value;
      else
        std::cout << setiosflags(std::ios::scientific) << std::setw(14) << value;

      std::cout << std::endl;
    }
    std::cout << "--------------------------------------" << std::endl;
    // Test CoinAssert
    std::cout << "If Clp compiled without NDEBUG below should give assert, if with NDEBUG or COIN_ASSERT CoinError"
              << std::endl;
    model = modelSave;
    model.deleteRows(2, del);
    // Deliberate error
    model.deleteColumns(1, del + 2);
    // Now use build +-1
    CoinBuild buildObject3;
    time1 = CoinCpuTime();
    for (k = 0; k < N_ROWS_TO_ADD; k++) {
      int row2Index[] = {0, 1, 2};
      double row2Value[] = {1.0, -1.0, 1.0};
      buildObject3.addRow(3, row2Index, row2Value, 1.0, 1.0);
    }
    model.addRows(buildObject3, true);
  } catch (CoinError e) {
    e.print();
    if (e.lineNumber() >= 0) std::cout << "This was from a CoinAssert" << std::endl;
  }
}

void test_small_LP() {
  ClpSimplex model;
  //  model.setLogLevel(1);

  // Objective - just nonzeros
  int objIndex[] = {0, 2};
  double objValue[] = {1.0, 4.0};
  // Upper bounds - as dense vector
  double upper[] = {2.0, COIN_DBL_MAX, 4.0};

  // Create space for 3 columns
  model.resize(0, 3);

  // Can also use maximumIterations
  int integerValue;
  model.getIntParam(ClpMaxNumIteration, integerValue);
  cout << "Value of ClpMaxNumIteration is " << integerValue << endl;
  model.setMaximumIterations(integerValue);

  // Fill in
  int i;
  // Virtuous way
  // First objective
  for (i = 0; i < 2; i++) model.setObjectiveCoefficient(objIndex[i], objValue[i]);
  // Now bounds (lower will be zero by default but do again)
  for (i = 0; i < 3; i++) {
    model.setColumnLower(i, 0.0);
    model.setColumnUpper(i, upper[i]);
  }
  /*
      We could also have done in non-virtuous way e.g.
      double * objective = model.objective();
      and then set directly
    */
  // Faster to add rows all at once - but this is easier to show
  // Now add row 1 as >= 2.0
  int row1Index[] = {0, 2};
  double row1Value[] = {1.0, 1.0};
  model.addRow(2, row1Index, row1Value, 2.0, COIN_DBL_MAX);
  // Now add row 2 as == 1.0
  int row2Index[] = {0, 1, 2};
  double row2Value[] = {1.0, -5.0, 1.0};
  model.addRow(3, row2Index, row2Value, 1.0, 1.0);

  int n = model.getdimVarXs();
  int m = model.getNumRows();
  cout << "Problem has " << n << " variables and " << m << " constraints.\n";

  // solve
  model.dual();

  // Check the solution
  if (model.isProvenOptimal()) {
    cout << "Found optimal solution!" << endl;
    cout << "Objective value is " << model.getObjValue() << endl;
    cout << "Model status is " << model.status() << " after " << model.numberIterations()
         << " iterations - objective is " << model.objectiveValue() << endl;
    const double *solution;
    solution = model.getColSolution();
    // We could then print the solution or examine it.
    cout << "Solution is: ";
    for (int i = 0; i < n; i++) cout << solution[i] << ", ";
    cout << endl;
  } else
    cout << "Didn’t find optimal solution." << endl;
}
#endif

int main() {
  cout << "Test LP Solvers (1 means ok, 0 means error)\n\n";

  {
    cout << "TEST QP OASES ON A SMALL 2-VARIABLE LP";
    /* Setup data of first LP. */
    real_t A[1 * 2] = {1.0, 1.0};
    real_t g[2] = {1.5, 1.0};
    real_t lb[2] = {0.5, -2.0};
    real_t ub[2] = {5.0, 2.0};
    real_t lbA[1] = {-1.0};
    real_t ubA[1] = {2.0};

    /* Setting up QProblem object with zero Hessian matrix. */
    QProblem example(2, 1, HST_ZERO);

    Options options;
    // options.setToMPC();
    example.setOptions(options);

    /* Solve first LP. */
    int nWSR = 10;
    int res = example.init(0, g, A, lb, ub, lbA, ubA, nWSR, 0);
    if (res == 0)
      cout << "[INFO] LP solved correctly\n";
    else
      cout << "[ERROR] QpOases could not solve the LP problem, error code: " << res << endl;
  }

  {
    cout << "\nTEST READ-WRITE METHODS OF SOLVER_LP_ABSTRACT\n";
    Solver_LP_abstract *solverOases = Solver_LP_abstract::getNewSolver(SOLVER_LP_QPOASES);
    const int n = 3;
    const int m = 4;
    const char *filename = "small_3_x_4_LP.dat";
    VectorX c = VectorX::Random(n);
    VectorX lb = -100 * VectorX::Ones(n);
    VectorX ub = 100 * VectorX::Ones(n);
    MatrixXX A = MatrixXX::Random(m, n);
    VectorX Alb = -100 * VectorX::Ones(m);
    VectorX Aub = 100 * VectorX::Ones(m);
    if (!solverOases->writeLpToFile(filename, c, lb, ub, A, Alb, Aub)) {
      SEND_ERROR_MSG("Error while writing LP to file");
      return -1;
    }
    VectorX c2, lb2, ub2, Alb2, Aub2;
    MatrixXX A2;
    if (!solverOases->readLpFromFile(filename, c2, lb2, ub2, A2, Alb2, Aub2)) {
      SEND_ERROR_MSG("Error while reading LP from file");
      return -1;
    }

    cout << "Check number of variables: " << (c.size() == c2.size()) << endl;
    cout << "Check number of constraints: " << (A.rows() == A2.rows()) << endl;
    cout << "Check gradient vector c: " << c.isApprox(c2) << endl;
    cout << "Check lower bound vector lb: " << lb.isApprox(lb2) << endl;
    cout << "Check upper bound vector ub: " << ub.isApprox(ub2) << endl;
    cout << "Check constraint matrix A: " << A.isApprox(A2) << endl;
    cout << "Check constraint lower bound vector Alb: " << Alb.isApprox(Alb2) << endl;
    cout << "Check constraint upper bound vector Aub: " << Aub.isApprox(Aub2) << endl;
  }

  {
    cout << "\nTEST QP OASES ON SOME LP PROBLEMS\n";
    string file_path = "../test_data/";
    Solver_LP_abstract *solverOases = Solver_LP_abstract::getNewSolver(SOLVER_LP_QPOASES);
    const int PROBLEM_NUMBER = 14;
    string problem_filenames[PROBLEM_NUMBER] = {
        "DLP_findExtremumOverLine20151103_112611",     "DLP_findExtremumOverLine20151103_115627",
        "DLP_findExtremumOverLine20151103_014022",     "DLP_findExtremumOverLine_32_generators",
        "DLP_findExtremumOverLine_64_generators",      "DLP_findExtremumOverLine_128_generators",
        "DLP_findExtremumOverLine_128_generators_bis", "LP_findExtremumOverLine20151103_112610",
        "LP_findExtremumOverLine20151103_112611",      "LP_findExtremumOverLine20151103_014022",
        "LP_findExtremumOverLine_32_generators",       "LP_findExtremumOverLine_64_generators",
        "LP_findExtremumOverLine_128_generators",      "LP_findExtremumOverLine_128_generators_bis"};
    VectorX c, lb, ub, Alb, Aub, realSol, sol;
    MatrixXX A;
    for (int i = 0; i < PROBLEM_NUMBER; i++) {
      string &problem_filename = problem_filenames[i];
      if (!solverOases->readLpFromFile(file_path + problem_filename + ".dat", c, lb, ub, A, Alb, Aub)) {
        SEND_ERROR_MSG("Error while reading LP from file " + problem_filename);
        return -1;
      }
      string solution_filename = problem_filename + "_solution";
      if (!readMatrixFromFile(file_path + solution_filename + ".dat", realSol)) {
        SEND_ERROR_MSG("Error while reading LP solution from file " + solution_filename);
        // return -1;
      }
      sol.resize(c.size());
      solverOases->solve(c, lb, ub, A, Alb, Aub, sol);
      if (sol.isApprox(realSol, EPS)) {
        cout << "[INFO] Solution of problem " << problem_filename << " (" << c.size() << " var, " << A.rows()
             << " constr) is equal to the expected value!\n";
      } else {
        if (fabs(c.dot(sol) - c.dot(realSol)) < EPS)
          cout << "[WARNING] Solution of problem " << problem_filename << " (" << c.size() << " var, " << A.rows()
               << " constr) is different from expected but it has the same cost\n";
        else {
          cout << "[ERROR] Solution of problem " << problem_filename << " (" << c.size() << " var, " << A.rows()
               << " constr) is different from the expected value:\n";
          cout << "\tSolution found    " << sol.transpose() << endl;
          cout << "\tExpected solution " << realSol.transpose() << endl;
          cout << "\tCost found    " << (c.dot(sol)) << endl;
          cout << "\tCost expected " << (c.dot(realSol)) << endl;
        }
      }
    }

    return 0;
  }

#ifdef CLP_FOUND
  test_addRows();
  test_small_LP();

  Solver_LP_abstract *solver = Solver_LP_abstract::getNewSolver(SOLVER_LP_CLP);
  Vector3 c, lb, ub, x;
  MatrixXX A(2, 3);
  Vector2 Alb, Aub;
  c << 1.0, 0.0, 4.0;
  lb << 0.0, 0.0, 0.0;
  ub << 2.0, COIN_DBL_MAX, 4.0;
  A << 1.0, 0.0, 1.0, 1.0, -5.0, 1.0;
  Alb << 2.0, 1.0;
  Aub << COIN_DBL_MAX, 1.0;
  if (solver->solve(c, lb, ub, A, Alb, Aub, x) == LP_STATUS_OPTIMAL) {
    cout << "solver_LP_clp solved the problem\n";
    cout << "The solution is " << x.transpose() << endl;
  } else
    cout << "solver_LP_clp failed to solve the problem\n";
#endif

  return 1;
}

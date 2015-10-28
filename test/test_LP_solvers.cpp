/*
 * Copyright 2015, LAAS-CNRS
 * Author: Andrea Del Prete
 */

#include "ClpSimplex.hpp"
#include "CoinTime.hpp"
#include "CoinBuild.hpp"
#include "CoinModel.hpp"

#include <robust-equilibrium-lib/solver_LP_clp.hh>

#include <iostream>
#include <iomanip>

using namespace std;
using namespace robust_equilibrium;

/** Example addRows.cpp */
void test_addRows()
{
  const int N_ROWS_TO_ADD = 30000;
  try
  {
    // Empty model
    ClpSimplex  model;

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
    for (i = 0; i < 2; i++)
      model.setObjectiveCoefficient(objIndex[i], objValue[i]);
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
    model.addRow(2, row1Index, row1Value,
                 2.0, COIN_DBL_MAX);
    // Now add row 2 as == 1.0
    int row2Index[] = {0, 1, 2};
    double row2Value[] = {1.0, -5.0, 1.0};
    model.addRow(3, row2Index, row2Value,
                 1.0, 1.0);
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
      model.addRow(3, row2Index, row2Value,
                   1.0, 1.0);
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
      buildObject.addRow(3, row2Index, row2Value,
                         1.0, 1.0);
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
      buildObject2.addRow(3, row2Index, row2Value,
                          1.0, 1.0);
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
      modelObject2.addRow(3, row2Index, row2Value,
                          1.0, 1.0);
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
      modelObject3.addRow(3, row2Index, row2Value,
                          1.0, 1.0);
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
      modelObject.addRow(3, row2Index, row2Value,
                         1.0, 1.0);
    }
    model.addRows(modelObject);
    printf("Time for 10000 addRow using CoinModel is %g\n", CoinCpuTime() - time1);
    model.dual();
    model.writeMps("b.mps");
    // Method using least memory - but most complicated
    time1 = CoinCpuTime();
    // Assumes we know exact size of model and matrix
    // Empty model
    ClpSimplex  model2;
    {
      // Create space for 3 columns and 10000 rows
      int numberRows = N_ROWS_TO_ADD;
      int numberColumns = 3;
      // This is fully dense - but would not normally be so
      int numberElements = numberRows * numberColumns;
      // Arrays will be set to default values
      model2.resize(numberRows, numberColumns);
      double * elements = new double [numberElements];
      CoinBigIndex * starts = new CoinBigIndex [numberColumns+1];
      int * rows = new int [numberElements];;
      int * lengths = new int[numberColumns];
      // Now fill in - totally unsafe but ....
      // no need as defaults to 0.0 double * columnLower = model2.columnLower();
      double * columnUpper = model2.columnUpper();
      double * objective = model2.objective();
      double * rowLower = model2.rowLower();
      double * rowUpper = model2.rowUpper();
      // Columns - objective was packed
      for (k = 0; k < 2; k++) {
        int iColumn = objIndex[k];
        objective[iColumn] = objValue[k];
      }
      for (k = 0; k < numberColumns; k++)
        columnUpper[k] = upper[k];
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
      CoinPackedMatrix * matrix = new CoinPackedMatrix(true, 0.0, 0.0);
      matrix->assignMatrix(true, numberRows, numberColumns, numberElements,
                           elements, rows, starts, lengths);
      ClpPackedMatrix * clpMatrix = new ClpPackedMatrix(matrix);
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
    double * columnPrimal = model.primalColumnSolution();
    // Alternatively getReducedCost()
    double * columnDual = model.dualColumnSolution();
    // Alternatively getColLower()
    double * columnLower = model.columnLower();
    // Alternatively getColUpper()
    double * columnUpper = model.columnUpper();
    // Alternatively getObjCoefficients()
    double * columnObjective = model.objective();

    int iColumn;

    std::cout << "               Primal          Dual         Lower         Upper          Cost"
              << std::endl;

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
    std::cout << "If Clp compiled without NDEBUG below should give assert, if with NDEBUG or COIN_ASSERT CoinError" << std::endl;
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
      buildObject3.addRow(3, row2Index, row2Value,
                          1.0, 1.0);
    }
    model.addRows(buildObject3, true);
  } catch (CoinError e) {
    e.print();
    if (e.lineNumber() >= 0)
      std::cout << "This was from a CoinAssert" << std::endl;
  }
}

void test_small_LP()
{
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
  for (i = 0; i < 2; i++)
    model.setObjectiveCoefficient(objIndex[i], objValue[i]);
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
  model.addRow(2, row1Index, row1Value,
               2.0, COIN_DBL_MAX);
  // Now add row 2 as == 1.0
  int row2Index[] = {0, 1, 2};
  double row2Value[] = {1.0, -5.0, 1.0};
  model.addRow(3, row2Index, row2Value,
               1.0, 1.0);

  int n = model.getNumCols();
  int m = model.getNumRows();
  cout<<"Problem has "<<n<<" variables and "<<m<<" constraints.\n";

  // solve
  model.dual();

  // Check the solution
  if ( model.isProvenOptimal() )
  {
    cout << "Found optimal solution!" << endl;
    cout << "Objective value is " << model.getObjValue() << endl;
    cout << "Model status is " << model.status() << " after "
              << model.numberIterations() << " iterations - objective is "
              << model.objectiveValue() << endl;
    const double *solution;
    solution = model.getColSolution();
    // We could then print the solution or examine it.
    cout<<"Solution is: ";
    for(int i=0; i<n; i++)
      cout<<solution[i]<<", ";
    cout<<endl;
  }
  else
    cout << "Didn’t find optimal solution." << endl;
}

int main()
{
  cout <<"Test LP Solvers\n";

//  test_addRows();
  test_small_LP();

  Solver_LP_abstract *solver = Solver_LP_abstract::getNewSolver(SOLVER_LP_CLP);
  Vector3 c, lb, ub, x;
  MatrixXX A(2,3);
  Vector2 Alb, Aub;
  c << 1.0, 0.0, 4.0;
  lb << 0.0, 0.0, 0.0;
  ub << 2.0, COIN_DBL_MAX, 4.0;
  A << 1.0, 0.0, 1.0,
      1.0, -5.0, 1.0;
  Alb << 2.0, 1.0;
  Aub << COIN_DBL_MAX, 1.0;
  if(solver->solve(c, lb, ub, A, Alb, Aub, x)==LP_STATUS_OPTIMAL)
  {
    cout<<"solver_LP_clp solved the problem\n";
    cout<<"The solution is "<<x.transpose()<<endl;
  }
  else
    cout<<"solver_LP_clp failed to solve the problem\n";

//  char x[81];
//  int iRow;
//  // get row copy
//  CoinPackedMatrix rowCopy = *model.matrix();
//  rowCopy.reverseOrdering();
//  const int * column = rowCopy.getIndices();
//  const int * rowLength = rowCopy.getVectorLengths();
//  const CoinBigIndex * rowStart = rowCopy.getVectorStarts();
//  x[n] = '\0';
//  for (iRow = 0; iRow < m; iRow++) {
//    memset(x, ' ', n);
//    for (int k = rowStart[iRow]; k < rowStart[iRow] + rowLength[iRow]; k++) {
//      int iColumn = column[k];
//      x[iColumn] = 'x';
//    }
//    cout<<x<<endl;
//  }
//  cout<<endl;

//  // Now matrix
//  CoinPackedMatrix * matrix = model.matrix();
//  const double * element = matrix->getElements();
//  const int * row = matrix->getIndices();
//  const int * start = matrix->getVectorStarts();
//  const int * length = matrix->getVectorLengths();
//  for (int iColumn = 0; iColumn < n; iColumn++)
//  {
//       cout << "Column " << iColumn;
//       int j;
//       for (j = start[iColumn]; j < start[iColumn] + length[iColumn]; j++)
//            cout << " ( " << row[j] << ", " << element[j] << ")";
//       cout << endl;
//  }

  return 1;
}

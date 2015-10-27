#include <vector>
#include <iostream>
#include <robust-equilibrium-lib/static_equilibrium.hh>

using namespace robust_equilibrium;
using namespace Eigen;
using namespace std;



int main()
{
  init_library();
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
  CONTACT_POINT_UPPER_BOUNDS << 0.5,  0.5,  0.0;
  double gamma = atan(mu);   // half friction cone angle
  RVector3 RPY_LOWER_BOUNDS, RPY_UPPER_BOUNDS;
  RPY_LOWER_BOUNDS << -0.0*gamma, -0.0*gamma, -M_PI;
  RPY_UPPER_BOUNDS << +0.0*gamma, +0.0*gamma, +M_PI;
  double X_MARG = 0.07;
  double Y_MARG = 0.07;
  const int GRID_SIZE = 25;

  StaticEquilibrium solver(mass, generatorsPerContact, SOLVER_LP_QPOASES);
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
//    while(true) // generate contact position
//    {
//      uniform(CONTACT_POINT_LOWER_BOUNDS, CONTACT_POINT_UPPER_BOUNDS, contact_pos.row(i));
//      if(i==0)
//        break;
//      collision = false;
//      for(unsigned int j=0; j<i-1; j++)
//        if((contact_pos.row(i)-contact_pos.row(j)).norm() < MIN_FEET_DISTANCE)
//          collision = true;
//      if(collision==false)
//        break;
//    }

//    // generate contact orientation
//    uniform(RPY_LOWER_BOUNDS, RPY_UPPER_BOUNDS, contact_rpy.row(i));

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

  if(!solver.setNewContacts(p, N, frictionCoefficients, STATIC_EQUILIBRIUM_ALGORITHM_PP))
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
      if(solver.checkRobustEquilibrium(com, 0.0))
      {
        equilibrium(i,j) = 1.0;
        printf("1 ");
      }
      else
      {
        equilibrium(i,j) = 0.0;
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

  cout<<"Max dist between contact points and grid points "<<minDistContactPoint.maxCoeff()<<"\n";

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

  return ret;
}

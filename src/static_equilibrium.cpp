#include <robust-equilibrium-lib/static_equilibrium.hh>
#include <iostream>
#include <vector>

using namespace std;

namespace robust_equilibrium
{


StaticEquilibrium::StaticEquilibrium(double mass, unsigned int generatorsPerContact, SolverLP solver_type)
{
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

  if(!computePolytopeProjection(m_G_centr))
    return false;

  m_HD = m_H * m_D;
  m_Hd = m_H * m_d;

  return true;
}

double StaticEquilibrium::computeEquilibriumRobustness(Cref_vector2 com)
{
  return 0.0;
}

bool StaticEquilibrium::checkRobustEquilibrium(Cref_vector2 com, double e_max)
{
  if(e_max!=0.0)
  {
    std::cerr<<"ERROR checkRobustEquilibrium with e_max!=0 not implemented yet!";
    return false;
  }

  VectorX res = m_HD * com + m_Hd;
  for(long i=0; i<res.size(); i++)
    if(res(i)>0.0)
      return false;
  return true;
}

double StaticEquilibrium::findExtremumOverLine(Cref_vector2 a, double b, double e_max)
{
  return 0.0;
}

double StaticEquilibrium::findExtremumInDirection(Cref_vector2 direction, double e_max)
{
  return 0.0;
}

bool StaticEquilibrium::computePolytopeProjection(Cref_matrix6X v)
{
  dd_MatrixPtr V = cone_span_eigen_to_cdd(v.transpose());
  dd_ErrorType error = dd_NoError;
  dd_PolyhedraPtr H_= dd_DDMatrix2Poly(V, &error);
  if(error != dd_NoError)
  {
//    if(dd_debug)
    cout << ("numerical instability in cddlib. ill formed polytope") << endl;
    return false;
  }

  dd_MatrixPtr b_A = dd_CopyInequalities(H_);
  // get equalities and add them as complementary inequality constraints
  std::vector<long> eq_rows;
  for(long elem=1;elem<=(long)(b_A->linset[0]);++elem)
  {
    if (set_member(elem,b_A->linset))
      eq_rows.push_back(elem);
  }
  int rowsize = (int)b_A->rowsize;
//  cout<<"Inequality matrix has "<<rowsize<<" rows and "<<b_A->colsize-1<<" columns\n";
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

  return true;
}

} // end namespace robust_equilibrium

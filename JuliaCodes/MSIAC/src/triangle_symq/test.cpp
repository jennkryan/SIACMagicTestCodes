# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <cstring>

using namespace std;

# include "triangle_symq_rule.hpp"

int main ( );
void test01 ( );
void test02 ( int degree, int numnodes, double vert1[], double vert2[], 
  double vert3[], string header  );
void test03 ( int degree, int numnodes, double vert1[], double vert2[], 
  double vert3[], string header );
void test04 ( int degree, int numnodes, double vert1[], double vert2[], 
  double vert3[], string header );
void test05 ( int degree, int numnodes, double vert1[], double vert2[], 
  double vert3[] );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TRIANGLE_SYMQ_RULE_TEST.
//
//  Discussion:
//
//    TRIANGLE_SYMQ_RULE_TEST tests the TRIANGLE_SYMQ_RULE library.
//
//  Licensing:
//
//    This code is distributed under the GNU GPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
{
  int degree;
  string header;
  int itype;
  int numnodes;
  double vert1[2];
  double vert2[2];
  double vert3[2];

  timestamp ( );
  cout << "\n";
  cout << "TRIANGLE_SYMQ_RULE_TEST\n";
  cout << "  C++ version\n";
  cout << "  Test the TRIANGLE_SYMQ_RULE library.\n";

  cout << "\n";
  cout << "  Region is the simplex (0,0),(1,0),(0,1).\n";
  vert1[0] = 1.0;
  vert1[1] = 0.0;
  vert2[0] = 4.0;
  vert2[1] = 4.0;
  vert3[0] = 0.0;
  vert3[1] = 3.0;
  header = "user08";
  degree = 3;
  cout << "\n";
  cout << "  Triangle:\n";
  cout << "\n";
  cout << vert1[0] << "  " << vert1[1] << "\n";
  cout << vert2[0] << "  " << vert2[1] << "\n";
  cout << vert3[0] << "  " << vert3[1] << "\n";

//
//  Determine the size of the rule.  Retrieve a rule and print it.
//
    numnodes = rule_full_size ( degree );
    test02 ( degree, numnodes, vert1, vert2, vert3, header );

    cout << "\n";
    cout << "TRIANGLE_SYMQ_RULE_TEST\n";
    cout << "  Normal end of execution.\n";
    cout << "\n";
    timestamp ( );  
    return 0;
}
//****************************************************************************80

void test02 ( int degree, int numnodes, double vert1[], double vert2[], 
  double vert3[], string header )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 calls TRIASYMQ for a quadrature rule of given order and region.
//
//  Licensing:
//
//    This code is distributed under the GNU GPL license.
//
//  Modified:
//
//    28 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int DEGREE, the desired total polynomial degree exactness
//    of the quadrature rule.  0 <= DEGREE <= 50.
//
//    Input, int NUMNODES, the number of nodes to be used by the rule.
//
//    Input, double VERT1[2], VERT2[2], VERT3[2], the
//    vertices of the triangle.
//
{
  double area;
  double d;
  int j;
  double *rnodes;
  double *weights;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  Symmetric quadrature rule for a triangle.\n";
  cout << "  Polynomial exactness degree DEGREE = " << degree << "\n";

  area = triangle_area ( vert1, vert2, vert3 );
//
//  Retrieve and print a symmetric quadrature rule.
//
  rnodes = new double[2*numnodes];
  weights = new double[numnodes];

  triasymq ( degree, vert1, vert2, vert3, rnodes, weights, numnodes );

  cout << "\n";
  cout << "  NUMNODES = " << numnodes << "\n";

  cout << "\n";
  cout << "     J      W               X               Y\n";
  cout << "\n";
  for ( j = 0; j < numnodes; j++ )
  {
    cout << j << "  "
         << weights[j] << "  "
         << rnodes[0+j*2] << "  "
         << rnodes[1+j*2] << "\n";
  }

  d = r8vec_sum ( numnodes, weights );

  cout << "   Sum  " << d << "\n";
  cout << "  Area  " << area << "\n";    
  triasymq_gnuplot ( vert1, vert2, vert3, numnodes, rnodes, header);
  
  delete [] rnodes;
  delete [] weights;

  return;
}
//****************************************************************************80

void test03 ( int degree, int numnodes, double vert1[], double vert2[], 
  double vert3[], string header )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 calls TRIASYMQ_GNUPLOT to generate graphics files.
//
//  Licensing:
//
//    This code is distributed under the GNU GPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int DEGREE, the desired total polynomial degree exactness
//    of the quadrature rule.  0 <= DEGREE <= 50.
//
//    Input, int NUMNODES, the number of nodes to be used by the rule.
//
//    Input, double VERT1[2], VERT2[2], VERT3[2], the
//    vertices of the triangle.
//
//    Input, string HEADER, an identifier for the graphics filenames.
//
{
  double *rnodes;
  double *weights;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  TRIASYMQ_GNUPLOT creates gnuplot graphics files.\n";
  cout << "  Polynomial exactness degree DEGREE = " << degree << "\n";

  rnodes = new double[2*numnodes];
  weights = new double[numnodes];

  triasymq ( degree, vert1, vert2, vert3, rnodes, weights, numnodes );

  cout << "  Number of nodes = " << numnodes << "\n";

  triasymq_gnuplot ( vert1, vert2, vert3, numnodes, rnodes, header );

  delete [] rnodes;
  delete [] weights;

  return;
}
//****************************************************************************80

void test04 ( int degree, int numnodes, double vert1[], double vert2[], 
  double vert3[], string header )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 gets a rule and writes it to a file.
//
//  Licensing:
//
//    This code is distributed under the GNU GPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int DEGREE, the desired total polynomial degree exactness
//    of the quadrature rule.  0 <= DEGREE <= 50.
//
//    Input, int NUMNODES, the number of nodes to be used by the rule.
//
//    Input, double VERT1[2], VERT2[2], VERT3[2], the
//    vertices of the triangle.
//
//    Input, string HEADER, an identifier for the filenames.
//
{
  int j;
  double *rnodes;
  ofstream rule_unit;
  string rule_filename;
  double *weights;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  Get a quadrature rule for a triangle.\n";
  cout << "  Then write it to a file.\n";
  cout << "  Polynomial exactness degree DEGREE = " << degree << "\n";
//
//  Retrieve a symmetric quadrature rule.
//
  rnodes = new double[2*numnodes];
  weights = new double[numnodes];

  triasymq ( degree, vert1, vert2, vert3, rnodes, weights, numnodes );
//
//  Write the points and weights to a file.
//
  rule_filename = header + ".txt";
 
  rule_unit.open ( rule_filename.c_str ( ) );
  for ( j = 0; j < numnodes; j++ )
  {
    rule_unit << rnodes[0+j*2] << "  "
              << rnodes[1+j*2] << "  "
              << weights[j] << "\n";
  }
  rule_unit.close ( );
  cout << "\n";
  cout << "  Quadrature rule written to file '" << rule_filename << "'\n";

  delete [] rnodes;
  delete [] weights;

  return;
}
//****************************************************************************80

void test05 ( int degree, int numnodes, double vert1[], double vert2[], 
  double vert3[] )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 calls TRIASYMQ for a quadrature rule of given order and region.
//
//  Licensing:
//
//    This code is distributed under the GNU GPL license.
//
//  Modified:
//
//    28 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int DEGREE, the desired total polynomial degree 
//    exactness of the quadrature rule.  0 <= DEGREE <= 50.
//
//    Input, int NUMNODES, the number of nodes to be used by the rule.
//
//    Input, double VERT1[2], VERT2[2], VERT3[2], the
//    vertices of the triangle.
//
{
  double area;
  double d;
  int i;
  int j;
  int npols;
  double *pols;
  double *r;
  double *rints;
  double *rnodes;
  double scale;
  double *weights;
  double z[2];

  cout << "\n";
  cout << "TEST05\n";
  cout << "  Compute a quadrature rule for a triangle.\n";
  cout << "  Check it by integrating orthonormal polynomials.\n";
  cout << "  Polynomial exactness degree DEGREE = " << degree << "\n";

  area = triangle_area ( vert1, vert2, vert3 );
//
//  Retrieve a symmetric quadrature rule.
//
  rnodes = new double[2*numnodes];
  weights = new double[numnodes];

  triasymq ( degree, vert1, vert2, vert3, rnodes, weights, numnodes );
//
//  Construct the matrix of values of the orthogonal polynomials
//  at the user-provided nodes
//
  npols = ( degree + 1 ) * ( degree + 2 ) / 2;
  rints = new double[npols];

  for ( j = 0; j < npols; j++ )
  {
    rints[j] = 0.0;
  }

  for ( i = 0; i < numnodes; i++ )
  {
    z[0] = rnodes[0+i*2];
    z[1] = rnodes[1+i*2];
    r = triangle_to_ref ( vert1, vert2, vert3, z );
    pols = ortho2eva ( degree, r );
    for ( j = 0; j < npols; j++ )
    {
      rints[j] = rints[j] + weights[i] * pols[j];
    }
    delete [] pols;
    delete [] r;
  }

  scale = sqrt ( sqrt ( 3.0 ) ) / sqrt ( area );

  for ( j = 0; j < npols; j++ )
  {
    rints[j] = rints[j] * scale;
  }

  d = pow ( rints[0] - sqrt ( area ), 2 );
  for ( j = 1; j < npols; j++ )
  {
    d = d + rints[j] * rints[j];
  }
  d = sqrt ( d ) / ( double ) ( npols );

  cout << "\n";
  cout << "  RMS integration error = " << d << "\n";

  delete [] rints;
  delete [] rnodes;
  delete [] weights;

  return;
}


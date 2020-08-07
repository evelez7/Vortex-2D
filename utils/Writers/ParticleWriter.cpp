#include "ParticleWriter.H"
#include "Proto_VisitWriter.H"
#include "ParticleSet.H"
#include <iostream>
#include <array>
#include <algorithm>
#include <tuple>
#include <memory>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
static FILE *fp = NULL;
static int useBinary = 0;
static int numInColumn = 0;

using namespace std;

array<array<double, DIM>, DIM> get_transpose(const array<array<double, DIM>, DIM>& matrix) {
  array<array<double, DIM>, DIM> transpose;

  for (int i=0; i < DIM; ++i)
  {
    for (int j=0; j < DIM; ++j)
    {
      transpose[i][j] = matrix[j][i];
    }
  }
  return transpose;
}

array<array<double, DIM>, DIM> multiply_matrices(const array<array<double, DIM>, DIM>& first, const array<array<double, DIM>, DIM>& second)
{
  array<array<double, DIM>, DIM> product;
  for (int i=0; i<DIM;++i)
  {
    for (int j=0; j<DIM;++j)
    {
      product[i][j] = 0;
      for (int k=0; k<DIM;++k)
      {
        product[i][j] += first[i][k] * second[k][j];
      }
    }
  }
  return product;
}

array<double, DIM> multiply_matrix_by_vector(const array<array<double, DIM>, DIM>& matrix, const array<double, DIM>& vec)
{
  array<double, DIM> new_vector;
  for (int i=0;i<DIM;++i) {
    double sum = 0;
    for (int j=0; j< matrix[i].size(); ++j) {
      sum +=  matrix[i][j] * vec[i];
    }
    new_vector[i] = sum;
  }
  return new_vector;
}

// only works for
array<array<double, DIM>, DIM> get_inverse(const array<array<double, DIM>, DIM>& matrix)
{
  array<array<double, DIM>, DIM> matrix_inverse;
  double determinant = 1/((matrix[0][0]*matrix[1][1]) - (matrix[0][1]*matrix[1][0]));
  for (int i=0; i<DIM; ++i)
  {
    for (int j=0; j<DIM; ++j)
    {
     matrix_inverse[i][j] = matrix[i][j] * determinant;
    }
  }
  return matrix_inverse;
}

void print_matrix(const array<array<double, DIM>, DIM>& matrix)
{
  for (int i=0; i<DIM; ++i)
  {
    for (int j=0; j<DIM; ++j)
    {
      cout << matrix[i][j] << " ";
    }
    cout << endl;
  }
  cout << endl;
}

// arr[0] = a, arr[1] = b, arr[2] = c, x=-b +- sqrt(b^2 -4ac) all over 2a
array<double, DIM> get_roots(const array<double, 3>& characteristic_polynomial)
{
  double a = characteristic_polynomial[0];
  double b = characteristic_polynomial[1];
  double c = characteristic_polynomial[2];
  double discriminant = pow(b, 2.0) - (4 * a * c);
  array<double, DIM> roots;
  if (discriminant > 0)
  {
    roots[0] = (-b + sqrt(discriminant))/(2*a);
    roots[1] = (-b - sqrt(discriminant))/(2*a);
  } else if (discriminant == 0)
  {
    roots[0] = -b/(2*a);
    roots[1]=roots[0];
  }
  return roots;
}

array<double, 3> get_characteristic_polynomial(const array<array<double, DIM>, DIM>& matrix)
{
  auto right_det = matrix[0][1] * matrix[1][0]; // the rhs of determinant, b*c
  auto constant = matrix[0][0] * matrix[1][1];
  auto lambda_coeff_1 = -matrix[0][0];
  auto lambda_coeff_2 = -matrix[1][1];

  array<double, 3> characteristic_polynomial;
  characteristic_polynomial[0] = 1;
  characteristic_polynomial[1] = lambda_coeff_1 + lambda_coeff_2;
  characteristic_polynomial[2] = constant - right_det;
  return characteristic_polynomial;
}

// prints out determinant for visual verification (should be close to 0)

array<double, DIM> multiply_vector_by_scalar(const array<double, DIM>& vec, const double& scalar)
{
  array<double, DIM> result;
  for (int i=0; i<DIM; ++i)
    result.at(i) = vec.at(i) * scalar;
  return result;
}

bool is_zero_matrix(const array<array<double, DIM>, DIM>& matrix)
{
  int threshold = 4;
  int count = 0;
  for (auto row : matrix)
    for (auto elem : row)
      if (elem == 0)
        return true;
  return false;
}

void verify_eigenvalues(const array<array<double, DIM>, DIM>& matrix, const array<double, DIM>& eigen_values)
{
  cout << "DETERMINANTS: ";
  for (auto lambda : eigen_values)
  {
    array<array<double, DIM>, DIM> to_check;
    to_check[0][0] = matrix[0][0] - lambda;
    to_check[0][1] = matrix[0][1];
    to_check[1][0] = matrix[1][0];
    to_check[1][1] = matrix[1][1] - lambda;

    double determinant = (to_check[0][0] * to_check[1][1]) - (to_check[0][1]*to_check[1][0]);
    cout << determinant << " ";
  }
  cout << endl;
}

// prints matrix for the eigenvectors, and the eigenvalues and eigenvectors
void verify_eigenvectors_verbose(const array<array<double, DIM>, DIM>& matrix, const array<array<double, DIM>, DIM>& eigenvectors, const array<double, DIM>& eigenvalues)
{
  cout << "VERIFYING EIGENVECTORS FOR MATRIX" << endl;
  print_matrix(matrix);
  for (int i=0; i<DIM; ++i)
  {
      auto rhs = multiply_vector_by_scalar(eigenvectors[i], eigenvalues[i]);
      auto lhs = multiply_matrix_by_vector(matrix, eigenvectors[i]);
      cout << "RHS: A times eigenvalue " << eigenvalues[i] << endl;
      for (auto elem : lhs)
        cout << elem << " ";
      cout << endl << "LHS: eigen vector: ";
      for (auto elem : eigenvectors[i])
        cout << elem << " ";
      cout << " times eigenvalue: " << eigenvalues[i] << endl;
      for (auto elem : rhs)
        cout << elem << " ";
      cout << endl << endl;;
  }
  cout << endl;
}

// push back errors of the components of the eigenvectors to later see the max and min errors between the components
void collect_eigenvector_errors(const array<array<double, DIM>, DIM>& matrix, const array<array<double, DIM>, DIM>& eigenvectors, const array<double, DIM>& eigenvalues, shared_ptr<vector<double>>& collection)
{
  for (int i=0; i<DIM; ++i)
  {
      auto rhs = multiply_vector_by_scalar(eigenvectors.at(i), eigenvalues.at(i));
      auto lhs = multiply_matrix_by_vector(matrix, eigenvectors.at(i));
      for (int j=0; j<DIM; ++j)
      {
        collection->push_back(lhs.at(j) - rhs.at(j));
      }
  }
}

// visually verify eigenvector accuracy by printing the error between the components
void verify_eigenvectors(const array<array<double, DIM>, DIM>& matrix, const array<array<double, DIM>, DIM>& eigenvectors, const array<double, DIM>& eigenvalues)
{
  cout << "verifying eigenvectors of matrix" << endl;
  print_matrix(matrix);
  for (int i=0; i<DIM; ++i)
  {
      auto rhs = multiply_vector_by_scalar(eigenvectors[i], eigenvalues[i]);
      auto lhs = multiply_matrix_by_vector(matrix, eigenvectors[i]);
      for (int j=0; j<DIM; ++j)
      {
        double error = lhs.at(j) - rhs.at(j);
        cout << "error: " << error << " ";
      }
      cout << endl;
  }
  cout << endl;
}

// iterative method to find eigenvectors
array<double, DIM> jacobi(const array<array<double, DIM>, DIM>& A, const array<double, DIM>& b, array<double, DIM>& x_guess, const double& error_max)
{
  array<double, DIM> x_estimates;
  array<double, DIM> error_vals;
  bool end_loop = false;
  do
  {
    for (int i=0; i<DIM; ++i)
    {
      double rhs_sum = 0;
      for (int j=0; j<DIM; ++j)
        if (j != i)
          rhs_sum += A.at(i).at(j)*x_guess.at(j);
      x_estimates.at(i) = (b.at(i) - rhs_sum) / A.at(i).at(i);
    }

    for (int i=0; i<error_vals.size(); ++i)
      error_vals[i] = fabs(x_guess.at(i) - x_estimates.at(i));

    for (int i=0; i<x_guess.size(); ++i)
      x_guess.at(i) = x_estimates.at(i);

    for (auto val : error_vals)
      if (val < error_max)
      {
        end_loop = true;
        break;
      }
  } while (!end_loop);
  return x_estimates;
}

// method that prepares matrix to send to jacobi iterative method
array<array<double, DIM>, DIM> find_eigenvectors(const array<array<double, DIM>, DIM>& matrix, const array<double, DIM>& eigen_values)
{
  array<array<double, DIM>, DIM> eigen_vectors;
  int count = 0;
  for (auto val : eigen_values)
  {
    array<array<double, DIM>, DIM> temp_matrix;
    for (int i=0; i<DIM; ++i)
      for (int j=0; j<DIM; ++j)
        if (i == j)
          temp_matrix[i][j] = matrix[i][j] - val;
        else
          temp_matrix[i][j] = matrix[i][j];

    array<double, DIM> guess;
    for (int i=0; i<DIM; ++i)
      guess[i] = 1;
    array<double, DIM> b;
    for (int i=0; i < DIM; ++i)
      b[i] = 0;
    eigen_vectors[count] = jacobi(temp_matrix, b, guess, 6);
    count++;
  }
  return eigen_vectors;
}

// equation 30-32
// returns the square root of the eigenvalues of the symmetric matrix R of the polar decomposition of the gradient
array<double, DIM> get_sym_eigenvalues(const array<array<double, DIM>, DIM>& matrix)
{
  // equation 31
  auto poly = get_characteristic_polynomial(matrix);
  auto roots = get_roots(poly);

  // visual verification where the determinants should be close to 0
  // verify_eigen_values(A_t_A, roots);

  array<double, DIM> eigen_diag;
  // equation 32
  for (int i=0; i<DIM; ++i)
    eigen_diag[i] = sqrt(roots[i]);

  return eigen_diag;
}

const char* PWrite(const ParticleSet* a_array)
{
  static int fileCount = 0;
  static char nameBuffer[10];
  if(a_array == NULL)
    {
      return nameBuffer;
    }
  sprintf(nameBuffer, "PART.%d",fileCount);
  PWrite(nameBuffer, a_array);
  fileCount++;
  return nameBuffer;
}

void PWrite(const char* a_filename, const ParticleSet* a_p)
{
  if(a_filename == NULL || a_p == NULL)
    {
      return;
    }
  vector<vector<double> > vars(5);
  unsigned int size = a_p->m_particles.size();
  std::vector<double> x(3*size);
  vars[0] = std::vector<double>(size);
  vars[1] = std::vector<double>(size);
  vars[2] = std::vector<double>(size);
  vars[3] = std::vector<double>(size);
  vars[4] = std::vector<double>(size);
  auto errors = make_shared<vector<double>>();
  for(unsigned int i=0; i<size; i++)
    {
      const Particle& p = a_p->m_particles[i];
      // equation 30
      auto A_t_A = multiply_matrices(get_transpose(p.m_gradx), p.m_gradx); // the symmetric and positive definite matrix
      //equation 30-32
      auto eigenvalues = get_sym_eigenvalues(A_t_A); // get the eigenvalues from the symmetric part of the polar decomp
      double eigen_product = eigenvalues[0] * eigenvalues[1];
      cout << eigen_product << endl;
      array<array<double, DIM>, DIM> eigenvectors;
      bool eigen;
      // avoid dividing by zero by checking for zeros in matrix, might be a bug but one or two zero matrices pop up
      if (!is_zero_matrix(p.m_gradx))
      {
        eigenvectors = find_eigenvectors(p.m_gradx, eigenvalues);
        eigen = true;
        collect_eigenvector_errors(p.m_gradx, eigenvectors, eigenvalues, errors);
      } else
      {
        // do not find eigenvectors of an empty matrix
        eigen = false;
      }

      double max_eigenvalue = -INFINITY;
      // find greatest eigenvalue from the decomp
      for (auto elem : eigenvalues)
        if (elem > max_eigenvalue)
          max_eigenvalue = elem;

      vars[0][i] = p.strength;
      vars[1][i] = p.m_alpha[0];
      vars[2][i] = p.m_alpha[1];
      vars[3][i] = max_eigenvalue;
      vars[4][i] = eigen_product;
      x[i*3] = p.m_x[0];
      x[i*3+1] = p.m_x[1];
#if DIM==3
      x[i*3+2] = p.m_x[2];
#else
      x[i*3+2] = 0.0;
#endif
    }
  double* varPtr[5];
  varPtr[0] = &vars[0][0];
  varPtr[1] = &vars[1][0];
  varPtr[2] = &vars[2][0];
  varPtr[3] = &vars[3][0];
  varPtr[4] = &vars[4][0];
  int vardim[5] = {1,1,1,1,1};
  const char* const varnames[] = {"strength","alpha1","alpha2","max_eigenvalue", "eigen_product"};

  write_point_mesh(a_filename, size,
		   &(x[0]), 5, vardim,
                   varnames, varPtr);
}
inline void write_point_mesh(const char* filename, int npts, double *pts,
                      int nvars, int *vardim, const char * const *varnames,
                      double **vars)
{
    FILE* fp = vtk_open_file(filename);
    int numInColumn = 0;

    int i;
    char  str[128];
    int  *centering = NULL;

    write_header(fp);

    write_string("DATASET UNSTRUCTURED_GRID\n",fp);
    sprintf(str, "POINTS %d double\n", npts);
    write_string(str,fp);
    for (i = 0 ; i < 3*npts ; i++)
    {
        write_double(pts[i],fp,numInColumn);
    }

    new_section(fp,numInColumn);
    sprintf(str, "CELLS %d %d\n", npts, 2*npts);
    write_string(str,fp);
    for (i = 0 ; i < npts ; i++)
    {
        write_int(1,fp,numInColumn);
        write_int(i,fp,numInColumn);
        end_line(fp,numInColumn);
    }

    centering = (int *) malloc(nvars*sizeof(int));
    for (i = 0 ; i < nvars ; i++)
        centering[i] = 1;
    write_variables(nvars, vardim, centering, varnames, vars, npts, npts, fp, numInColumn);
    free(centering);

    vtk_close_file(fp,numInColumn);
}

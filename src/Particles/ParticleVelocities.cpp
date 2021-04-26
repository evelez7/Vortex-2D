#include <iostream>
#include <cassert>
#include <cmath>
#include <memory>
#include <vector>
#include <fstream>
#include "Hockney.H"
#include "ParticleSet.H"
#include "ParticleVelocities.H"
#include "Proto_Point.H"
#include "Proto_Box.H"
#include "Proto_BoxData.H"
#include "Proto_WriteBoxData.H"
//#include "VisitWriter.H"
#include "Proto_Timer.H"
#include "interpolate.h"

using namespace std;
using namespace Proto;
double second_diff(const int &, const double &, Point, const BoxData<double> &);
double second_diff_xy(const double &, Point, const BoxData<double> &);
double G_deriv(const array<int, DIM>&, const Point&, const double&, const BoxData<double>&);
double G_deriv_original(const array<int, DIM>&, const Point&, const double&, const BoxData<double>&);
void print_matrix_here(const array<array<double, DIM>, DIM>&);
void write_curve_file(const BoxData<double> [DIM][DIM], const int&, const double&, const double&, const double&, const string&);
void write_curve_file(shared_ptr<array<vector<double>, 4>>, const int&, const double&, const double&, const double&, const string&);
void add_to_components(shared_ptr<array<vector<double>, 4>>, const array<array<double, DIM>, DIM>&);

bool add_vorticity = false;

double f_smooth_example(double z)
{
  return (12 * pow(z, 10) - 75*pow(z, 8) + 200*pow(z, 6) - 300 * pow(z, 4) + 300*pow(z, 2) - 137)/1200;
}
void GWrite(BoxData<double> G[DIM][DIM], int fileCount)
{
  string g1 = string("G1_finite");
  string g2 = string("G2_finite");
  string g3 = string("G3_finite");
  string g4 = string("G4_finite");
  WriteData(G[0][0], fileCount, 1, g1, g1);
  WriteData(G[0][1], fileCount, 1, g2, g2);
  WriteData(G[1][0], fileCount, 1, g3, g3);
  WriteData(G[1][1], fileCount, 1, g4, g4);
  // END WRITE ALL INDICES AS SEPARATE BOXES`
}

ParticleVelocities:: ParticleVelocities(){};
void ParticleVelocities::operator()
  (ParticleShift& a_k,const double& a_time,const double& a_dt,
   ParticleSet& a_state)
{
  PR_TIMERS("ParticleVelocities");
  PR_TIMER("getRHS",t1);
  PR_TIMER("Hockney",t2);
  PR_TIMER("Fields",t3);
  BoxData<double> rhs(a_state.m_box); // Psi, the field
  Point ipos;
  double xpos[DIM];
  double weight;
  Point e0 = Point::Basis(0);
  Point e1 = Point::Basis(1);
  rhs.setVal(0.);
  PR_START(t1);
  int fileCount = 0;

  int N = a_state.m_box.high()[1];

  vector<Particle >& oldPart = a_state.m_particles;
  vector<DX >& dPart = a_k.m_particles;
  // loop evaluates Psi, the velocity field (rhs)
  for (int k = 0; k < a_state.m_particles.size(); k++)
	{
	  for (int l = 0; l < DIM; l++)
		{
		  double newpos = oldPart[k].m_x[l] + dPart[k].m_x[l];
		  ipos[l] = newpos/a_state.m_dx;
		  xpos[l] = (newpos - ipos[l]*a_state.m_dx)/a_state.m_dx;
		}
	  //Point pt(ipos);
	  int kind = ipos[0] + ipos[1]*(N + 1);
	  for (int l0=0; l0 < DIM;l0++)
		{
		  for (int l1=0;l1 < DIM ; l1++)
			{
			  rhs(ipos + e0*l0 + e1*l1) += (1.-xpos[0] + (2*xpos[0] - 1.)*l0)*
				(1.- xpos[1] + (2*xpos[1] - 1.)*l1)*oldPart[k].strength;
			}
		}
	}
  string filename = "rhs";
  WriteData(rhs, 0, a_state.m_dx, filename, filename);

  // Anayltic solutions for velocity field
  #if EXAMPLE_CASE == 1
  double r0 = 0.5;
  for (auto it = rhs.box().begin(); it != rhs.box().end(); ++it)
  {
    auto current_point = *it;
    double r_sq = pow((current_point[0] * a_state.m_dx) - 0.5, 2.) + pow((current_point[1] * a_state.m_dx) - 0.5, 2.);
    if (r_sq < r0*r0)
    {
      rhs(current_point) = r_sq/2.;
    }
    else if (r_sq > r0*r0)
    {
      rhs(current_point) = pow(r0, 2.) * log(r_sq/(r0*r0)) * 0.5;
    }
  }
  #elif EXAMPLE_CASE == 2
    double r0=0.25;
    for (auto it = rhs.box().begin(); it != rhs.box().end(); ++it)
    {
      auto current_point = *it;

      double r = sqrt(pow(current_point[0] * a_state.m_dx  - 0.5, 2) + pow(current_point[1]*a_state.m_dx - 0.5, 2));
      if (r < r0)
      {
        rhs(current_point) = 5*pow(r0, 2)*f_smooth_example(r/r0);
      } else if (r > r0)
      {
        rhs(current_point) = (pow(r0, 2)/2) * log(r/r0);
      }
    }
  #endif
  string filename_rhs = "rhs";
  WriteData(rhs, 0, a_state.m_dx, filename_rhs, filename_rhs);

  PR_STOP(t1);
  PR_START(t2);
  // a_state.m_hockney.convolve(rhs);
  PR_STOP(t2);
  PR_START(t3);
  Box gradBox = a_state.m_box.grow(-2);
  BoxData<double > field[DIM]; // velocity field on the grid
  for (int l=0 ; l < DIM; l++)
	{
	  field[l].define(gradBox);
	  Point evec = Point::Basis(l);
	  for (auto it = gradBox.begin();!it.done();++it)
		{
		  Point pt = *it;
		  field[l](pt) = (rhs(pt + evec) - rhs(pt - evec))/(2*a_state.m_dx);
		}
	}

  // Equation 63 Calculate G_i matrix
  BoxData<double> G_i_data[DIM][DIM];
  for (int i = 0; i < DIM; ++i)
  {
    for (int j = 0; j < DIM; ++j)
    {
      G_i_data[i][j].define(gradBox);
      for (auto it = G_i_data[i][j].box().begin(); it != G_i_data[i][j].box().end(); ++it)
      {
        auto current_point = *it;
        array<int, DIM> index = {i, j};
        // Uncomment lines here to switch schemes to compute the finite difference matrix G
        double val = G_deriv_original(index, current_point, a_state.m_dx, rhs);
        // double val = G_deriv(index, current_point, a_state.m_dx, rhs);
        G_i_data[i][j](current_point) = val;
      }
    }
  }
  Point low_corner = gradBox.low();
  Point high_corner = gradBox.high();
  Point eval_point = (high_corner * 7)/16;

  // analytic solution @ grid points - interpolated solution
  // for (auto it = rhs.box().begin(); it != rhs.box().end(); ++it)
  // {
  // }

  if (a_state.RK_count == 0)
  {
    string finite_G_identifier = "finite";
    write_curve_file(G_i_data, a_state.m_particles.size(), a_time, a_state.m_hp, a_state.m_dx, finite_G_identifier);
    GWrite(G_i_data, a_state.file_count);
    a_state.file_count++;
  }

  array<array<double, DIM>, DIM> G_k;
  auto G_interpolation_components= make_shared<array<vector<double>, 4>>();

  for (int k = 0; k < a_state.m_particles.size(); ++k)
  {
    // equation 70, pass values to interpolation
    auto omega_k = a_state.m_particles[k].strength * pow(a_state.m_dx / a_state.m_hp, 2.0);
    G_k = interpolate_array(G_i_data, a_state.m_particles[k], a_state.m_dx);

    string interpolation_G_identifier = "interpolation";
    add_to_components(G_interpolation_components, G_k);
    if (add_vorticity)
    {
      // equation 70, add omega_k
      G_k[0][1] += 0.5 * omega_k;
      G_k[1][0] += 0.5 *(-omega_k);
    }

    array<array<double, DIM>, DIM> combined;
    for (int i=0;i<DIM;++i)
      for (int j=0; j<DIM; ++j)
        combined[i][j] = a_state.m_particles[k].m_gradx[i][j] + a_k.m_particles[k].m_gradx[i][j];

    // Equation 62, right hand side, evolve gradient
    for (int i = 0; i < DIM; ++i)
    {
      for (int j = 0; j < DIM; ++j)
      {
        // extra loop for matrix multiplication
        dPart[k].m_gradx[i][j] = 0;
        for (int z = 0; z < DIM; ++z)
          dPart[k].m_gradx[i][j] += (a_dt * G_k[i][z]) * combined[z][j];
      }
    }
    a_state.m_particles[k].G = G_k;
  }
  // write_curve_file(G_interpolation_components, a_state.m_particles.size(), a_time, a_state.m_hp, a_state.m_dx, interpolation_G_identifier);
  // end equation 62, rhs

  for (int k = 0; k < a_state.m_particles.size(); k++)
	{

	  for (int l = 0; l < DIM; l++)
		{
		  double newpos = a_state.m_particles[k].m_x[l] + a_k.m_particles[k].m_x[l];
		  ipos[l] = newpos/a_state.m_dx;
		  xpos[l] = (newpos - ipos[l]*a_state.m_dx)/a_state.m_dx;
		}
	  Point pt(ipos);
	  assert(a_state.m_box.contains(pt));
	  double vel[DIM];
	  for (int dir = 0 ; dir < DIM ; dir++)
		{
		  int dirperp = (dir + 1)%DIM;
		  int sdir = 1-2*dir;
		  double field00 = field[dirperp](pt);
		  double field10 = field[dirperp](pt + e0);
		  double field01 = field[dirperp](pt + e1);
		  double field11 = field[dirperp](pt + e0 + e1);

		  vel[dir] =
			(field00*(1.-xpos[0])*(1. - xpos[1]) +
			 field10*(xpos[0])*(1.-xpos[1]) +
			 field01*(xpos[1])*(1.-xpos[0]) +
			 field11*xpos[0]*xpos[1])*sdir;
		  dPart[k].m_x[dir] = vel[dir]*a_dt;
		}
	  double velmag = sqrt(vel[0]*vel[0] + vel[1]*vel[1]);
	  // cout << "velmag at " << k << " = " << velmag <<endl;
	  //cout << "vel at particle " << k << " = " << vel[0] << " , " << vel[1] << endl;
	}
  PR_STOP(t3);
};

double G_deriv(const array<int, DIM>& index, const Point& current_point, const double& m_dx, const BoxData<double>& rhs)
{
  // determine which derivative to calculate based on the current iterations
  double val;
  int i = index[0]; int j = index[1];
  add_vorticity = false;
  if (i == 0 && j == 0)
  {
    val = -second_diff_xy(m_dx, current_point, rhs);
  }
  else if ((i == 0 && j == 1) || (i==1 && j==0))
  {
    val = 0.5 * (second_diff(0, m_dx, current_point, rhs) - second_diff(1, m_dx, current_point, rhs));
  }
  else if (i == 1 && j == 1)
  {
    val = second_diff_xy(m_dx, current_point, rhs);
  }

  return val;
}

double G_deriv_original(const array<int, DIM>& index, const Point& current_point,  const double& m_dx, const BoxData<double>& rhs)
{
  int i = index[0]; int j = index[1];
  double val;
  add_vorticity = false;
  if (i == 0 && j == 0)
  {
    val = -second_diff_xy(m_dx, current_point, rhs);
  }
  else if (i == 0 && j == 1)
  {
    val=-second_diff(1, m_dx, current_point, rhs);
  }
  else if (i == 1 && j == 0)
  {
    val = second_diff(0, m_dx, current_point, rhs);
  }
  else if (i == 1 && j == 1)
  {
    val = second_diff_xy(m_dx, current_point, rhs);
  }
  return val;
}

/**
 * The finite difference for a second derivative with a respect to a single variable
 *
 * \param axis the variable (x or y) of the derivative
 */
double second_diff(const int &axis, const double &dx, const Point i, const BoxData<double> &function_data) {
  double sum;
	// axis 0 = x, axis 1 = y
	if (axis == 0) {
		// the third element is just i
		Point first(i[0] + 1, i[1]);
		Point third(i[0] - 1, i[1]);
		sum = function_data(first) + (-2.0 * function_data(i)) + function_data(third);
	} else if (axis == 1) {
		Point first(i[0], i[1] + 1);
		Point third(i[0], i[1] - 1);
		sum = function_data(first) + (-2.0 * function_data(i)) + function_data(third);
	} else {
		throw runtime_error("axis not in 2D");
	}
  return (sum / pow(dx, 2.0));
}

/**
 * The finite difference of a partial derivative with respect to x then y
 */
double second_diff_xy(const double &dx, Point i, const BoxData<double> &function_data) {
	Point first(i[0] + 1, i[1] +1);
  Point second(i[0] + 1, i[1] - 1);
  Point third(i[0] - 1, i[1] + 1);
  Point fourth(i[0] - 1, i[1] - 1);
  double sum = function_data(first) - function_data(second) - function_data(third) + function_data(fourth);
  return sum / (4. * pow(dx, 2.));
}

void print_matrix_here(const array<array<double, DIM>, DIM>& matrix)
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

void write_curve_file(const BoxData<double> G_data[DIM][DIM], const int& num_of_particles, const double& time, const double& hp, const double& h, const string& G_type)
{
  int count = 0;
  array<double, 4> maxes;
  array<double, 4> mins;
  for (int i=0; i<DIM; ++i)
  {
    for (int j=0; j<DIM; ++j)
    {
      mins[count] = G_data[i][j].min();
      maxes[count] = G_data[i][j].max();
      count++;
    }
  }
  for (int i=0; i<4; ++i)
  {
    string G_ident;
    switch(i)
    {
      case 0:
        G_ident = "G00_";
        break;
      case 1:
        G_ident = "G01_";
        break;
      case 2:
        G_ident = "G10_";
        break;
      case 3:
        G_ident = "G11_";
        break;
      default:
        break;
    }
    string max_file_name = "max_";
    string min_file_name = "min_";
    max_file_name.append(G_type);
    min_file_name.append(G_type);
    max_file_name.append(G_ident);
    min_file_name.append(G_ident);
    max_file_name.append(to_string(num_of_particles));
    min_file_name.append(to_string(num_of_particles));
    max_file_name.append("_");
    min_file_name.append("_");
    max_file_name.append(to_string(hp));
    min_file_name.append(to_string(hp));
    max_file_name.append("_");
    min_file_name.append("_");
    max_file_name.append(to_string(h));
    min_file_name.append(to_string(h));
    max_file_name.append(".curve");
    min_file_name.append(".curve");
    ofstream max_curve_file(max_file_name, std::ios_base::app);
    ofstream min_curve_file(min_file_name, std::ios_base::app);
    max_curve_file << time << " " << maxes[i] << endl;
    min_curve_file << time << " " << mins[i] << endl;
  }
}

void write_curve_file(shared_ptr<array<vector<double>, 4>> G_components, const int& num_of_particles, const double& time, const double& hp, const double& h, const string& G_type)
{
  int count = 0;
  array<double, 4> maxes;
  array<double, 4> mins;

  for (int i=0; i < G_components->size(); ++i)
  {
    maxes[i] = *max_element(G_components->at(i).begin(), G_components->at(i).end());
    mins[i] = *min_element(G_components->at(i).begin(), G_components->at(i).end());
  }

  for (int i=0; i<4; ++i)
  {
    string G_ident;
    switch(i)
    {
      case 0:
        G_ident = "G00_";
        break;
      case 1:
        G_ident = "G01_";
        break;
      case 2:
        G_ident = "G10_";
        break;
      case 3:
        G_ident = "G11_";
        break;
      default:
        break;
    }
    string max_file_name = "max_";
    string min_file_name = "min_";
    max_file_name.append(G_type);
    min_file_name.append(G_type);
    max_file_name.append(G_ident);
    min_file_name.append(G_ident);
    max_file_name.append(to_string(num_of_particles));
    min_file_name.append(to_string(num_of_particles));
    max_file_name.append("_");
    min_file_name.append("_");
    max_file_name.append(to_string(hp));
    min_file_name.append(to_string(hp));
    max_file_name.append("_");
    min_file_name.append("_");
    max_file_name.append(to_string(h));
    min_file_name.append(to_string(h));
    max_file_name.append(".curve");
    min_file_name.append(".curve");
    ofstream max_curve_file(max_file_name, std::ios_base::app);
    ofstream min_curve_file(min_file_name, std::ios_base::app);
    max_curve_file << time << " " << maxes[i] << endl;
    min_curve_file << time << " " << mins[i] << endl;
  }
}

void add_to_components(shared_ptr<array<vector<double>, 4>> G_k_components, const array<array<double, DIM>, DIM>& G_k)
{
  int count = 0;
  for (int i=0; i<DIM; ++i)
  {
    for (int j=0; j<DIM; ++j)
    {
      G_k_components->at(count).push_back(G_k[i][j]);
    }
  }
}
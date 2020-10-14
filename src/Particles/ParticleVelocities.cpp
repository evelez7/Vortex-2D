#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
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

  int N = a_state.m_box.high()[1];

  const vector<Particle >& oldPart = a_state.m_particles;
  vector<DX >& dPart = a_k.m_particles;
  // loop evaluates Psi, the field (rhs)
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

  PR_STOP(t1);
  PR_START(t2);
  a_state.m_hockney.convolve(rhs);
  PR_STOP(t2);
  PR_START(t3);
  Box gradBox = a_state.m_box.grow(-1);
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
        double val = G_deriv(index, current_point, a_state.m_dx, rhs);
        G_i_data[i][j](current_point) = val;
      }
    }
  }

  for (int k = 0; k < a_state.m_particles.size(); ++k)
  {
    // equation 70, pass values to interpolation
    auto omega_k = a_state.m_particles[k].strength * pow(a_state.m_dx / a_state.m_hp, 2.0);
    auto G_k = interpolate_array(G_i_data, a_state.m_particles[k], a_state.m_dx, omega_k);
    // equation 70, add omega_k
    G_k[0][1] += 0.5 * omega_k;
    G_k[1][0] += 0.5 *(-omega_k);

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
          // dPart[k].m_gradx[i][j] += G_k[i][z] * a_state.m_particles[k].m_gradx[z][j];
          dPart[k].m_gradx[i][j] += (a_dt * G_k[i][z]) * combined[z][j];
        // dPart[k].m_gradx[i][j] *= a_dt;
      }
    }
    print_matrix_here(dPart[k].m_gradx);
  }
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
	Point first(i[0] + 1, i[0] +1);
  Point second(i[0] + 1, i[1] - 1);
  Point third(i[0] - 1, i[1] + 1);
  Point fourth(i[0] - 1, i[1] - 1);
  double sum = function_data(first) - function_data(second) - function_data(third) + function_data(fourth);
  return sum / (4 * pow(dx, 2.0));
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
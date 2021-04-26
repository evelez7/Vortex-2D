#ifndef INTERPOLATE_H
#define INTERPOLATE_H
#include "Particle.H"
#include "Proto_BoxData.H"
#include "Proto_Box.H"
#include "Proto_Point.H"
#include <array>

std::array<std::array<double, DIM>, DIM> interpolate(const Proto::BoxData<double> [DIM][DIM], const Particle&, const double&);
std::array<std::array<double, DIM>, DIM> interpolate_array(const Proto::BoxData<double> [DIM][DIM], const Particle&, const double&);

#endif
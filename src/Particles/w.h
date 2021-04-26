#ifndef W_H
#define W_H

#include <string>
#include <memory>
#include <vector>
#include <array>
#include "Proto_Point.H"

using namespace Proto;

double W_2(double);

double W_3(double);

double W_4(double);

double W_6(double);

double W(const Point&, const double&, double (*)(double));
double W(const std::array<double, DIM>&, const double&, double (*)(double));
std::string get_W_name(double (*)(double));
#endif
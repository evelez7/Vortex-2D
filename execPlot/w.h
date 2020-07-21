#ifndef W_H
#define W_H

#include <string>
#include <memory>
#include <vector>

double W_2(double);

double W_3(double);

double W_4(double);

double W_6(double);

double W(std::shared_ptr<std::vector<double> >, double, double (*)(double));

std::string get_W_name(double (*)(double));
#endif
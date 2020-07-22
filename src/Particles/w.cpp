#include <cmath>
#include <iostream>
#include <map>
#include "w.h"
#include "Proto_Point.H"

// forward declarations
double w_2_helper(double);
double w_3_helper(double);
double w_4_helper(double);
double w_6_helper(double);
// end forward declarations

double W_2(double x) {
    auto abs_x = abs(x);
    if (abs_x >= 0 && abs_x <= 1) {

        return 1-abs_x;
    }    

    return 0;
}

double W_3(double x) {
    auto abs_x = abs(x);

    if (abs_x >= 0 && abs_x <= 1) {
        return 1.0 - ( (5.0/2.0) * pow(abs_x, 2.0) ) + ( (3.0/2.0) * pow(abs_x, 3.0));
    } else if (abs_x >= 1 && abs_x <= 2) {
        return ( (1.0/2.0) * pow(2.0 - abs_x, 2.0) ) * (1.0 - abs_x);
    }

    return 0;
}

double W_4(double x) {
    auto abs_x = abs(x);

    if (abs_x >= 0.0 && abs_x <= 1.0) {
        return 1.0 - (abs_x / 2.0) - pow(abs_x, 2.0) + (pow(abs_x, 3.0) / 2.0); 
    } else if (abs_x >= 1.0 && abs_x <= 2.0) {
        return 1.0 - ( (11.0 * abs_x) / 6.0) + pow(abs_x, 2.0) - (pow(abs_x, 3.0) / 6.0);
    } 

    return 0;
}

double W_6(double x) {
    auto abs_x = abs(x);

    if (abs_x >= 0.0 && abs_x <= 1.0) {
        return 1.0 - (abs_x / 3.0) - ( ( 5.0*pow(abs_x, 2.0) ) / 4.0) + ( (5.0 * pow(abs_x, 3.0)) / 12.0) + (pow(abs_x, 4.0) / 4.0) - (pow(abs_x, 5.0) / 12.0);
    } else if (abs_x >= 1.0 && abs_x <= 2.0) {
        return 1.0 - ( (13.0 * abs_x) / 12.0) - ( (5.0 * pow(abs_x, 2.0)) / 8.0) + ( (25.0 * pow(abs_x, 3.0)) / 24.0) - ( (3.0 * pow(abs_x, 4.0)) / 8.0) + (pow(abs_x, 5.0) / 24.0);
    } else if (abs_x >= 2.0 && abs_x <= 3.0) {
        return 1.0 - ( (137.0 * abs_x) / 60.0 ) + ( (15.0 * pow(abs_x, 2.0)) / 8.0 ) - ( (17.0 * pow(abs_x, 3.0)) / 24.0 ) + (pow(abs_x, 4.0) / 8.0) - (pow(abs_x, 5.0) / 120.0); 
    }

    return 0;
}

double W(const Point &z, const double &h_g, double (*w_ptr)(double)) {
    double product = 1.0;
    for (int i = 0; i < DIM; ++i) {
        product*=w_ptr(z[i] / h_g);
    }
    return product / (pow(h_g, DIM));
}

std::string get_W_name(double (*w_script)(double)) {
    if (w_script == W_2) {
        return "W2";
    } else if (w_script == W_3) {
        return "W3";
    } else if (w_script == W_4) {
        return "W4";
    } else if (w_script == W_6) {
    }
    return NULL;
}
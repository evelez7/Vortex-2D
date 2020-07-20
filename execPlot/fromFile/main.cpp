#include "CutoffKernel.H"
#include "ParticleSet.H"
#include "ParticleVelocities.H"
#include "Proto_Box.H"
#include "Proto_BoxData.H"
#include "Proto_Point.H"
#include "Proto_WriteBoxData.H"
#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
//#include "VisitWriter.H"
#include "ParticleWriter.H"
#include "Proto_RK4.H"
using namespace std;
namespace fs = std::filesystem;

int main(int argc, char *argv[]) {
    fs::path current_path = fs::current_path();
    fs::path input_path = current_path / "input"; // the directory to files must be called input

    if (!fs::is_directory(input_path)) {
        throw fs::filesystem_error("Programs path does not exist", input_path, std::error_code());
    }

    for (auto &dir_entry : fs::directory_iterator(input_path)) {
        if (dir_entry.path().filename() == ".DS_Store") {
            printf("Skipping .DS_Store\n");
            continue;
        }
        cout << "Opening " << dir_entry.path().filename() << endl;
        ifstream input_file(dir_entry.path().c_str(), ios::binary);
        cout << "Parsing " << dir_entry.path().filename() << endl;

        string line;
        vector<vector<double>> interpolation_values;
        int rows = 0;
        while (getline(input_file, line)) {
            interpolation_values.push_back(vector<double>()); // a new line means a new row
            istringstream tokenized_line(line);
            string token;
            while (tokenized_line >> token) {
                if (tokenized_line.peek() == ' ') {
                    tokenized_line.ignore();
                }
                double interpolation_value = stod(token);
                interpolation_values.at(rows).push_back(interpolation_value);
            }
            rows++;
        }

        if (interpolation_values.size() == 0) {
            cout << "Skipping empty file" << endl << endl;
            continue;
        }

        auto N = interpolation_values.size(); // the rows of the matrix are determined by the iterations of the interpolation
        auto interpolation_box_data = BoxData<double>(Box::Cube(N));
        for (int i = 0; i <= N - 1; ++i) {
            for (int j = 0; j <= N - 1; ++j) {
                interpolation_box_data(Point(i, j)) = interpolation_values.at(i).at(j);
            }
        }
        string file_name = "inteprolation_data";
        string var_name = "./output/" + dir_entry.path().filename().string();
        double h_g = 2.0*(1.0 /(2.0 * static_cast<double>(N)));
        WriteData(interpolation_box_data, 0, h_g, file_name, var_name);
        cout << endl;
    }
    return 0;
}
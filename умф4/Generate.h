#pragma once
#include <vector>
#include <fstream>

void Create_grid_from_points(std::string path, double lambda = 1, double sigma = 1, double hi = 1);
void Create_material_file(std::string path, double lambda = 1, double sigma = 1, double hi = 1);
void Create_time_grid();
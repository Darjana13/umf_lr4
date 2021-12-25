#pragma once
#include <stdio.h>
#include <fstream>
#include <vector>
#include <string>
#include <set>
#include "LU_solver.h"
#include "Solver.h"

using namespace std;

struct node
{
	double x, y, z;
};

struct elem
{
	vector<int> nodes;
	int material;
};

struct material
{
	double lambda, sigma, hi;
};

class FEM
{
	vector<node> all_nodes; 
	vector<elem> all_elems;
	vector<material> all_materials;
	vector<vector<int>> S1;

	vector<vector<double>> M1, G1;
	vector<double> b;
	vector<int> ia, ja;
	vector<double> al, au, di;

	double omega = 1;

	double func_s(double x, double y, double z);
	double func_c(double x, double y, double z);
	double func_S_s(double x, double y, double z, int s_id);
	double func_S_c(double x, double y, double z, int s_id);
	double func_true_s(double x, double y, double z, int s_id);
	double func_true_c(double x, double y, double z, int s_id);

	void Set_const_matrix();
	int Get_M_Loc(double hx, double hy, double hz, double gam, vector<vector<double>>& M_loc);
	int Get_G_Loc(double hx, double hy, double hz, double lam, vector<vector<double>>& G_loc);
	int Get_Loc(vector<vector<double>>& p_loc, vector<vector<double>>& c_loc, vector<double>& fc_loc, vector<double>& fs_loc, int id);
	
	int GeneratePortrait();
	int AddLocal(int el_id, vector<vector<double>>& p_loc, vector<vector<double>>& c_loc, vector<double>& fc_loc, vector<double>& fs_loc);
	
	int Create_grid(string file);
	int SetS1();
public:
	vector<double> q0;
	int Make_SLAU(string file, double om = 1);
	void Solve_SLAU(bool LU);
	double Calc_error();
};
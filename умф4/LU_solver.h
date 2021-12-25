#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
using namespace std;
const double eps = 1E-20;

class LU_solver
{
	int N = 0; // Размер матрицы
	vector<double> di; // Диагональ
	vector<int> ia; // Портрет матрицы
	vector<double> au; // верхний треугольник
	vector<double> al; // нижний треугольник
	vector<double> b; // Вектор

	bool LU();
	void revers();
	void forward_stroke(); // Прямой ход Ly = F

	int input(string path);
	void output();
public:
	LU_solver();
	LU_solver(string path);
	LU_solver(std::vector<int>& _ia, std::vector<int>& _ja, std::vector<double>& _di,
		std::vector<double>& _al, std::vector<double>& _au,
		std::vector<double>& _b);

	int Solve_task();
	void getx0(vector<double> &q);

};
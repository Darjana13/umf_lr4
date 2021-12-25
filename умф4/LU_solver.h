#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
using namespace std;
const double eps = 1E-20;

class LU_solver
{
	int N = 0; // ������ �������
	vector<double> di; // ���������
	vector<int> ia; // ������� �������
	vector<double> au; // ������� �����������
	vector<double> al; // ������ �����������
	vector<double> b; // ������

	bool LU();
	void revers();
	void forward_stroke(); // ������ ��� Ly = F

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
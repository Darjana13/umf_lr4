#include "FEM.h"
#include "Generate.h"
#include <time.h> 

int main()
{
	double om[] = { 0.001, 1, 100, 1E+6, 1E+9 },
		sig[] = { 0, 1000, 1000000 },
		lamb[] = {100, 1000, 10000 },
		hii[] = {8.81E-11, 1E-10, 1E-9};
	int iter = 5*3*3*3;

	cout << "start make file grid\n";
	Create_grid_from_points("");
	cout << "end make file grid\n";

	ofstream fileLU("result_LU.txt"), fileIter("result_LOS.txt"), file("table.txt");
	fileLU.imbue((locale)"Rus");
	fileIter.imbue((locale)"Rus");

	
	
	// one task
	/*FEM task;
	task.Make_SLAU("");
	clock_t start = clock();
	task.Solve_SLAU(true);
	clock_t end = clock();
	file << task.Calc_error() << "\t" << (double)(end - start) / CLOCKS_PER_SEC << "\t";
	for (int i = 0; i < task.q0.size(); i++)
	{
		fileLU << task.q0[i] << "\n";
	}

	start = clock();
	task.Solve_SLAU(false);
	end = clock();
	file << task.Calc_error() << "\t" << (double)(end - start) / CLOCKS_PER_SEC << "\n";
	for (int i = 0; i < task.q0.size(); i++)
	{
		fileIter << task.q0[i] << "\n";
	}*/

	// tasks with different params
	for (int o = 0; o < sizeof(om) / 8; o++)
	{
		for (int s = 0; s < sizeof(sig) / 8; s++)
		{
			for (int l = 0; l < sizeof(lamb) / 8; l++)
			{
				for (int h = 0; h < sizeof(hii) / 8; h++)
				{
					FEM task;
					Create_material_file("", lamb[l], sig[s], hii[h]);
					cout << "start make SLAU\n";
					task.Make_SLAU("", om[o]);
					cout << "end make SLAU\n";

					file << om[o] << "\t" << sig[s] << "\t" << lamb[l] << "\t" << hii[h] << "\t";

					clock_t start = clock();
					//task.Solve_SLAU(true);
					clock_t end = clock();
					//file << task.Calc_error() << "\t" << (double)(end - start) / CLOCKS_PER_SEC << "\t";

					//for (int i = 0; i < task.q0.size(); i++)
					//{
					//	fileLU << task.q0[i] << "\n";
					//}

					start = clock();
					task.Solve_SLAU(false);
					end = clock();
					file << task.Calc_error() << "\t" << (double)(end - start) / CLOCKS_PER_SEC << "\n";
					//for (int i = 0; i < task.q0.size(); i++)
					//{
					//	fileIter << task.q0[i] << "\n";
					//}

					cout << "\n--------------------------------\nNEED PARAM ITER " << iter << "\n";
					iter--;
				}
			}
		}
	}
	
	return 0;
}
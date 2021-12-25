#include "FEM.h"

double FEM::func_s(double x, double y, double z)
{
	// один материал везде
	// div grad u = 0
	return -omega * all_materials[0].sigma* func_true_c(x,y,z,0)- omega * omega * all_materials[0].hi* func_true_s(x, y, z, 0);
}
double FEM::func_c(double x, double y, double z)
{
	// один материал везде
    // div grad u = 0
	return  omega * all_materials[0].sigma * func_true_s(x, y, z, 0) - omega * omega * all_materials[0].hi * func_true_c(x, y, z, 0);
}
double FEM::func_true_c(double x, double y, double z, int id)
{
	switch (id)
	{
	case 0:
		return //3
			2 * x + y - 3 * z
			;
	default:
		break;
	}
	return 1;
}
double FEM::func_true_s(double x, double y, double z, int id)
{
	switch (id)
	{
	case 0:
		return //5
			4 * x + 5 * y - 7 * z
			;
	default:
		break;
	}
	return 1;
}
double FEM::func_S_c(double x, double y, double z, int s_id)
{
	return func_true_c(x, y, z, s_id);
}
double FEM::func_S_s(double x, double y, double z, int s_id)
{
	return func_true_s(x, y, z, s_id);
}

int FEM::Create_grid(string file)
{
	int N, Nmat, Kel, NS1, NS;
	std::ifstream in;

	in.open("info.txt");
	in >> N >> Nmat >> Kel >> NS1;
	in.close();

	in.open("xyz.txt");
	all_nodes.resize(N);
	for (int i = 0; i < N; i++)
	{
		in >> all_nodes[i].x >> all_nodes[i].y >> all_nodes[i].z;
	}
	in.close();

	q0.resize(2*N);

	in.open("S1.txt");
	S1.resize(NS1);
	for (int i = 0; i < NS1; i++)
	{
		int size = 0;
		in >> size;
		S1[i].resize(size);
		for (int j = 0; j < size; j++)
		{
			in >> S1[i][j];
		}
	}
	in.close();

	

	in.open("material.txt");
	all_materials.resize(Nmat);
	for (int i = 0; i < Nmat; i++)
	{
		in >> all_materials[i].lambda >> all_materials[i].sigma >> all_materials[i].hi;
	}
	in.close();

	in.open("elem.txt");
	all_elems.resize(Kel);
	for (int i = 0; i < Kel; i++)
	{
		all_elems[i].nodes.resize(8);
		for (int j = 0; j < 8; j++)
		{
			in >> all_elems[i].nodes[j];
		} 	
		in >> all_elems[i].material;
	}
	in.close();



	return 0;
}

void FEM::Set_const_matrix()
{
	M1.resize(2);
	G1.resize(2);

	for (int i = 0; i < 2; i++)
	{
		M1[i].resize(2);
		G1[i].resize(2);
	}
	M1[0][0] = 1. / 3;
	M1[0][1] = 1. / 6;
	M1[1][0] = 1. / 6;
	M1[1][1] = 1. / 3;

	G1[0][0] = 1;
	G1[0][1] = -1;
	G1[1][0] = -1;
	G1[1][1] = 1;
}

int FEM::Get_M_Loc(double hx, double hy, double hz, double gam, vector<vector<double>>& M_loc)
{

	for (int i = 0; i < 8; i++)
	{
		int x_i = (i) % 2,
			y_i = ((i) / 2) % 2,
			z_i = ((i) / 4);
		//cout << i << "\t" << x_i << y_i << z_i << endl;
		for (int j = 0; j < 8; j++) //for (int j = 0; j < i; j++) M_loc[j][i] = M_loc[i][j];
		{
			int x_j = (j ) % 2,
				y_j = ((j ) / 2) % 2,
				z_j = ((j ) / 4);
			M_loc[i][j] = /*gam **/ hx * M1[x_i][x_j] * hy * M1[y_i][y_j] * hz * M1[z_i][z_j];
		}
	}
	/*M_loc[0][0] = gam * hx * M1[0][0] * hy * M1[0][0] * hz * M1[0][0];
	M_loc[0][1] = gam * hx * M1[0][1] * hy * M1[0][0] * hz * M1[0][0];
	M_loc[0][2] = gam * hx * M1[0][0] * hy * M1[0][1] * hz * M1[0][0];
	M_loc[0][3] = gam * hx * M1[0][1] * hy * M1[0][1] * hz * M1[0][0];
	M_loc[0][4] = gam * hx * M1[0][0] * hy * M1[0][0] * hz * M1[0][1];
	M_loc[0][5] = gam * hx * M1[0][1] * hy * M1[0][0] * hz * M1[0][1];
	M_loc[0][6] = gam * hx * M1[0][0] * hy * M1[0][1] * hz * M1[0][1];
	M_loc[0][7] = gam * hx * M1[0][1] * hy * M1[0][1] * hz * M1[0][1];*/
	return 0;
}


int FEM::Get_G_Loc(double hx, double hy, double hz, double lam, vector<vector<double>>& G_loc)
{
	for (int i = 0; i < 8; i++)
	{
		int x_i = (i)  % 2,
			y_i = ((i) / 2) % 2,
			z_i = ((i) / 4);
		for (int j = 0; j < 8; j++)
		{
			int x_j = (j) % 2,
				y_j = ((j ) / 2) % 2,
				z_j = ((j ) / 4);
			G_loc[i][j] = lam * (
				hx * G1[x_i][x_j] * hy * M1[y_i][y_j] * hz * M1[z_i][z_j] +
				hy * M1[x_i][x_j] * hy * G1[y_i][y_j] * hz * M1[z_i][z_j] + 
				hz * M1[x_i][x_j] * hy * M1[y_i][y_j] * hz * G1[z_i][z_j]);
		}
	}
	return 0;
}

/*int FEM::Get_b_Loc(double hx, double hy, double hz, vector<double>& b_loc, vector<double>& f_loc)
{
	vector<vector<double>> M_loc(8, vector<double>(8));
	for (int i = 0; i < 8; i++)
	{
		int x_i = (i - 1) % 2,
			y_i = ((i - 1) / 2) % 2,
			z_i = ((i - 1) / 4);
		for (int j = 0; j < 8; j++)
		{
			int x_j = (j - 1) % 2,
				y_j = ((j - 1) / 2) % 2,
				z_j = ((j - 1) / 4);
			M_loc[i][j] = hx * M1[x_i][x_j] * hy * M1[y_i][y_j] * hz * M1[z_i][z_j];
		}
	}
	fill(b_loc.begin(), b_loc.end(), 0);
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 8; j++)
		{
			b_loc[i] += M_loc[i][j] * f_loc[j];
		}
	}
}*/

int FEM::Get_Loc(vector<vector<double>>& p_loc, vector<vector<double>>& c_loc, vector<double>& fc_loc, vector<double>& fs_loc, int id)
{
	vector<vector<double>> M_loc(8, vector<double>(8));
	vector<vector<double>> G_loc(8, vector<double>(8));
	vector<double> f_c(8), f_s(8);
	double 
		x1 = all_nodes[all_elems[id].nodes[0]].x,
		x2 = all_nodes[all_elems[id].nodes[1]].x,
		y1 = all_nodes[all_elems[id].nodes[0]].y,
		y2 = all_nodes[all_elems[id].nodes[3]].y,
		z1 = all_nodes[all_elems[id].nodes[0]].z,
		z2 = all_nodes[all_elems[id].nodes[4]].z;
	double 
		hx = x2 - x1,
		hy = y2 - y1,
		hz = z2 - z1;
	Get_M_Loc(hx, hy, hz, 1, M_loc);
	Get_G_Loc(hx, hy, hz, all_materials[all_elems[id].material].lambda, G_loc);
	double 
		sigma = all_materials[all_elems[id].material].sigma,
		hi = all_materials[all_elems[id].material].hi;
	for (int i = 0; i < 8; i++)
	{
		node temp = all_nodes[all_elems[id].nodes[i]];
		fc_loc[i] = 0;
		fs_loc[i] = 0;
		f_c[i] = func_c(temp.x, temp.y, temp.z);
		f_s[i] = func_s(temp.x, temp.y, temp.z);
	}
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 8; j++)
		{
			p_loc[i][j] = G_loc[i][j] - omega * omega * hi * M_loc[i][j];
			c_loc[i][j] = omega * sigma * M_loc[i][j];
			if (i != j)
			{
				p_loc[j][i] = p_loc[i][j];
				c_loc[j][i] = c_loc[i][j];
			}

			fc_loc[i] += f_c[j] * M_loc[i][j];
			fs_loc[i] += f_s[j] * M_loc[i][j];
		}
	}
	return 0;
}

int FEM::AddLocal(int el_id, vector<vector<double>>& p_loc, vector<vector<double>>& c_loc, vector<double>& fc_loc, vector<double>& fs_loc)
{
	vector<int> global_num(8);
	for (int i = 0; i < 8; i++)
		global_num[i] = 2 * all_elems[el_id].nodes[i];

	for (int i = 0; i < 8; i++) 
	{
		di[global_num[i]] += p_loc[i][i];
		di[global_num[i] + 1] += p_loc[i][i];

		int end0 = ia[global_num[i] + 2] - 1;
		int ind0 = end0;
		al[ind0] += c_loc[i][i]; //-c11
		au[ind0] -= c_loc[i][i];

		// блоки нижнего треугольника

		// верхняя строка блока
		int beg = ia[global_num[i]];
		for (int j = 0; j < i; j++, beg++) 
		{
			int end = ia[global_num[i] + 1] - 1;
			int ind = beg;
			while (ja[ind] != global_num[j]) 
			{
				ind++;
			}
			al[ind] += p_loc[i][j]; //p21
			al[ind + 1] -= c_loc[i][j]; //-c21
		}
		// нижняя строка блока
		int beg2 = ia[global_num[i] + 1];
		for (int j = 0; j < i; j++, beg++) 
		{
			int end2 = ia[global_num[i] + 2] - 1;
			int ind2 = beg2;
			while (ja[ind2] != global_num[j]) 
			{
				ind2++;
			}
			al[ind2] += c_loc[i][j]; //c21
			al[ind2 + 1] += p_loc[i][j]; //p21
		}

		// блоки верхнего треугольника
		// левый столбец
		beg = ia[global_num[i]];
		for (int j = 0; j < i; j++, beg++) 
		{
			int end = ia[global_num[i] + 1] - 1;
			int ind = beg;
			while (ja[ind] != global_num[j]) 
			{
				ind++;
			}
			//au[ind] += p_loc[i][j]; //p21
			//au[ind + 1] = c_loc[i][j]; //-c21
			au[ind] += p_loc[i][j];     // p21
			au[ind + 1] += c_loc[i][j]; // c21
		}
		// правый столбец
		beg2 = ia[global_num[i] + 1];
		for (int j = 0; j < i; j++, beg++) {
			int end2 = ia[global_num[i] + 2] - 1;
			int ind2 = beg2;
			while (ja[ind2] != global_num[j]) {
				ind2++;
			}
			au[ind2] -= c_loc[i][j];     // -c21
			au[ind2 + 1] += p_loc[i][j]; // p21
		}

		b[global_num[i]] += fs_loc[i];
		b[global_num[i] + 1] += fc_loc[i];
	}
	return 0;
}


int FEM::GeneratePortrait()
{
	int m = all_nodes.size();
	vector<set<int>> list(m);
	for (int elem_id = 0; elem_id < all_elems.size(); elem_id++) 
	{
		for (int i = 0; i < 8; i++) 
		{
			for (int j = i + 1; j < 8; j++) 
			{
				int ind1 = all_elems[elem_id].nodes[i];
				int ind2 = all_elems[elem_id].nodes[j];

				if (ind1 < ind2) 
					swap(ind1, ind2);
				list[ind1].insert(ind2);	
			}	
		}	
	}
	ia.resize(2*m + 1);

	//создание портрета по списку
	ia[0] = 0;
	ia[1] = 0;
	ia[2] = 1;

	for (int i = 1; i < m; i++) 
	{
		ia[2 * i + 1] = ia[2 * i] + list[i].size() * 2;
		ia[2 * (i + 1)] = ia[2 * i + 1] + list[i].size() * 2 + 1;
		//if (i % 2 == 1) ig[i]++;
	}
	ja.resize(ia[2 * m]);
	au.resize(ia[2 * m]);
	al.resize(ia[2 * m]);
	di.resize(2 * m);
	b.resize(2 * m);
	/*ggu_new.resize(ig[2 * m]);
	ggl_new.resize(ig[2 * m]);
	di_new.resize(2 * m);*/


	for (int i = 1, k = 1; i < m; i++) 
	{
		for (int j : list[i]) 
		{
			ja[k] = 2 * j;
			ja[k + 1] = 2 * j + 1;
			k += 2;
		}
		for (int j : list[i]) {

			ja[k] = 2 * j;
			ja[k + 1] = 2 * j + 1;
			k += 2;
		}
		ja[k] = 2 * i;
		k++;
	}
	return 0;
}

//int FEM::GeneratePortrait() // генерация портрета
//{
//	int N = all_nodes.size();
//	int Kel = all_elems.size();
//	diM.resize(N);
//	diG.resize(N);
//	ia_temp.resize(N + 1);
//	ja_temp.resize(64 * Kel);
//	std::vector<int> temp_list1(64 * Kel),
//		temp_list2(64 * Kel);
//	std::vector<int> listbeg(N);
//	int listsize = 0;
//	for (int i = 0; i < N; i++)
//	{
//		listbeg[i] = 0;
//	}
//	for (int ielem = 0; ielem < Kel; ielem++)
//	{
//		for (int i = 0; i < 8; i++)
//		{
//			int k = all_elems[ielem].nodes[i];
//			for (int j = i + 1; j < 8; j++)
//			{
//				int ind1 = k;
//				int ind2 = all_elems[ielem].nodes[j];
//				if (ind2 < ind1)
//				{
//					ind1 = ind2;
//					ind2 = k;
//				}
//				int iaddr = listbeg[ind2];
//				if (iaddr == 0)
//				{
//					listsize++;
//					listbeg[ind2] = listsize;
//					temp_list1[listsize] = ind1;
//					temp_list2[listsize] = 0;
//				}
//				else
//				{
//					while (temp_list1[iaddr] < ind1 && temp_list2[iaddr] > 0)
//					{
//						iaddr = temp_list2[iaddr];
//					}
//					if (temp_list1[iaddr] > ind1)
//					{
//						listsize++;
//						temp_list1[listsize] = temp_list1[iaddr];
//						temp_list2[listsize] = temp_list2[iaddr];
//						temp_list1[iaddr] = ind1;
//						temp_list2[iaddr] = listsize;
//					}
//					else if (temp_list1[iaddr] < ind1)
//					{
//						listsize++;
//						temp_list2[iaddr] = listsize;
//						temp_list1[listsize] = ind1;
//						temp_list2[listsize] = 0;
//					}
//				}
//			}
//		}
//	}
//
//	ia_temp[0] = 0;
//	for (int i = 0; i < N; i++)
//	{
//		ia_temp[i + 1] = ia_temp[i];
//		int iaddr = listbeg[i];
//		while (iaddr != 0)
//		{
//			ja_temp[ia_temp[i + 1]] = temp_list1[iaddr];
//			ia_temp[i + 1]++;
//			iaddr = temp_list2[iaddr];
//		}
//	}
//
//	ja_temp.resize(ia_temp[N]);
//	alM.resize(ia_temp[N]);
//	alG.resize(ia_temp[N]);
//
//	return 0;
//}
//
//int FEM::AddLocalM(int el_id, std::vector<vector<double>>& M_loc)
//// внесение локальных A в глобальную СЛАУ
//{
//	std::vector<int> L = all_elems[el_id].nodes;
//	int k = all_elems[el_id].nodes.size(); // размерность локальной матрицы
//	for (int i = 0; i < 8; i++)
//	{
//		diM[L[i]] += M_loc[i][i];
//	}
//
//	for (int i = 0; i < 8; i++)
//	{
//		int temp = ia_temp[L[i]];
//		for (int j = 0; j < i; j++)
//		{
//			for (int k = temp; k < ia_temp[L[i] + 1]; k++)
//			{
//				if (ja_temp[k] == L[j])
//				{
//					alM[k] += M_loc[i][j];
//					k++;
//					break;
//				}
//			}
//		}
//	}
//
//	return 0;
//}
//int FEM::AddLocalG(int el_id, std::vector<vector<double>>& G_loc)
//// внесение локальных A в глобальную СЛАУ
//{
//	std::vector<int> L = all_elems[el_id].nodes;
//	int k = all_elems[el_id].nodes.size(); // размерность локальной матрицы
//	for (int i = 0; i < 8; i++)
//	{
//		diG[L[i]] += G_loc[i][i];
//	}
//
//	for (int i = 0; i < 8; i++)
//	{
//		int temp = ia_temp[L[i]];
//		for (int j = 0; j < i; j++)
//		{
//			for (int k = temp; k < ia_temp[L[i] + 1]; k++)
//			{
//				if (ja_temp[k] == L[j])
//				{
//					alG[k] += G_loc[i][j];
//					k++;
//					break;
//				}
//			}
//		}
//	}
//
//	return 0;
//}


/*int AddLocal_b(std::vector<double>& b, std::vector<double>& b_loc, int el_id)
// внесение локальных b  в глобальную СЛАУ
{
    std::vector<int> L = all_elems[el_id].node_loc;
    int k = all_elems[el_id].node_loc.size(); // размерность локальной матрицы
    for (int i = 0; i < k; i++)
    {
        b[L[i]] += b_loc[i];
    }
    return 0;
}*/

int FEM::SetS1() // учет первых краевых
{
	int NS1 = S1.size();
	for (int i = 0; i < NS1; i++)
	{
		int s1_id = i;
		for (int j = 0; j < S1[i].size(); j++)
		{
			int node_id = 2 * S1[i][j],
				node_glob_num = S1[i][j];
			//cout << "S1 on node " << S1[i][j] << "\n";

			// first
			di[node_id] = 1;
			b[node_id] = func_S_s(all_nodes[node_glob_num].x, all_nodes[node_glob_num].y, all_nodes[node_glob_num].z, s1_id);
			for (int k = ia[node_id]; k < ia[node_id + 1]; k++)
			{
				al[k] = 0;
			}
			for (int k = 0; k < ja.size(); k++)
			{
				if (ja[k] == node_id)
				{
					au[k] = 0;
				}
			}
			// second
			node_id++;
			di[node_id] = 1;
			b[node_id] = func_S_c(all_nodes[node_glob_num].x, all_nodes[node_glob_num].y, all_nodes[node_glob_num].z, s1_id);
			for (int k = ia[node_id]; k < ia[node_id + 1]; k++)
			{
				al[k] = 0;
			}
			for (int k = 0; k < ja.size(); k++)
			{
				if (ja[k] == node_id)
				{
					au[k] = 0;
				}
			}
		}
	}
	return 0;
}

int FEM::Make_SLAU(string file, double om)
{
	omega = om;

	Create_grid("");
	Set_const_matrix();
	GeneratePortrait();
	vector<vector<double>> p_loc(8, vector<double>(8)),
		c_loc(8, vector<double>(8));
	vector<double> fc_loc(8), fs_loc(8);
	for (int i = 0; i < all_elems.size(); i++)
	{
		Get_Loc(p_loc, c_loc, fc_loc, fs_loc, i);
		AddLocal(i, p_loc, c_loc, fc_loc, fs_loc);
		//cout << "on elem " << i << "\n";
	}
	SetS1();

	return 0;
}

void FEM::Solve_SLAU(bool LU)
{
	if (!LU)
	{
		Solver t(ia, ja, di, al, au, b);
		//t.LOS_LU();
		//t.CGM_LU();
		//t.BSG_LU();
		t.BSG();

		t.getx0(q0);
	}
	else
	{
		LU_solver t2(ia, ja, di, al, au, b);
		t2.Solve_task();
		t2.getx0(q0);
	}
}

double FEM::Calc_error()
{
	double res = 0, f_s, f_c, norm = 0;
	for (int i = 0; i < all_nodes.size(); i++)
	{
		f_s = func_true_s(all_nodes[i].x, all_nodes[i].y, all_nodes[i].z, 0);
		f_c = func_true_c(all_nodes[i].x, all_nodes[i].y, all_nodes[i].z, 0);

		res += (f_s - q0[2 * i]) * (f_s - q0[2 * i]);
		res += (f_c - q0[2 * i + 1]) * (f_c - q0[2 * i + 1]);

		norm += f_s * f_s + f_c * f_c;
	}
	return sqrt(res / norm);
}
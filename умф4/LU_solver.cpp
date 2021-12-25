#include "LU_solver.h"

LU_solver::LU_solver()
{
}
LU_solver::LU_solver(string path)
{
    input(path);
}
LU_solver::LU_solver(std::vector<int>& _ia, std::vector<int>& _ja, std::vector<double>& _di,
    std::vector<double>& _al, std::vector<double>& _au,
    std::vector<double>& _b)
{
    N = _di.size();
    di = _di;
    b = _b;
    ia.resize(N + 1);
    al.reserve(N * 1000);
    au.reserve(N * 1000);

    ia[0] = 0;
    for (int i = 1; i < N; i++)
    {
        if (_ia[i] != _ia[i + 1])
        {
            int current = _ja[_ia[i]];
                //end = i - current; // кол-во элементов в профиле
            int k = _ia[i];
            for (; current < i; current++) // по всем элементам строки
            {
                if (k < _ia[i+1] && _ja[k] == current)
                {
                    al.push_back(_al[k]);
                    au.push_back(_au[k]);
                    k++;
                }
                else
                {
                    al.push_back(0);
                    au.push_back(0);
                }
            }
        }
        ia[i + 1] = al.size();

    }
}
int  LU_solver::input(string path)
{
    ifstream fin;
    fin.open(path + "IN.txt");

    fin >> N;
    di.resize(N);
    for (int i = 0; i < N; i++)
        fin >> di[i];

    ia.resize(N + 1);
    for (int i = 0; i <= N; i++)
        fin >> ia[i];

    au.resize(ia[N]);
    for (int i = 0; i < ia[N]; i++)
        fin >> au[i];

    al.resize(ia[N]);
    for (int i = 0; i < ia[N]; i++)
        fin >> al[i];

    b.resize(N);
    for (int i = 0; i < N; i++)
        fin >> b[i];

    fin.close();
    return 0;
}

void LU_solver::output()
{
    ofstream fout;
    fout.open("out.txt");
    fout.precision(2 * sizeof(double) - 1);

    for (int i = 0; i < N; i++)
        fout << b[i] << "\n";
    fout.close();
}
void LU_solver::getx0(vector<double> &q)
{
    if (q.size() != N)
        q.resize(N);
    for (int i = 0; i < N; i++)
        q[i] = b[i];
}
void LU_solver::forward_stroke() // Прямой ход Ly = F
{
    double res = 0;
    for (int i = 0; i < N; i++)
    {
        int count = i - (ia[i + 1] - ia[i]);
        res = 0;
        for (int j = ia[i]; j < ia[i + 1]; j++, count++)
        {
            res += b[count] * al[j];
        }
        b[i] = (b[i] - res) / di[i];
    }
}

void LU_solver::revers()
{

    for (int i = N - 1; i >= 0; i--)
    {
        int j = i - (ia[i + 1] - ia[i]);
        for (int k = ia[i]; k < ia[i + 1]; k++, j++)
        {
            b[j] -= au[k] * b[i];
        }
    }

}

bool LU_solver::LU()
{
    for (int i = 1; i < N; i++)
    {
        int j0 = i - (ia[i + 1] - ia[i]);
        for (int ii = ia[i]; ii < ia[i + 1]; ii++)
        {
            int j = ii - ia[i] + j0;
            double sum_l = 0, sum_u = 0;
            if (ia[j] < ia[j + 1])
            {
                int j0j = j - (ia[j + 1] - ia[j]);
                int jjbeg = j0 < j0j ? j0j : j0; // max (j0, j0j)
                int jjend = j < i - 1 ? j : i - 1; // min (j, i - 1)
                for (int k = 0; k < jjend - jjbeg; k++)
                {
                    int ind_prev = ia[j] + jjbeg - j0j + k;
                    int ind_now = ia[i] + jjbeg - j0 + k;
                    sum_l += au[ind_prev] * al[ind_now];
                    sum_u += au[ind_now] * al[ind_prev];
                }
            }
            al[ii] -= sum_l;
            au[ii] -= sum_u;
            if (abs(di[j]) < eps) // matrix hasn't LU
            {
                cout << "di[" << j << "] = " << di[j] << endl;
                return false;
            }
            au[ii] /= di[j];
            di[i] -= al[ii] * au[ii];
        }
    }
    return true;
}

int LU_solver::Solve_task()
{
    cout << "start LU decomposition" << endl;

    if (!LU())
    {
        cout << "The matrix has no LU decomposition" << endl;
       return 1;
    }
    cout << "start  forward_stroke" << endl;
    forward_stroke();
    cout << "start  revers" << endl;
    revers();
    cout << "end LU" << endl;

    return 0;
}
#include "Generate.h"

void Create_material_file(std::string path, double lambda, double sigma, double hi)
{
    std::ofstream out("material.txt");
    out << lambda << " " << sigma << " " << hi;
    out.close();
}

void Create_grid_from_points(std::string path, double lambda, double sigma, double hi)
{
    std::ofstream out;
    out.precision(15);  
    std::vector<double> all_X, all_Y, all_Z;
    std::ifstream in(path + "grid.txt");
    double X, Y, Z, kx, ky, kz;
    int Nx, Ny, Nz;
    int count_x, count_y, count_z;
    in >> count_x >> count_y >> count_z;
    all_X.resize(count_x);
    all_Y.resize(count_y);
    all_Z.resize(count_z);

    in >> all_X[0] >> all_Y[0] >> all_Z[0];
    for (int curr_count_x = 0; curr_count_x < count_x - 1; )
    {
        in >> X >> Nx >> kx;
        double hx;
        if (kx == 1)
        {
            hx = (X - all_X[curr_count_x]) / Nx;
            for (int p = 1; p < Nx; p++)
            {
                all_X[curr_count_x + p] = all_X[curr_count_x] + hx * p;
            }
            curr_count_x += Nx;
        }
        else
        {
            hx = (X - all_X[curr_count_x]) * (kx - 1) / (pow(kx, Nx) - 1);
            for (int p = 0; p < Nx - 1; curr_count_x++, p++)
            {
                all_X[curr_count_x + 1] = all_X[curr_count_x] + hx * pow(kx, p);
            }
            curr_count_x++;
        }
        all_X[curr_count_x] = X;
    }
    for (int curr_count_y = 0; curr_count_y < count_y - 1; )
    {
        in >> Y >> Ny >> ky;
        double hy;
        if (ky == 1)
        {
            hy = (Y - all_Y[curr_count_y]) / Ny;
            for (int p = 1; p < Ny; p++)
            {
                all_Y[curr_count_y + p] = all_Y[curr_count_y] + hy * p;
            }
            curr_count_y += Ny;
        }
        else
        {
            hy = (Y - all_Y[curr_count_y]) * (ky - 1) / (pow(ky, Ny) - 1);
            for (int p = 0; p < Ny - 1; curr_count_y++, p++)
            {
                all_Y[curr_count_y + 1] = all_Y[curr_count_y] + hy * pow(ky, p);
            }
            curr_count_y++;
        }
        all_Y[curr_count_y] = Y;
    }
    for (int curr_count_z = 0; curr_count_z < count_z - 1; )
    {
        in >> Z >> Nz >> kz;
        double hy;
        if (kz == 1)
        {
            hy = (Z - all_Z[curr_count_z]) / Nz;
            for (int p = 1; p < Nz; p++)
            {
                all_Z[curr_count_z + p] = all_Z[curr_count_z] + hy * p;
            }
            curr_count_z += Nz;
        }
        else
        {
            hy = (Z - all_Z[curr_count_z]) * (kz - 1) / (pow(kz, Nz) - 1);
            for (int p = 0; p < Nz - 1; curr_count_z++, p++)
            {
                all_Z[curr_count_z + 1] = all_Z[curr_count_z] + hy * pow(kz, p);
            }
            curr_count_z++;
        }
        all_Z[curr_count_z] = Z;
    }
    in.close();
    out.open("xyz.txt");
    for (int i = 0; i < count_z; i++)
    {
        for (int k = 0; k < count_y; k++)
        {
            for (int j = 0; j < count_x; j++)
            {
                out << all_X[j] << "\t" << all_Y[k] << "\t" << all_Z[i] << "\n\n";
            }
        }
    }
    out.close();

    // input area
    out.open("elem.txt");
    for (int i = 0; i < count_z - 1; i++)
    {
        for (int k = 0; k < count_y - 1; k++)
        {
            for (int j = 0; j < count_x - 1; j++)
            {
                out << i*count_y* count_x + k * count_x + j << " " << i * count_y * count_x + k * count_x + j + 1 << " "
                    << i * count_y * count_x + (k + 1) * count_x + j << " " << i * count_y * count_x + (k + 1) * count_x + j + 1 << " "
                    << (i+1) * count_y * count_x + k * count_x + j << " " << (i + 1) * count_y * count_x + k * count_x + j + 1 << " "
                    << (i + 1) * count_y * count_x + (k + 1) * count_x + j << " " << (i + 1) * count_y * count_x + (k + 1) * count_x + j + 1
                    << " 0\n";
            }
        }
    }
    out.close();

    out.open("info.txt");
    out << count_x * count_y * count_z << " 1 " << (count_x - 1) * (count_y - 1) * (count_z - 1) << " 1 ";
    out.close();

    Create_material_file(path, lambda, sigma, hi);

    out.open("S1.txt");
    out << 2*(count_y * count_x + (count_z - 2) * count_y + (count_z - 2) * (count_x - 2)) << "\n";
    // верх/низ
    for (int k = 0; k < count_y; k++)
    {
        for (int j = 0; j < count_x; j++)
        {
            out << k * count_x + j << " " << 
                (count_z - 1) * count_y * count_x + k * count_x + j << " ";
        }
    }
    // зад/перед
    for (int i = 1; i < count_z - 1; i++)
    {
        for (int k = 0; k < count_y; k++)
        {
            out << i * count_y * count_x + k * count_x << " "
                << i * count_y * count_x + k * count_x + (count_x - 1)
                << " ";
        }
    }
    // право/лево
    for (int i = 1; i < count_z - 1; i++)
    {
            for (int j = 1; j < count_x - 1; j++)
            {
                out << i * count_y * count_x + j << " " 
                    << i * count_y * count_x + (count_y - 1) * count_x + j
                    << " ";
            }
    }
}


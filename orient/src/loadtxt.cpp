#include <algorithm>
#include <cstdlib>
#include <fstream>
#include "loadtxt.h"

//TODO if cols is an empty vector, get all columns from the file
//TODO error checking

Loadtxt::Loadtxt(std::string file_name, std::vector<int> cols, int skiprows)
{
    std::sort(cols.begin(), cols.end());
    n_cols = cols.size();
    const int tmp_number_of_rows = 16384;
    int n_rows_alloc = tmp_number_of_rows;
    buffer = (double*)malloc(n_cols * sizeof(double) * n_rows_alloc);
    std::ifstream file(file_name);
    if (!file.good()) throw std::runtime_error("Input file not found (" + file_name + ")");
    std::string line;
    for (int row = 0; row < skiprows; row++) getline(file, line);
    int row = -1;
    while (getline(file, line)) {
        if (line[line.find_first_not_of(whitespaces)]=='#') continue;
        if (++row >= n_rows_alloc) {
            n_rows_alloc += tmp_number_of_rows;
            buffer = (double*)realloc(buffer, n_cols * sizeof(double) * n_rows_alloc);
        }
        line_to_buf(cols, line, buffer + row*n_cols);
    }
    file.close();
    buffer = (double*)realloc(buffer, n_cols * sizeof(double) * (++row));
    n_rows = row;
}

Loadtxt::~Loadtxt()
{
    free(buffer);
}

std::vector<std::vector<double>> Loadtxt::get_cols()
{
    std::vector<std::vector<double>> data(n_cols);
    for (int col=0; col<n_cols; col++) data[col] = std::vector<double>(n_rows);
    for (int row=0; row<n_rows; row++) {
        for (int col=0; col<n_cols; col++) {
            data[col][row] = buffer[row*n_cols + col];
        }
    }
    return data;
}

void Loadtxt::line_to_buf(std::vector<int> cols, std::string line, double *buffer)
{
    int n_cols = cols.size();
    line = line.substr(line.find_first_not_of(whitespaces));
    auto pos = line.find_first_of(whitespaces, 1);
    int col=0, idx=0;
    std::vector<double> data;
    while (pos != std::string::npos) {
        std::string num_as_string;
        num_as_string = line.substr(0, pos);
        pos = line.find_first_not_of(whitespaces, pos);
        line = line.substr(pos, std::string::npos);
        pos = line.find_first_of(whitespaces, 1);
        if (col++ == cols[idx]) buffer[idx++] = std::stod(num_as_string);
    }
    if (col++ == cols[idx]) buffer[idx++] = std::stod(line);
    if (idx < n_cols) throw std::runtime_error("Read fewer columns than expected");
}

// Below is a deomonstration. The file has multiple columns but we are only
// interested in the second and fourth. We pass the file name and the column
// vector {1, 3} (since numbering starts at zero) and immediately call the
// get_cols() member, which returns a vector of vectors. We then manually assign
// the data members into named vectors.
//
// #include <iostream>
// int main()
// {

//     auto data = Loadtxt("file.dat", {1, 3}).get_cols();
//     auto& time  = data[0];
//     auto& value = data[1];
//     for (size_t i=0; i<time.size(); i++) {
//         std::cout << "time: " << time[i] << "   value: " << value[i] << std::endl;
//     }
// }
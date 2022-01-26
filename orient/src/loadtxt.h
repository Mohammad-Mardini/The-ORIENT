#pragma once
#include <string>
#include <vector>
class Loadtxt {
public:
    Loadtxt(std::string file_name, std::vector<int> cols, int skiprows=0);
    ~Loadtxt();
    std::vector<std::vector<double>> get_cols();
private:
    const char *whitespaces = " \t";
    void line_to_buf(std::vector<int> cols, std::string line, double *buffer);
    double *buffer;
    int n_rows, n_cols;
};
#pragma once

#include <array>
#include <gsl/gsl_spline.h>
#include <vector>

void set_data_dir(const std::string &path);

class Interp {
    // This is a wrapper around GSL spline interpolation. I tried to use
    // boost::math::interpolators but as of version 1.72 none were suitable:
    // barycentric_rational is the one suitalbe for non-uniform sampling but it
    // is very slow. I also tried to resample the data uniformly using
    // barycentric rational interpolation and then using cardinal_cubic_b_spline
    // on the result, but was still slower than GSL.
public:
    Interp(std::vector<double>& x, std::vector<double>& y);
    Interp() {}
    inline double operator()(const double x) const;
    void freeze(const double x);
private:
    gsl_interp_accel *acc;
    gsl_spline *spline;
    double x_min, x_max;
    bool frozen;
    double frozen_value;
};

class Galaxy {
public:
    constexpr static unsigned num_params = 7;
    Galaxy(std::string file_name);
    void func(const std::array<double, 6> &y, std::array<double, 6> &f, const double t) const;
    std::tuple<double,double> get_time_limits() {
        return { t_min, t_max };
    }
    std::array<double, num_params> get_fit_params(const double t);
    std::vector<std::array<double, num_params>> get_fit_params(const std::vector<double> t);
    void freeze_parameter(const int i, const double t);
    void freeze_parameter(const std::string& name, const double t);
    void freeze_all_parameters(const double t);
private:
    double t_min, t_max;
    std::array<Interp, num_params> interp;
    const std::unordered_map<std::string, int> parameter_idx_by_name;
};

std::vector<std::array<double,6>> integrate(const Galaxy &galaxy, const std::array<double,6> y0, const double t_min, const double t_max, const double stride_size, const size_t max_size);

#include <algorithm>
#include <boost/numeric/odeint.hpp>
#include <Eigen/Geometry> 
#include <fstream>

#ifdef __APPLE__
#warning I don't know anything about MacOS but on Mojave this doesn't seem to be working
#else
#include <filesystem>
#endif

#include <iostream>
#include <numeric>
#include <string>
#include <stdexcept>

#include "loadtxt.h"
#include "replay.h"

using double3 = Eigen::Vector3d;

// Our units are {kiloparsec, solar mass, gigayear}
constexpr double G = 4.498317481097514e-06;

std::string data_dir; //TODO explain

void set_data_dir(const std::string &path)
{
    data_dir = path;
}

Interp::Interp(std::vector<double>& x, std::vector<double>& y)
{
    acc    = gsl_interp_accel_alloc();
    spline = gsl_spline_alloc(gsl_interp_cspline, x.size());
    gsl_spline_init(spline, x.data(), y.data(), x.size());
    x_min = x.front();
    x_max = x.back();
}

inline double Interp::operator()(double x) const
{
    return gsl_spline_eval(spline, x, acc);
}

double3 plummer(const double M, const double b, const double3& pos)
{
    double r2 = pos.squaredNorm() + b*b;
    double r  = sqrt(r2);
    double r3_inv = 1/(r*r2);
    return -G*M*pos*r3_inv;
}

double3 nfw(const double rho_0, const double b, const double3& pos)
{
    double r2 = pos.squaredNorm();
    double r  = sqrt(r2);
    double tmp = -4*M_PI*G*rho_0*b*b*b*(log1p(r/b) - r/(b+r))/(r2*r);
    return pos*tmp;
}

double3 miyamoto_nagai(const double M, const double a, const double b, const double phi, const double theta, const double3& pos)
{
    // Construct new z-axis from the angles
    double cos_theta = cos(theta);
    double sin_theta = sqrt(1-cos_theta*cos_theta);
    double cos_phi   = cos(phi);
    double sin_phi   = sqrt(1-cos_phi*cos_phi);
    Eigen::Vector3d old_z_axis = {cos_phi*sin_theta, sin_phi*sin_theta, cos_theta};

    // Construct rotation
    auto rot = Eigen::Quaternion<double>::FromTwoVectors(old_z_axis, Eigen::Vector3d::UnitZ());

    // Rotate the position vector
    auto new_pos = rot * pos;

    // Calculate acceleration in new frame
    Eigen::Vector3d acc_in_rotated_frame;
    auto z_tmp  = sqrt(new_pos[2]*new_pos[2] + b*b);
    auto r2_tmp = new_pos[0]*new_pos[0] + new_pos[1]*new_pos[1] + (z_tmp + a)*(z_tmp + a);
    auto r_tmp  = sqrt(r2_tmp);
    auto tmp = G*M / (r2_tmp*r_tmp);
    acc_in_rotated_frame[0] = - tmp * new_pos[0];
    acc_in_rotated_frame[1] = - tmp * new_pos[1];
    acc_in_rotated_frame[2] = - tmp * new_pos[2] * (z_tmp + a)/z_tmp;

    // Return to original frame
    return rot.inverse() * acc_in_rotated_frame;
}

Galaxy::Galaxy(std::string file_name)
{
    #ifndef __APPLE__
    namespace fs = std::filesystem;
    if (!fs::is_regular_file(file_name)) {
        file_name = fs::path(data_dir) / fs::path(file_name);
        if (!fs::is_regular_file(file_name))
            throw std::runtime_error("File not found");
    }
    #endif

    auto data = Loadtxt(file_name, {1, 2, 3, 4, 5, 6, 7, 8}).get_cols();
    auto& t_data = data[0];
    t_min = t_data.front();
    t_max = t_data.back();
    for (int i=0; i<7; i++) interp[i] = Interp(t_data, data[i+1]);
}

void Galaxy::func(const std::array<double, 6> &y, std::array<double, 6> &f, const double t) const
{
    f[0] = y[3]; // vx -> x'
    f[1] = y[4]; // vy -> y'
    f[2] = y[5]; // vz -> z'
    double phi         = interp[0](t);
    double theta       = interp[1](t);
    double m_disk      = interp[2](t);
    double a_disk      = interp[3](t);
    double b_disk      = interp[4](t);
    double rho_0_nfw   = interp[5](t);
    double b_nfw       = interp[6](t);

    double3 pos(y.data());
    double3 acc_disk = miyamoto_nagai(m_disk, a_disk, b_disk, phi, theta, pos);
    double3 acc_halo = nfw(rho_0_nfw, b_nfw, pos);
    for (int i=0; i<3; i++) f[3+i] = acc_disk[i] + acc_halo[i];
};


std::vector<std::array<double,6>> integrate(const Galaxy &galaxy, const std::array<double,6> y0, const double t_min, const double t_max, const double stride_size, const size_t max_size = 0)
{
    using namespace boost::numeric::odeint;
    using Coordinates = std::array<double, 6>;
    auto stepper = bulirsch_stoer<Coordinates>(1E-7, 0);
    auto function_wrapper = [&galaxy](const Coordinates &x, Coordinates &dxdt, const double t) { return galaxy.func(x, dxdt, t); };
    size_t stride_count = (t_max-t_min) / stride_size;
    if (stride_count*stride_size > t_max-t_min) stride_count++;
    if (stride_count > max_size) stride_count = max_size;
    std::vector<Coordinates> y(stride_count);
    Coordinates y_current;
    std::copy(begin(y0), end(y0), begin(y_current));
    std::copy(begin(y0), end(y0), begin(y[0]));
    double t = t_min;
    const double h = std::min(1./4096., stride_size);
    for (size_t i=1; i<stride_count; i++) {
        integrate_adaptive(stepper, function_wrapper, y_current, t, t+stride_size, h);
        // NOTE h here is just a recommended initial step size for the stepper,
        // the actual step is adapted as needed. Since the integration is
        // interrupted here in order for the data to be saved, the result
        // somewhat depends on stride_size.
        std::copy(begin(y_current), end(y_current), begin(y[i]));
        t += stride_size;
    }
    return y;
}

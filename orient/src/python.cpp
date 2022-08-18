#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "replay.h"

extern std::string data_dir;

namespace py = pybind11;

PYBIND11_MODULE(orient, m) {
    m.doc() = "cOsmologically deRIved timE-varyiNg Galactic poTentials";

    py::class_<Galaxy>(m, "Galaxy")
        .def(py::init<std::string>())
        .def("get_time_limits", &Galaxy::get_time_limits)
        .def("get_fit_params", py::overload_cast<const double>(&Galaxy::get_fit_params))
        .def("get_fit_params", py::overload_cast<const std::vector<double>>(&Galaxy::get_fit_params))
        .def("freeze_parameter", py::overload_cast<const int, const double>(&Galaxy::freeze_parameter))
        .def("freeze_parameter", py::overload_cast<const std::string&, const double>(&Galaxy::freeze_parameter))
        .def("freeze_all_parameters", &Galaxy::freeze_all_parameters);

    m.def("integrate", [](const Galaxy &galaxy, const std::array<double,6> y0, const double t_min, const double t_max, const double stride_size, const size_t max_size)
        {
            return (py::array)py::cast(integrate(galaxy, y0, t_min, t_max, stride_size, max_size));
        }, "just a test");
    
    m.def("set_data_dir", &set_data_dir);
}
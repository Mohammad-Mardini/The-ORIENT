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
        .def("get_fit_params", &Galaxy::get_fit_params);

    m.def("integrate", [](const Galaxy &galaxy, const std::array<double,6> y0, const double t_min, const double t_max, const double stride_size, const size_t max_size)
        {
            return (py::array)py::cast(integrate(galaxy, y0, t_min, t_max, stride_size, max_size));
        }, "just a test");
    
    m.def("set_data_dir", &set_data_dir);
}
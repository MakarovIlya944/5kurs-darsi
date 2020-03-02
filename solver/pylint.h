#pragma once
#include "net.hpp"
#include "matrix.hpp"


#ifdef _DEBUG

#else

namespace py = pybind11;

PYBIND11_MODULE(cpp_solver, m) {
	py::class_<NetGenerator>(m, "NetGenerator")
		.def(py::init<>())
		.def("GenerateFromFiles", &NetGenerator::GenerateFromFiles, R"pbdoc(FNameN - net.txt, FNameB - border.txt, FNameT - time.txt)pbdoc")
		.def("test", &NetGenerator::test, R"pbdoc(FNameN - net.txt, FNameB - border.txt, FNameT - time.txt)pbdoc")
		.def("Generate", &NetGenerator::Generate, R"pbdoc(
Points at r
Points at z
Default step
Step and koefficient per elemens at r
Step and koefficient per elemens at z
Materials border
)pbdoc")
//.def_readwrite("isLoging", &NetGenerator::isLoging)
;

	py::class_<Net>(m, "Net")
		.def(py::init<>())
		.def_readwrite("rNum", &Net::rNum, R"pbdoc(rNum - num of points at r)pbdoc")
		.def_readwrite("zNum", &Net::zNum, R"pbdoc(zNum - num of points at z)pbdoc")
		.def_readwrite("Num", &Net::Num, R"pbdoc(Num - num of points)pbdoc")
		.def_readwrite("NumT", &Net::NumT, R"pbdoc(NumT - num of time layers)pbdoc")
		.def_readwrite("NumEl", &Net::NumEl, R"pbdoc(NumEl - num of elements)pbdoc")
		.def_readwrite("GlobalNet", &Net::GlobalNet, R"pbdoc(GlobalNet - koords of all points (r,z))pbdoc")
		.def_readwrite("GlobalTime", &Net::GlobalTime, R"pbdoc(GlobalTime - all time layers)pbdoc")
		.def_readwrite("nvtr", &Net::nvtr, R"pbdoc(nvtr - all elements)pbdoc")
		.def_readwrite("nvcat", &Net::nvcat, R"pbdoc(nvcat - all materials)pbdoc")
		.def_readwrite("nvr1", &Net::nvr1, R"pbdoc(first border condition)pbdoc")
		.def_readwrite("nvr2", &Net::nvr2, R"pbdoc(second border condition)pbdoc")
		.def_readwrite("nvr3", &Net::nvr3, R"pbdoc(third border condition)pbdoc");

	py::class_<Pointd>(m, "Pointd")
		.def(py::init<>())
		.def_readwrite("r", &Pointd::r)
		.def_readwrite("z", &Pointd::z)
		.def("+", &Pointd::operator+)
		.def("-", &Pointd::operator-)
		;

	py::class_<Pointi>(m, "Pointi")
		.def(py::init<>())
		.def_readwrite("r", &Pointi::r)
		.def_readwrite("z", &Pointi::z)
		;

#ifdef VERSION_INFO
	m.attr("__version__") = VERSION_INFO;
#else
	m.attr("__version__") = "dev";
#endif
}


#endif 
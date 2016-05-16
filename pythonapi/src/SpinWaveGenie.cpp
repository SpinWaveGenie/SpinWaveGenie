#include "SpinWaveGenie/Containers/Cell.h"
#include "SpinWaveGenie/Containers/Sublattice.h"
#include "SpinWaveGenie/Genie/Neighbors.h"
#include "pybind11/pybind11.h"
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

namespace py = pybind11;
using namespace SpinWaveGenie;

PYBIND11_PLUGIN(PySpinWaveGenie)
{
  py::module m("PySpinWaveGenie", "Python Bindings for the SpinWaveGenie Library");

  py::class_<Sublattice>(m, "Sublattice")
      .def(py::init<>())
      .def("setName", &Sublattice::setName, "set name to describe sublattice")
      .def("getName", &Sublattice::getName, "returns name of a given sublattice")
      .def("setType", &Sublattice::setType,
           "set type to describe magnetic form factor used in the calculation of intensities")
      .def("getType", &Sublattice::getType, "returns name of a given sublattice")
      .def("setMoment", &Sublattice::setMoment, "set moment in spherical coordinates r,theta,phi")
      .def("getMoment", &Sublattice::getMoment, "coordinate r")
      .def("getTheta", &Sublattice::getTheta, "coordinate theta in radians")
      .def("getPhi", &Sublattice::getPhi, "coordinate phi in radians")
      .def("getRotationMatrix", &Sublattice::getRotationMatrix, "returns rotation matrix as an Eigen::Matrix3d object")
      .def("getInverseMatrix", &Sublattice::getInverseMatrix,
           "returns inverse rotation matrix as an Eigen::Matrix3d object")
      .def("addAtom", &Sublattice::addAtom, "add atom to the sublattice")
      .def("__iter__", [](const Sublattice &s) { return py::make_iterator(s.begin(), s.end()); },
           py::keep_alive<0, 1>())
      .def("__len__", &Sublattice::size, "returns the number of atoms in this sublattice");

  py::class_<Cell>(m, "Cell")
      .def(py::init<Cell>())
      .def("setBasisVectors",
           static_cast<void (Cell::*)(double, double, double, double, double, double)>(&Cell::setBasisVectors),
           "Set basis vectors from parameters")
      .def("setBasisVectors", static_cast<void (Cell::*)(double, const Matrix3 &)>(&Cell::setBasisVectors),
           "Set basis vectors from parameters")
      .def("getBasisVectors", &Cell::getBasisVectors, "get basis vectors")
      .def("getReciprocalVectors", &Cell::getReciprocalVectors, "get basis vectors")
      .def("addSublattice", &Cell::addSublattice, "Add sublattice to cell")
      .def("getSublattice", static_cast<Sublattice &(Cell::*)(const std::string &name)>(&Cell::getSublattice),
           "Returns sublattice object")
      .def("__getitem__",
           [](const Cell &s, size_t i) {
             if (i >= s.size())
               throw py::index_error();
             return s[i];
           })
      .def("addAtom", &Cell::addAtom, "Add atom to sublattice name at a given position")
      .def("getPosition", &Cell::getPosition, "Returns the position where sublattice name is stored.")
      .def("__len__", &Cell::size, "Returns the number of sublattices in the cell")
      .def("__iter__", [](const Cell &s) { return py::make_iterator(s.begin(), s.end()); }, py::keep_alive<0, 1>());

  py::class_<Neighbors>(m, "Neighbors")
      .def(py::init<Neighbors>())
      .def("empty", &Neighbors::empty, "Returns whether of not neighbors have been calculated previously")
      .def("findNeighbors", &Neighbors::findNeighbors,
           "Finds neighbors of two sublattices between distances min and max.")
      .def("__len__", &Neighbors::size, " Get the number of neighbors")
      .def("getGamma", &Neighbors::getGamma, "Get variable")
      .def("__iter__", [](const Neighbors &s) { return py::make_iterator(s.begin(), s.end()); },
           py::keep_alive<0, 1>());

  return m.ptr();
}

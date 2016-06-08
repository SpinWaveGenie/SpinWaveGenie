#include "SpinWaveGenie/Containers/Cell.h"
#include "SpinWaveGenie/Containers/Results.h"
#include "SpinWaveGenie/Containers/Sublattice.h"
#include "SpinWaveGenie/Genie/Neighbors.h"
#include "SpinWaveGenie/Genie/SpinWave.h"
#include "SpinWaveGenie/Genie/SpinWaveBuilder.h"
#include "SpinWaveGenie/Interactions/Interaction.h"
#include "SpinWaveGenie/Interactions/InteractionFactory.h"
#include "pybind11/eigen.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

namespace py = pybind11;
using namespace SpinWaveGenie;

PYBIND11_PLUGIN(python_SpinWaveGenie)
{
  py::module m("python_SpinWaveGenie", "Python Bindings for the SpinWaveGenie Library");

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
      .def(py::init<>())
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
      .def(py::init<>())
      .def("empty", &Neighbors::empty, "Returns whether of not neighbors have been calculated previously")
      .def("findNeighbors", &Neighbors::findNeighbors,
           "Finds neighbors of two sublattices between distances min and max.")
      .def("__len__", &Neighbors::size, " Get the number of neighbors")
      .def("getGamma", &Neighbors::getGamma, "Get variable")
      .def("__iter__", [](const Neighbors &s) { return py::make_iterator(s.begin(), s.end()); },
           py::keep_alive<0, 1>());

  py::class_<Interaction>(m, "Interaction")
      .def("sublattices",&Interaction::sublattices,"get Sublattices associated with this Interaction.")
      .def("getName",&Interaction::getName,"get Name associated with this Interaction.")
      .def("updateValue",&Interaction::updateValue,"update Values associated with this Interaction.");

  py::class_<InteractionFactory>(m, "InteractionFactory")
      .def(py::init<>())
      .def("getExchange",&InteractionFactory::getExchange,"Construct exchange interaction term.")
      .def("getDzyaloshinskiiMoriya",&InteractionFactory::getDzyaloshinskiiMoriya,"Construct Dzyaloshinskii-Moriya interaction term")
      .def("getAnisotropy",&InteractionFactory::getAnisotropy,"Construct single-ion anisotropy term.")
      .def("getMagneticField",&InteractionFactory::getMagneticField,"Construct magnetic field term.");

  py::class_<Point>(m, "Point")
      .def_readonly("frequency", &Point::frequency, "Frequency associated with a given excitation (in meV).")
      .def_readonly("intensity", &Point::intensity, " Measurable intensity of a given excitation (arb. units).");

  py::class_<Results>(m, "Results")
      .def("insert", &Results::insert, "Insert Point struct into container.")
      .def("__len__", &Results::size, " size of Results container.")
      .def("__iter__", [](const Results &s) { return py::make_iterator(s.begin(), s.end()); }, py::keep_alive<0, 1>())
      .def("sort", &Results::sort, "Sort Results by frequency.")
      .def("clear", &Results::clear, "Clear results container so that the size is zero.")
      .def("uniqueSolutions", &Results::uniqueSolutions, "Filter Points by combining those with identical frequencies.")
      .def("significantSolutions", &Results::significantSolutions,
           "Filter Points by removing those without significant intensity.");

  py::class_<SpinWaveBuilder>(m, "SpinWaveBuilder")
      .def(py::init<>())
      .def(py::init<const Cell &>())
      .def("updateCell", &SpinWaveBuilder::updateCell, "")
      .def("addInteraction",
           static_cast<void (SpinWaveBuilder::*)(const Interaction &)>(&SpinWaveBuilder::addInteraction),
           "Add a magnetic interaction to the model.")
      .def("updateInteraction", &SpinWaveBuilder::updateInteraction, "Update interaction parameter.")
      .def("getEnergy", &SpinWaveBuilder::getEnergy, "Get the classical energy associated with the model.")
      .def("getFirstOrderTerms", &SpinWaveBuilder::getFirstOrderTerms,
           "Get coefficients associated with linear terms. For a ground state, each must be zero.")
      .def("createElement", &SpinWaveBuilder::createElement,
           "Constructs a SpinWave object to calculate the spin-wave dispersion.");

  py::class_<SpinWave>(m, "SpinWave")
      .def("clearMatrix", &SpinWave::clearMatrix,
           "Clear dynamical matrix and any previously calculated frequencies and intensities.")
      .def("getCell", &SpinWave::getCell, "Access the cell used in this calculation.")
      .def("createMatrix", &SpinWave::createMatrix, "Set the location to calculate in reciprocal space.")
      .def("calculate", &SpinWave::calculate, "Calculate spin-wave frequencies and intensities.")
      .def("getPoints", &SpinWave::getPoints, " Get calculated frequencies and intensities.");

  return m.ptr();
}

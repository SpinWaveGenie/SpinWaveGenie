#include "SpinWaveGenie/Containers/Cell.h"
#include "SpinWaveGenie/Containers/Energies.h"
#include "SpinWaveGenie/Containers/PointsAlongLine.h"
#include "SpinWaveGenie/Containers/Results.h"
#include "SpinWaveGenie/Containers/Sublattice.h"
#include "SpinWaveGenie/Containers/ThreeVectors.h"
#include "SpinWaveGenie/Containers/UniqueThreeVectors.h"
#include "SpinWaveGenie/Genie/Neighbors.h"
#include "SpinWaveGenie/Genie/SpinWave.h"
#include "SpinWaveGenie/Genie/SpinWaveBuilder.h"
#include "SpinWaveGenie/Interactions/Interaction.h"
#include "SpinWaveGenie/Interactions/InteractionFactory.h"
#include "SpinWaveGenie/Plot/Plot.h"
#include "pybind11/eigen.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
// #include "SpinWaveGenie/Plot/IntegrateEnergy.h "
// #include "SpinWaveGenie/Plot/IntegrateThetaPhi.h"
// #include "SpinWaveGenie/Plot/SpinWaveDispersion.h"
// #include "SpinWaveGenie/Plot/SpinWavePlot.h"
// #include "SpinWaveGenie/Plot/TwoDGaussian.h"
// #include "SpinWaveGenie/Plot/TwoDimensionalGaussian.h"

namespace py = pybind11;
using namespace SpinWaveGenie;

PYBIND11_MODULE(python_SpinWaveGenie, m)
{
  m.doc() = "Python Bindings for the SpinWaveGenie Library";

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
      .def("setBasisVectors", py::overload_cast<double, double, double, double, double, double>(&Cell::setBasisVectors),
           "Set basis vectors from parameters")
      .def("setBasisVectors", py::overload_cast<double, const Eigen::Matrix3d &>(&Cell::setBasisVectors),
           "Set basis vectors from parameters")
      .def("getBasisVectors", &Cell::getBasisVectors, "get basis vectors")
      .def("getReciprocalVectors", &Cell::getReciprocalVectors, "get basis vectors")
      .def("addSublattice", &Cell::addSublattice, "Add sublattice to cell")
      .def("getSublattice", py::overload_cast<const std::string &>(&Cell::getSublattice), "Returns sublattice object")
      .def("__getitem__",
           [](const Cell &s, size_t i) {
             if (i >= s.size())
             {
               throw py::index_error();
             }
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
      .def("sublattices", &Interaction::sublattices, "get Sublattices associated with this Interaction.")
      .def("getName", &Interaction::getName, "get Name associated with this Interaction.")
      .def("updateValue", &Interaction::updateValue, "update Values associated with this Interaction.");

  py::class_<InteractionFactory>(m, "InteractionFactory")
      .def(py::init<>())
      .def("getExchange", &InteractionFactory::getExchange, "Construct exchange interaction term.")
      .def("getDzyaloshinskiiMoriya", &InteractionFactory::getDzyaloshinskiiMoriya,
           "Construct Dzyaloshinskii-Moriya interaction term")
      .def("getAnisotropy",
           py::overload_cast<const std::string &, double, const Eigen::Vector3d &, const std::string &>(
               &InteractionFactory::getAnisotropy),
           "Construct single-ion anisotropy term.")
      .def("getAnisotropy",
           py::overload_cast<const std::string &, const Eigen::Matrix3d &, const std::string &>(
               &InteractionFactory::getAnisotropy),
           "Construct single-ion anisotropy term.")
      .def("getMagneticField", &InteractionFactory::getMagneticField, "Construct magnetic field term.");

  py::class_<Point>(m, "Point")
      .def_readonly("frequency", &Point::frequency, "Frequency associated with a given excitation (in meV).")
      .def_readonly("intensity", &Point::intensity, " Measurable intensity of a given excitation (arb. units).");

  py::class_<PointsAlongLine>(m, "PointsAlongLine")
      .def(py::init<>())
      .def("setFirstPoint", &PointsAlongLine::setFirstPoint, "Set the first point in a line of points")
      .def("setFinalPoint", &PointsAlongLine::setFinalPoint, "Set the last point in a line of points")
      .def("setNumberPoints", &PointsAlongLine::setNumberPoints, "Set the Number of points")
      .def("getPoints", &PointsAlongLine::getPoints, "Get the points as described above");

  py::class_<Energies>(m, "Energies").def(py::init<>()).def(py::init<double, double, std::size_t>());

  py::class_<Results>(m, "Results")
      .def("insert", &Results::insert, "Insert Point struct into container.")
      .def("__len__", &Results::size, " size of Results container.")
      .def("__iter__", [](const Results &s) { return py::make_iterator(s.begin(), s.end()); }, py::keep_alive<0, 1>())
      .def("__getitem__",
           [](const Results &s, size_t i) {
             if (i >= s.size())
             {
               throw py::index_error();
             }
             return s[i];
           })
      .def("sort", &Results::sort, "Sort Results by frequency.")
      .def("clear", &Results::clear, "Clear results container so that the size is zero.")
      .def("uniqueSolutions", &Results::uniqueSolutions, "Filter Points by combining those with identical frequencies.")
      .def("significantSolutions", &Results::significantSolutions,
           "Filter Points by removing those without significant intensity.");

  py::class_<SpinWaveBuilder>(m, "SpinWaveBuilder")
      .def(py::init<>())
      .def(py::init<const Cell &>())
      .def("updateCell", &SpinWaveBuilder::updateCell, "")
      .def("addInteraction", py::overload_cast<const Interaction &>(&SpinWaveBuilder::addInteraction),
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
      .def("getPoints", &SpinWave::getPoints, "Get calculated frequencies and intensities.");

  py::class_<ThreeVectors<double>>(m, "ThreeVectors")
      .def(py::init<>())
      .def("clear", &ThreeVectors<double>::clear, "Clear container so that the size is zero.")
      .def("empty", &ThreeVectors<double>::empty, "check if there are no elements in this container.")
      .def("insert", &ThreeVectors<double>::insert, "insert an array of 3 elements into container.")
      .def("__len__", &ThreeVectors<double>::size, "returns the size of the container.")
      .def("__iter__", [](const ThreeVectors<double> &s) { return py::make_iterator(s.begin(), s.end()); },
           py::keep_alive<0, 1>());

  py::class_<ThreeVectors<std::complex<double>>>(m, "ComplexThreeVectors")
      .def(py::init<>())
      .def("clear", &ThreeVectors<std::complex<double>>::clear, "Clear container so that the size is zero.")
      .def("empty", &ThreeVectors<std::complex<double>>::empty, "check if there are no elements in this container.")
      .def("insert", &ThreeVectors<std::complex<double>>::insert, "insert an array of 3 elements into container.")
      .def("__len__", &ThreeVectors<std::complex<double>>::size, "returns the size of the container.")
      .def("__iter__",
           [](const ThreeVectors<std::complex<double>> &s) { return py::make_iterator(s.begin(), s.end()); },
           py::keep_alive<0, 1>());

  py::class_<UniqueThreeVectors<double>, ThreeVectors<double>>(m, "UniqueThreeVectors")
      .def(py::init<>())
      .def("insert", &UniqueThreeVectors<double>::insert, "insert an array of 3 elements into container.")
      .def("__eq__", &UniqueThreeVectors<double>::operator==, "check if two arrays contain the same elements");

  py::class_<OneDimensionalShapes>(m, "OneDimensionalShapes")
      .def("getFunction", &OneDimensionalShapes::getFunction, "Get the function used");

  py::class_<OneDimensionalFactory>(m, "OneDimensionalFactory")
      .def(py::init<>())
      .def("getGaussian", py::overload_cast<double, double>(&OneDimensionalFactory::getGaussian),
           "Given a FWHM and a tolerance provide a 1D Gaussian object")
      .def("getGaussian", py::overload_cast<const std::array<double, 4> &, double>(&OneDimensionalFactory::getGaussian),
           "Given a FWHM and a tolerance provide a 1D Gaussian object")
      .def("getLorentzian", &OneDimensionalFactory::getLorentzian,
           "Given a FWHM and a tolerance provide a 1D Lorentzian object")
      .def("getPseudoVoigt", &OneDimensionalFactory::getPseudoVoigt,
           "Given an eta, a FWHM ,and a tolerance provide a 1D PseudoVoigh object");

  py::class_<SpinWavePlot>(m, "SpinWavePlot")
      .def("getCell", &SpinWavePlot::getCell, py::return_value_policy::reference_internal, "retrieve Cell")
      .def("getCut", &SpinWavePlot::getCut, "retrieve the Cut")
      .def("getEnergies", &SpinWavePlot::getEnergies, py::return_value_policy::reference_internal, "retrieve energies")
      .def("setEnergies", &SpinWavePlot::setEnergies, "set energies");

  py::class_<EnergyResolutionFunction, SpinWavePlot>(m, "EnergyResolutionFunction")
      .def(py::init<>())
      .def(py::init<const OneDimensionalShapes &, const SpinWave &, const Energies &>())
      .def("setResolutionFunction",
           py::overload_cast<const OneDimensionalShapes &>(&EnergyResolutionFunction::setResolutionFunction),
           "Set the resolution function")
      .def("setSpinWave", &EnergyResolutionFunction::setSpinWave, "Set the spin wave model");

  py::class_<IntegrateEnergy, SpinWavePlot>(m, "IntegrateEnegy")
      .def(py::init<const SpinWavePlot &, const Energies &, double, double, int>());

  py::class_<TwoDimensionalCut>(m, "TwoDimensionalCut")
      .def(py::init<>())
      .def("setFilename", &TwoDimensionalCut::setFilename, "Set filename to save results of cut")
      .def("setPoints", &TwoDimensionalCut::setPoints, "")
      .def("setEnergyPoints", &TwoDimensionalCut::setEnergyPoints, "")
      .def("setPlotObject", py::overload_cast<const SpinWavePlot &>(&TwoDimensionalCut::setPlotObject), "")
      .def("getMatrix", &TwoDimensionalCut::getMatrix, "Get cut as a 2D array")
      .def("save", &TwoDimensionalCut::save, "Calculate and the save the result to the specified filename");

  py::class_<IntegrateThetaPhi, SpinWavePlot>(m, "IntegrateThetaPhi")
      .def(py::init<const SpinWavePlot &, double, int>());
}

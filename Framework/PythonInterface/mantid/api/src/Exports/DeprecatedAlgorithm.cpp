#include "MantidAPI/DeprecatedAlgorithm.h"
#include <boost/python/class.hpp>

using namespace Mantid::API;
using namespace boost::python;

void export_DeprecatedAlgorithm() {

  class_<DeprecatedAlgorithm, boost::noncopyable>(
      "DeprecatedAlgorithm", "Base class of deprecated algorithms")

      .def("useAlgorithm", &DeprecatedAlgorithm::useAlgorithm,
           (arg("self"), arg("alg"), arg("version")),
           "Set the name of the algorithm to use instead")

      .def("deprecatedDate", &DeprecatedAlgorithm::deprecatedDate,
           (arg("self"), arg("date")), "Set the deprecation date");
}

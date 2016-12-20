#include "MantidAPI/DeprecatedAlgorithm.h"
#include <boost/python/class.hpp>

using namespace Mantid::API;
using namespace boost::python;

void export_DeprecatedAlgorithm() {
  class_<DeprecatedAlgorithm, bases<Algorithm>, boost::noncopyable>(
      "DeprecatedAlgorithm", "Base class of deprecated algorithms")
      .def("useAlgorithm", &DeprecatedAlgorithm::useAlgorithm,
           (arg("self"), arg("alg")),
           "Set the name of the algorithm to use instead")

      .def("deprecationMsg", &DeprecatedAlgorithm::deprecationMsg,
           (arg("self"), arg("message")), "Set the deprecation message")

      .def("deprecatedDate", &DeprecatedAlgorithm::deprecatedDate,
           (arg("self"), arg("date")), "Set the deprecation date");
}

#include "MantidAPI/DeprecatedAlgorithm.h"
#include "MantidPythonInterface/kernel/GetPointer.h"
#include <boost/python/class.hpp>

using namespace Mantid::API;
using namespace boost::python;

void export_DeprecatedAlgorithm() {

  class_<DeprecatedAlgorithm, boost::noncopyable>(
      "DeprecatedAlgorithm", "Base class of deprecated algorithms")

      .def("useAlgorithm", &DeprecatedAlgorithm::useAlgorithm,
           (arg("self"), arg("alg"), arg("version")),
           "Set the name of the algorithm to use instead")

      .def("deprecationMsg", &DeprecatedAlgorithm::deprecationMsg,
           (arg("self"), arg("ialg")), "Get the deprecation message")

      .def("deprecatedDate", &DeprecatedAlgorithm::deprecatedDate,
           (arg("self"), arg("date")), "Set the deprecation date");
}

//#include "MantidPythonInterface/api/PythonAlgorithm/DeprecatedAdapter.h"
//#include <boost/python/register_ptr_to_python.hpp>
//using Mantid::PythonInterface::DeprecatedAdapter;
//GET_POINTER_SPECIALIZATION(DeprecatedAlgorithm)
//register_ptr_to_python<boost::shared_ptr<DeprecatedAlgorithm>>();
//.add_property("deprecatedDate", &DeprecatedAlgorithm::deprecatedDate)
//.add_property("useAlgorithm", &DeprecatedAlgorithm::useAlgorithm)

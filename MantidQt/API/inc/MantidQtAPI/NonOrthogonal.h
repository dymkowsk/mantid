#ifndef MANTID_MANTIDQT_API_NON_ORTHOGONAL_H_
#define MANTID_MANTIDQT_API_NON_ORTHOGONAL_H_

#include "MantidQtAPI/DllOption.h"
#include "MantidAPI/IMDWorkspace.h"
#include "MantidKernel/Matrix.h"

namespace MantidQt {
namespace API {
void EXPORT_OPT_MANTIDQT_API
provideSkewMatrix(Mantid::Kernel::DblMatrix &skewMatrix,
                  Mantid::API::IMDWorkspace_const_sptr workspace);

bool EXPORT_OPT_MANTIDQT_API
requiresSkewMatrix(Mantid::API::IMDWorkspace_const_sptr workspace);

bool EXPORT_OPT_MANTIDQT_API
isHKLDimensions(Mantid::API::IMDWorkspace_const_sptr workspace, size_t dimX,
                size_t dimY);

size_t EXPORT_OPT_MANTIDQT_API
getMissingHKLDimensionIndex(Mantid::API::IMDWorkspace_const_sptr workspace, size_t dimX,
	size_t dimY);
}
}

#endif
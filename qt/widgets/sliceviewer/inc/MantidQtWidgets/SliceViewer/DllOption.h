#ifndef MANTIDQT_MANTIDSLICEVIEWER_DLLOPTION_H_
#define MANTIDQT_MANTIDSLICEVIEWER_DLLOPTION_H_

#include <MantidKernel/System.h>

#ifdef IN_MANTIDQT_SLICEVIEWER
#define EXPORT_OPT_MANTIDQT_SLICEVIEWER DLLExport
#else
#define EXPORT_OPT_MANTIDQT_SLICEVIEWER DLLImport
#endif

#endif // MANTIDQT_MANTIDSLICEVIEWER_DLLOPTION_H_
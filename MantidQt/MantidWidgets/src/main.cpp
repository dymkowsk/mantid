#include "MantidQtMantidWidgets/DataSelector.h"
#include <QApplication>

using namespace MantidQt::MantidWidgets;

int main(int argc, char **argv) {
 QApplication app (argc, argv);
 DataSelector d1;
 d1.show();
 DataSelector d2;
 d2.allowMultipleFiles(true);
 d2.show();
 DataSelector d3;
 d3.isForDirectory(true);
 d3.show();
 DataSelector d4;
 d4.setLabelText("Custom Label");
 d4.show();
 DataSelector d5;
 d5.setLoadBtnText("RUN");
 d5.show();
 DataSelector d6;
 d6.setAutoLoad(false);
 d6.show();
 return app.exec();
}

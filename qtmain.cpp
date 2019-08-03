#include "newickdemomenudialog.h"

#pragma GCC diagnostic push

#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#pragma GCC diagnostic ignored "-Wunused-but-set-parameter"
#include <QApplication>
#include "qtnewickdemodialog.h"
#pragma GCC diagnostic pop


int main(int argc, char *argv[])
{
  const std::vector<std::string> args { ribi::MenuDialog::ConvertArguments(argc,argv) };
  return ribi::TestNewickMenuDialog().Execute(args);
}



#include "newickdemomenudialog.h"

#pragma GCC diagnostic push



#include <QApplication>
#include "qtnewickdemodialog.h"
#pragma GCC diagnostic pop


int main(int argc, char *argv[])
{
  const std::vector<std::string> args { ribi::MenuDialog::ConvertArguments(argc,argv) };
  return ribi::TestNewickMenuDialog().Execute(args);
}



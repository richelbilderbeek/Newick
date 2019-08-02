
# C++14
CONFIG += c++14
QMAKE_CXXFLAGS += -std=c++14

QMAKE_CXXFLAGS += -Wall -Wextra -Weffc++ -Werror -std=c++14


# Debug and release mode
CONFIG += debug_and_release
CONFIG(release, debug|release) {
  message(Release mode)
  DEFINES += NDEBUG
}

# Qt
QT       += core gui
greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

include(../RibiLibraries/Boost.pri)
include(../RibiLibraries/BigInteger.pri)

include(../RibiClasses/CppAbout/CppAbout.pri)
include(../RibiClasses/CppFileIo/CppFileIo.pri)
include(../RibiClasses/CppHelp/CppHelp.pri)
include(../RibiClasses/CppMenuDialog/CppMenuDialog.pri)
include(../RibiClasses/CppQtAboutDialog/CppQtAboutDialog.pri)
include(../RibiClasses/CppQtHideAndShowDialog/CppQtHideAndShowDialog.pri)
include(../RibiClasses/CppContainer/CppContainer.pri)
include(../RibiClasses/CppFuzzy_equal_to/CppFuzzy_equal_to.pri)
include(../ManyDigitNewick/ManyDigitNewick.pri)
include(../RibiClasses/CppMultiVector/CppMultiVector.pri)
include(../BinaryNewickVector/BinaryNewickVector.pri)
include(Newick.pri)
include(../NewickVector/NewickVector.pri)
include(../SortedBinaryNewickVector/SortedBinaryNewickVector.pri)
include(../TwoDigitNewick/TwoDigitNewick.pri)

include(NewickDemoConsole.pri)
include(NewickDemoDesktop.pri)

SOURCES += \
    qtmain.cpp\

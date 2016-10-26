CONFIG += debug_and_release

include(../RibiLibraries/DesktopApplication.pri)
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

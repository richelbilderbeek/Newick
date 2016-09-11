include(../RibiLibraries/DesktopApplication.pri)
include(../RibiLibraries/Boost.pri)
include(../RibiLibraries/BigInteger.pri)
include(../RibiLibraries/GeneralConsole.pri)
include(../RibiLibraries/GeneralDesktop.pri)

include(../RibiClasses/CppContainer/CppContainer.pri)
include(../RibiClasses/CppFuzzy_equal_to/CppFuzzy_equal_to.pri)
include(../ManyDigitNewick/ManyDigitNewick.pri)
include(../RibiClasses/CppMultiVector/CppMultiVector.pri)
include(../BinaryNewickVector/BinaryNewickVector.pri)
include(Newick.pri)
include(../NewickVector/NewickVector.pri)
include(../RibiClasses/CppSortedBinaryNewickVector/CppSortedBinaryNewickVector.pri)
include(../TwoDigitNewick/TwoDigitNewick.pri)

include(TestNewickDesktop.pri)

SOURCES += \
    qtmain.cpp\

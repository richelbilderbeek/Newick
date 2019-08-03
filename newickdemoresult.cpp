
#include <iostream>

#pragma GCC diagnostic push



#include "newickdemoresult.h"
#pragma GCC diagnostic pop

std::ostream& operator<<(std::ostream& os, const TestNewickResult& r)
{
  os
    << r.newick << '\t'
    << r.theta << '\t'
    << r.test_name << '\t'
    << r.probability << '\t'
    << r.time;
  return os;
}



#ifndef NEWICKCPP98_H
#define NEWICKCPP98_H

#include <string>
#include <vector>

#pragma GCC diagnostic push

#include <boost/tuple/tuple.hpp>
#pragma GCC diagnostic pop

namespace ribi {

struct NewickCpp98
{
  NewickCpp98();

  ///Functions that do not use the C++11 standard
  std::vector<std::string> CreateValidTrinaryNewicks() noexcept;

  ///Functions that do not use the C++11 standard
  ///Thus, cannot use std::tuple
  std::vector<boost::tuple<std::string,double,double> > GetKnownProbabilities() noexcept;

  ///Functions that do not use the C++11 standard
  std::vector<std::pair<std::vector<int>,int> > GetSimplerNewicksFrequencyPairs(
    const std::vector<int>& n
  );

  ///Functions that do not use the C++11 standard
  std::vector<std::vector<int> > GetRootBranches(const std::vector<int>& n);
};

} //~namespace ribi

#endif // NEWICKCPP98_H

//---------------------------------------------------------------------------
/*
TestNewick, test the Newick classes and functions
Copyright (C) 2009-2015 Richel Bilderbeek

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.If not, see <http://www.gnu.org/licenses/>.
*/
//---------------------------------------------------------------------------
//From http://www.richelbilderbeek.nl/ToolTestNewick.htm
//---------------------------------------------------------------------------
#include <fstream>
#include <vector>





#include <boost/numeric/conversion/cast.hpp>
#include <boost/timer.hpp>

#include "binarynewickvector.h"
#include "manydigitnewick.h"
#include "newick.h"
#include "newickvector.h"
#include "sortedbinarynewickvector.h"
#include "twodigitnewick.h"
#include "newickdemo.h"



ribi::TestNewick::TestNewick()
  : m_time(0.0),
    m_probability(0.0)
{
  //Do nothing
}

///CreateAllTests creates all tests
const std::vector<boost::shared_ptr<ribi::TestNewick>> ribi::TestNewick::CreateTests(const int flags)
{
  std::vector<boost::shared_ptr<TestNewick> > v;
  if (flags & m_flag_binary_newick_vector)
    v.push_back(boost::shared_ptr<TestNewick>(new TestBinaryNewickVector));
  if (flags & m_flag_many_digit_newick)
    v.push_back(boost::shared_ptr<TestNewick>(new TestManyDigitNewick));
  if (flags & m_flag_newick_vector)
    v.push_back(boost::shared_ptr<TestNewick>(new TestNewickVector));
  if (flags & m_flag_sorted_binary_newick_vector)
    v.push_back(boost::shared_ptr<TestNewick>(new TestSortedBinaryNewickVector));
  if (flags & m_flag_two_digit_newick)
    v.push_back(boost::shared_ptr<TestNewick>(new TestTwoDigitNewick));
  return v;
}

void ribi::TestNewick::SetProbability(const double probability)
{
  m_probability = probability;
}

void ribi::TestNewick::SetTime(const double time)
{
  m_time = time;
}

bool ribi::TestBinaryNewickVector::CanCalculate(const std::string& newick_str, const double theta)
{
  if (theta <= 0.0) return false;
  if (!newick::IsNewick(newick_str))
    return false;
  const std::vector<int> newick = newick::StringToNewick(newick_str);
  if (!newick::IsUnaryNewick(newick) && !newick::IsBinaryNewick(newick))
    return false;
  return true;
}

void ribi::TestBinaryNewickVector::Calculate(const std::string& newick_str, const double theta)
{
  assert(CanCalculate(newick_str,theta));
  boost::timer t;
  const double p
    = BinaryNewickVector::CalculateProbability(
      newick_str,theta);
  SetTime(t.elapsed());
  SetProbability(p);
}

bool ribi::TestManyDigitNewick::CanCalculate(const std::string& newick_str, const double theta)
{
  //TestManyDigitNewick gives incorrect results!
  return false;

  if (theta <= 0.0) return false;
  if (!newick::IsNewick(newick_str)) return false;
  return true;
}

void ribi::TestManyDigitNewick::Calculate(const std::string& newick_str, const double theta)
{
  assert(CanCalculate(newick_str,theta));
  boost::timer t;
  const double p
    = ManyDigitNewick::CalculateProbability(
      newick_str,theta);
  SetTime(t.elapsed());
  SetProbability(p);
}

bool ribi::TestNewickVector::CanCalculate(const std::string& newick_str, const double theta)
{
  if (theta <= 0.0) return false;
  if (!newick::IsNewick(newick_str)) return false;
  return true;
}

void ribi::TestNewickVector::Calculate(const std::string& newick_str, const double theta)
{
  assert(CanCalculate(newick_str,theta));
  boost::timer t;
  const double p
    = CalculateProbabilityNewickVector(
      newick_str,theta);
  SetTime(t.elapsed());
  SetProbability(p);
}

bool ribi::TestSortedBinaryNewickVector::CanCalculate(const std::string& newick_str, const double theta)
{
  if (theta <= 0.0) return false;
  if (!newick::IsNewick(newick_str))
    return false;
  const std::vector<int> newick = newick::StringToNewick(newick_str);
  if (!newick::IsUnaryNewick(newick) && !newick::IsBinaryNewick(newick))
    return false;
  return true;
}

void ribi::TestSortedBinaryNewickVector::Calculate(const std::string& newick_str, const double theta)
{
  assert(CanCalculate(newick_str,theta));
  boost::timer t;
  const double p
    = SortedBinaryNewickVector::CalculateProbability(
      newick_str,
      theta);
  SetTime(t.elapsed());
  SetProbability(p);
}

bool ribi::TestTwoDigitNewick::CanCalculate(const std::string& newick_str, const double theta)
{
  if (theta <= 0.0) return false;
  if (!newick::IsNewick(newick_str))
    return false;
  const std::vector<int> newick = newick::StringToNewick(newick_str);
  if (!newick::IsUnaryNewick(newick) && !newick::IsBinaryNewick(newick))
    return false;
  return true;
}

void ribi::TestTwoDigitNewick::Calculate(const std::string& newick_str, const double theta)
{
  assert(CanCalculate(newick_str,theta));
  boost::timer t;
  const double p{
    CalculateProbabilityTwoDigitNewick(
      newick_str,
      theta
    )
  };
  SetTime(t.elapsed());
  SetProbability(p);
}


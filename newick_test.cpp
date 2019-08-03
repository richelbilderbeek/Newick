#pragma GCC diagnostic push

#include "newick.h"

#include <algorithm>
#include <deque>
#include <iostream>
#include <functional>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <sstream>
#include <tuple>
#include <vector>

#include "fuzzy_equal_to.h"
#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include "BigIntegerLibrary.hh"

#include "newickcpp98.h"
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/test/unit_test.hpp>

using namespace ribi;

BOOST_AUTO_TEST_CASE(ribi_newick_test_tests_with_newickcpp98)
{
  //Check difference between C++98 and C++0x
  BOOST_CHECK(ribi::newick::CreateValidTrinaryNewicks() == NewickCpp98().CreateValidTrinaryNewicks());
  BOOST_CHECK(ribi::newick::GetKnownProbabilities() == NewickCpp98().GetKnownProbabilities());
}

BOOST_AUTO_TEST_CASE(ribi_newick_str_to_vector_1)
{
  //Check conversions from std::string to std::vector #1
  const std::vector<int> v = ribi::newick::StringToNewick("(11,(22,33))");
  BOOST_CHECK(v.size() == 7);
  BOOST_CHECK(v[0]==ribi::newick::bracket_open);
  BOOST_CHECK(v[1]==11);
  BOOST_CHECK(v[2]==ribi::newick::bracket_open);
  BOOST_CHECK(v[3]==22);
  BOOST_CHECK(v[4]==33);
  BOOST_CHECK(v[5]==ribi::newick::bracket_close);
  BOOST_CHECK(v[6]==ribi::newick::bracket_close);
}

BOOST_AUTO_TEST_CASE(ribi_newick_str_to_vector_2)
{
  //Check conversions from std::string to std::vector #2
  const std::vector<int> v = ribi::newick::StringToNewick("((11,22),33)");
  BOOST_CHECK(v.size() == 7);
  BOOST_CHECK(v[0]==ribi::newick::bracket_open);
  BOOST_CHECK(v[1]==ribi::newick::bracket_open);
  BOOST_CHECK(v[2]==11);
  BOOST_CHECK(v[3]==22);
  BOOST_CHECK(v[4]==ribi::newick::bracket_close);
  BOOST_CHECK(v[5]==33);
  BOOST_CHECK(v[6]==ribi::newick::bracket_close);
}

BOOST_AUTO_TEST_CASE(ribi_newick_well_formed_newicks_must_be_accepted)
{
  //Check if well-formed Newicks are accepted
  const std::vector<std::string> v = ribi::newick::CreateValidNewicks();
  for(const std::string& s: v)
  {
    BOOST_CHECK(ribi::newick::IsNewick(s));
    const std::vector<int> v = ribi::newick::StringToNewick(s);
    BOOST_CHECK(ribi::newick::IsNewick(v));
  }
}

BOOST_AUTO_TEST_CASE(ribi_newick_ill_formed_newicks_must_be_rejected)
{
  //Check if ill-formed Newicks are rejected
  const std::vector<std::string> v = ribi::newick::CreateInvalidNewicks();
  for(const std::string& s: v)
  {
    #ifdef TRACE_REJECTED_NEWICKS
    const std::string debug = "I must be rejected: " + s;
    TRACE(debug);
    #endif
    BOOST_CHECK(!ribi::newick::IsNewick(s));
    //Cannot test if std::vector<int> versions are rejected,
    //because ribi::newick::StringToNewick assumes a valid Newick
    //const std::vector<int> v = ribi::newick::StringToNewick(s);
    //BOOST_CHECK(!ribi::newick::IsNewick(v));
  }
}

BOOST_AUTO_TEST_CASE(ribi_newick_CalcNumOfSymmetriesBinary)
{
  using namespace ribi::newick;
  BOOST_CHECK(CalcNumOfSymmetriesBinary(StringToNewick("(1,(3,1))"))==0);
  BOOST_CHECK(CalcNumOfSymmetriesBinary(StringToNewick("(3,(1,1))"))==1);
  BOOST_CHECK(CalcNumOfSymmetriesBinary(StringToNewick("(1,((1,1),(1,1)))"))==3);
  BOOST_CHECK(CalcNumOfSymmetriesBinary(StringToNewick("(1,((1,1),(2,2)))"))==2);
  BOOST_CHECK(CalcNumOfSymmetriesBinary(StringToNewick("(1,(2,3))"))==0);
  BOOST_CHECK(CalcNumOfSymmetriesBinary(StringToNewick("(99,99)"))==1);
  BOOST_CHECK(CalcNumOfSymmetriesBinary(StringToNewick("(3,(2,2))"))==1);
  BOOST_CHECK(CalcNumOfSymmetriesBinary(StringToNewick("(2,(2,2))"))==1);
  BOOST_CHECK(CalcNumOfSymmetriesBinary(StringToNewick("((3,3),(2,2))"))==2);
  BOOST_CHECK(CalcNumOfSymmetriesBinary(StringToNewick("((3,3),(3,3))"))==3);
  BOOST_CHECK(CalcNumOfSymmetriesBinary(StringToNewick("((3,3),(3,4))"))==1);
  BOOST_CHECK(CalcNumOfSymmetriesBinary(StringToNewick("(((3,3),(4,4)),5)"))==2);
  BOOST_CHECK(CalcNumOfSymmetriesBinary(StringToNewick("(((3,3),(5,5)),5)"))==2);
  BOOST_CHECK(CalcNumOfSymmetriesBinary(StringToNewick("(((5,5),(5,5)),5)"))==3);
  BOOST_CHECK(CalcNumOfSymmetriesBinary(StringToNewick("(((5,5),(5,5)),(4,4))"))==4);
  BOOST_CHECK(CalcNumOfSymmetriesBinary(StringToNewick("(((5,5),(4,4)),(4,4))"))==3);
  BOOST_CHECK(CalcNumOfSymmetriesBinary(StringToNewick("(((4,4),(4,4)),(4,4))"))==4);
  BOOST_CHECK(CalcNumOfCombinationsBinary(StringToNewick("(3,(1,1))"))==10);
  BOOST_CHECK(CalcNumOfCombinationsBinary(StringToNewick("(1,(3,1))"))==20);
  BOOST_CHECK(CalcNumOfCombinationsBinary(StringToNewick("(1,(1,(1,(1,1))))"))==60);
  BOOST_CHECK(CalcNumOfCombinationsBinary(StringToNewick("(1,((1,1),(1,1)))"))==15);
}

BOOST_AUTO_TEST_CASE(ribi_newick_CreateSortedPair)
{
  //Supplied in order
  {
    const auto p = ribi::newick::CreateSortedPair(3,14);
    BOOST_CHECK_EQUAL(p.first, 3);
    BOOST_CHECK_EQUAL(p.second, 14);
  }
  //Supplied in reverse order
  {
    const auto p = ribi::newick::CreateSortedPair(14,3);
    BOOST_CHECK_EQUAL(p.first, 3);
    BOOST_CHECK_EQUAL(p.second, 14);
  }
}

BOOST_AUTO_TEST_CASE(ribi_newick_bigIntegerToString)
{
  BOOST_CHECK(bigIntegerToString(ribi::newick::FactorialBigInt(1))=="1");
  BOOST_CHECK(bigIntegerToString(ribi::newick::FactorialBigInt(2))=="2");
  BOOST_CHECK(bigIntegerToString(ribi::newick::FactorialBigInt(3))=="6");
  BOOST_CHECK(bigIntegerToString(ribi::newick::FactorialBigInt(4))=="24");
  BOOST_CHECK(bigIntegerToString(ribi::newick::FactorialBigInt(5))=="120");
  BOOST_CHECK(bigIntegerToString(ribi::newick::FactorialBigInt(6))=="720");

}

BOOST_AUTO_TEST_CASE(ribi_newick_GetLeafMaxArity)
{
  BOOST_CHECK(ribi::newick::GetLeafMaxArity(ribi::newick::StringToNewick("(1)"))   == 1);
  BOOST_CHECK(ribi::newick::GetLeafMaxArity(ribi::newick::StringToNewick("(12)"))  == 1);
  BOOST_CHECK(ribi::newick::GetLeafMaxArity(ribi::newick::StringToNewick("(123)")) == 1);
  BOOST_CHECK(ribi::newick::GetLeafMaxArity(ribi::newick::StringToNewick("(1,2)"))   == 2);
  BOOST_CHECK(ribi::newick::GetLeafMaxArity(ribi::newick::StringToNewick("(12,2)"))  == 2);
  BOOST_CHECK(ribi::newick::GetLeafMaxArity(ribi::newick::StringToNewick("(123,2)")) == 2);
  BOOST_CHECK(ribi::newick::GetLeafMaxArity(ribi::newick::StringToNewick("(1,(1,2))"))   == 2);
  BOOST_CHECK(ribi::newick::GetLeafMaxArity(ribi::newick::StringToNewick("(1,(12,2))"))  == 2);
  BOOST_CHECK(ribi::newick::GetLeafMaxArity(ribi::newick::StringToNewick("(1,(123,2))")) == 2);
  BOOST_CHECK(ribi::newick::GetLeafMaxArity(ribi::newick::StringToNewick("((1,2),3)"))   == 2);
  BOOST_CHECK(ribi::newick::GetLeafMaxArity(ribi::newick::StringToNewick("((12,2),3)"))  == 2);
  BOOST_CHECK(ribi::newick::GetLeafMaxArity(ribi::newick::StringToNewick("((123,2),3)")) == 2);
  BOOST_CHECK(ribi::newick::GetLeafMaxArity(ribi::newick::StringToNewick("(1,2,3)"))   == 3);
  BOOST_CHECK(ribi::newick::GetLeafMaxArity(ribi::newick::StringToNewick("(12,2,3)"))  == 3);
  BOOST_CHECK(ribi::newick::GetLeafMaxArity(ribi::newick::StringToNewick("(123,2,3)")) == 3);
  BOOST_CHECK(ribi::newick::GetLeafMaxArity(ribi::newick::StringToNewick("(1,(1,2,3))"))   == 3);
  BOOST_CHECK(ribi::newick::GetLeafMaxArity(ribi::newick::StringToNewick("(1,(12,2,3))"))  == 3);
  BOOST_CHECK(ribi::newick::GetLeafMaxArity(ribi::newick::StringToNewick("(1,(123,2,3))")) == 3);
  BOOST_CHECK(ribi::newick::GetLeafMaxArity(ribi::newick::StringToNewick("((1,2,3),4)"))   == 3);
  BOOST_CHECK(ribi::newick::GetLeafMaxArity(ribi::newick::StringToNewick("((12,2,3),4)"))  == 3);
  BOOST_CHECK(ribi::newick::GetLeafMaxArity(ribi::newick::StringToNewick("((123,2,3),4)")) == 3);
}

BOOST_AUTO_TEST_CASE(ribi_newick_CalcDenominator)
{
  ribi::fuzzy_equal_to f;
  BOOST_CHECK(f(  2.0,ribi::newick::CalcDenominator(ribi::newick::StringToNewick("(1,1)"),10.0)));
  BOOST_CHECK(f(  6.0,ribi::newick::CalcDenominator(ribi::newick::StringToNewick("((1,1),1)"),10.0)));
  BOOST_CHECK(f( 26.0,ribi::newick::CalcDenominator(ribi::newick::StringToNewick("(1,2)"),10.0)));
  BOOST_CHECK(f( 32.0,ribi::newick::CalcDenominator(ribi::newick::StringToNewick("((1,1),2)"),10.0)));
  BOOST_CHECK(f( 32.0,ribi::newick::CalcDenominator(ribi::newick::StringToNewick("(2,(1,1))"),10.0)));
  BOOST_CHECK(f( 50.0,ribi::newick::CalcDenominator(ribi::newick::StringToNewick("((1,1),3)"),10.0)));
  BOOST_CHECK(f( 80.0,ribi::newick::CalcDenominator(ribi::newick::StringToNewick("((1,2),3)"),10.0)));
  BOOST_CHECK(f( 80.0,ribi::newick::CalcDenominator(ribi::newick::StringToNewick("((3,1),2)"),10.0)));
  BOOST_CHECK(f( 80.0,ribi::newick::CalcDenominator(ribi::newick::StringToNewick("((2,3),1)"),10.0)));
  BOOST_CHECK(f(102.0,ribi::newick::CalcDenominator(ribi::newick::StringToNewick("((2,1),4)"),10.0)));
  BOOST_CHECK(f(152.0,ribi::newick::CalcDenominator(ribi::newick::StringToNewick("(2,(1,(3,3)))"),10.0)));
  BOOST_CHECK(f(162.0,ribi::newick::CalcDenominator(ribi::newick::StringToNewick("((2,3),4)"),10.0)));
  BOOST_CHECK(f(180.0,ribi::newick::CalcDenominator(ribi::newick::StringToNewick("((1,2),(3,4))"),10.0)));
  BOOST_CHECK(f(180.0,ribi::newick::CalcDenominator(ribi::newick::StringToNewick("((4,1),(2,3))"),10.0)));
  BOOST_CHECK(f(180.0,ribi::newick::CalcDenominator(ribi::newick::StringToNewick("((3,4),(1,2))"),10.0)));
  BOOST_CHECK(f(180.0,ribi::newick::CalcDenominator(ribi::newick::StringToNewick("((2,3),(4,1))"),10.0)));
}

BOOST_AUTO_TEST_CASE(ribi_newick_FindPosAfter_and_FindPosBefore)
{
  const std::vector<int> v = { 0,1,2,3,4,5,6 };
  BOOST_CHECK(ribi::newick::FindPosAfter(v,3,2)==3);
  BOOST_CHECK(ribi::newick::FindPosAfter(v,4,2)==4);
  BOOST_CHECK(ribi::newick::FindPosAfter(v,5,2)==5);
  BOOST_CHECK(ribi::newick::FindPosAfter(v,6,2)==6);
  BOOST_CHECK(ribi::newick::FindPosBefore(v,3,4)==3);
  BOOST_CHECK(ribi::newick::FindPosBefore(v,2,4)==2);
  BOOST_CHECK(ribi::newick::FindPosBefore(v,1,4)==1);
  BOOST_CHECK(ribi::newick::FindPosBefore(v,0,4)==0);
}

BOOST_AUTO_TEST_CASE(ribi_newick_GetDepth)
{
  {
    const std::vector<int> v = ribi::newick::GetDepth(ribi::newick::StringToNewick("(1,(2,2))"));
    const std::vector<int> w = ribi::newick::GetDepth(ribi::newick::StringToNewick("(9,(9,9))"));
    const std::vector<int> x = { 0,0,1,1,1,1,0 };
    BOOST_CHECK(v == x);
    BOOST_CHECK(w == x);
  }
  {
    const std::vector<int> v = ribi::newick::GetDepth(ribi::newick::StringToNewick("((2,2),1)"));
    const std::vector<int> w = { 0,1,1,1,1,0,0 };
    BOOST_CHECK(v == w);
  }
  {
    const std::vector<int> v = ribi::newick::GetDepth(ribi::newick::StringToNewick("(1,(2,2),1)"));
    const std::vector<int> w = { 0,0,1,1,1,1,0,0 };
    BOOST_CHECK(v == w);
  }
  {
    const std::vector<int> v = ribi::newick::GetDepth(ribi::newick::StringToNewick("(1,(2,3),4,(5,6))"));
    const std::vector<int> w = { 0,0,1,1,1,1,0,1,1,1,1,0 };
    BOOST_CHECK(v == w);
  }
  {
    const std::vector<int> v = ribi::newick::GetDepth(ribi::newick::StringToNewick("(1,(2,3),(5,6))"));
    const std::vector<int> w = { 0,0,1,1,1,1,1,1,1,1,0 };
    BOOST_CHECK(v == w);
  }
  {
    const std::vector<int> v = ribi::newick::GetDepth(ribi::newick::StringToNewick("(1,(2,(3,4)),((5,6),7))"));
    const std::vector<int> w = { 0,0,1,1,2,2,2,2,1,1,2,2,2,2,1,1,0 };
    BOOST_CHECK(v == w);
  }
}

BOOST_AUTO_TEST_CASE(ribi_newick_GetRootBranches)
{
  {
    const std::vector<std::vector<int> > v = ribi::newick::GetRootBranches(ribi::newick::StringToNewick("(1,2)"));
    BOOST_CHECK(v.size() == 2);
    BOOST_CHECK(std::find(v.begin(),v.end(),
      ribi::newick::StringToNewick("(1)")) != v.end());
    BOOST_CHECK(std::find(v.begin(),v.end(),
      ribi::newick::StringToNewick("(2)")) != v.end());
  }
  {
    const std::vector<std::vector<int> > v = ribi::newick::GetRootBranches(ribi::newick::StringToNewick("(1,(2,3))"));
    BOOST_CHECK(v.size() == 2);
    BOOST_CHECK(std::find(v.begin(),v.end(),
      ribi::newick::StringToNewick("(1)")) != v.end());
    BOOST_CHECK(std::find(v.begin(),v.end(),
      ribi::newick::StringToNewick("(2,3)")) != v.end());
  }
  {
    const std::vector<std::vector<int> > v = ribi::newick::GetRootBranches(ribi::newick::StringToNewick("(1,2,(3,4))"));
    BOOST_CHECK(v.size() == 3);
    BOOST_CHECK(std::find(v.begin(),v.end(),
      ribi::newick::StringToNewick("(1)")) != v.end());
    BOOST_CHECK(std::find(v.begin(),v.end(),
      ribi::newick::StringToNewick("(2)")) != v.end());
    BOOST_CHECK(std::find(v.begin(),v.end(),
      ribi::newick::StringToNewick("(3,4)")) != v.end());
  }
}

BOOST_AUTO_TEST_CASE(ribi_newick_GetRootBranches_cpp98_and_cpp11_must_have_same)
{
  //Compare C++98 and C++0x version
  const std::vector<std::string> v = ribi::newick::CreateValidBinaryNewicks();
  for(const std::string& s: v)
  {
    const std::vector<int> n = ribi::newick::StringToNewick(s);
    BOOST_CHECK(ribi::newick::GetRootBranches(n) == NewickCpp98().GetRootBranches(n));
  }
}

BOOST_AUTO_TEST_CASE(ribi_newick_IsUnaryNewick_on_CreateValidUnaryNewicks)
{
  //Check if unary Newicks are detected correctly
  {
    const std::vector<std::string> v = ribi::newick::CreateValidUnaryNewicks();
    for(const std::string& s: v)
    {
      const std::vector<int> n = ribi::newick::StringToNewick(s);
      BOOST_CHECK( ribi::newick::GetLeafMaxArity(n)<=1);
      BOOST_CHECK( ribi::newick::IsUnaryNewick(n));
      BOOST_CHECK(!ribi::newick::IsBinaryNewick(n));
      BOOST_CHECK(!ribi::newick::IsTrinaryNewick(n));
    }
  }
}

BOOST_AUTO_TEST_CASE(ribi_newick_IsBinaryNewick_on_CreateValidBinaryNewicks)
{
  //Check if binary Newicks are detected correctly
  {
    const std::vector<std::string> v = ribi::newick::CreateValidBinaryNewicks();
    for(const std::string& s: v)
    {
      const std::vector<int> n = ribi::newick::StringToNewick(s);
      BOOST_CHECK( ribi::newick::GetLeafMaxArity(n)<=2);
      BOOST_CHECK(!ribi::newick::IsUnaryNewick(n));
      BOOST_CHECK( ribi::newick::IsBinaryNewick(n));
      BOOST_CHECK(!ribi::newick::IsTrinaryNewick(n));
    }
  }
}

BOOST_AUTO_TEST_CASE(ribi_newick_IsTrinaryNewick_on_CreateValidTrinaryNewicks)
{
  //Check if trinary Newicks are detected correctly
  {
    const std::vector<std::string> v = ribi::newick::CreateValidTrinaryNewicks();
    for(const std::string& s: v)
    {
      //TRACE(s);
      const std::vector<int> n = ribi::newick::StringToNewick(s);
      BOOST_CHECK( ribi::newick::GetLeafMaxArity(n)<=3);
      BOOST_CHECK(!ribi::newick::IsUnaryNewick(n));
      BOOST_CHECK(!ribi::newick::IsBinaryNewick(n));
      BOOST_CHECK( ribi::newick::IsTrinaryNewick(n));
    }
  }
}

BOOST_AUTO_TEST_CASE(ribi_newick_GetSimplerNewicks)
{
  {
    const std::string s("(1,(2,3))");
    const std::vector<std::vector<int>> n{
      newick::GetSimplerNewicks(newick::StringToNewick(s))
    };
    BOOST_CHECK(n.size() == 2);
    BOOST_CHECK(std::find(n.begin(),n.end(),newick::StringToNewick("(1,(1,3))"))
      != n.end());
    BOOST_CHECK(std::find(n.begin(),n.end(),newick::StringToNewick("(1,(2,2))"))
      != n.end());
  }
  {
    const std::string s("(1,(2,3,4))");
    const std::vector<std::vector<int>> n{
      newick::GetSimplerNewicks(newick::StringToNewick(s))
    };
    BOOST_CHECK(n.size() == 3);
    BOOST_CHECK(std::find(n.begin(),n.end(),newick::StringToNewick("(1,(1,3,4))"))
      != n.end());
    BOOST_CHECK(std::find(n.begin(),n.end(),newick::StringToNewick("(1,(2,2,4))"))
      != n.end());
    BOOST_CHECK(std::find(n.begin(),n.end(),newick::StringToNewick("(1,(2,3,3))"))
      != n.end());
  }
  {
    const std::string s("(1,(1,3,4))");
    const std::vector<std::vector<int> > n{
      newick::GetSimplerNewicks(newick::StringToNewick(s))
    };
    BOOST_CHECK(n.size() == 4);
    BOOST_CHECK(std::find(n.begin(),n.end(),newick::StringToNewick("(1,(4,4))"))
      != n.end());
    BOOST_CHECK(std::find(n.begin(),n.end(),newick::StringToNewick("(1,(3,5))"))
      != n.end());
    BOOST_CHECK(std::find(n.begin(),n.end(),newick::StringToNewick("(1,(1,2,4))"))
      != n.end());
    BOOST_CHECK(std::find(n.begin(),n.end(),newick::StringToNewick("(1,(1,3,3))"))
      != n.end());
  }
}

BOOST_AUTO_TEST_CASE(ribi_newick_GetSimplerNewicksFrequencyPairs)
{
  {
    const std::string s("(1,(1,3,4))");
    const std::vector<std::pair<std::vector<int>,int> > n
      = ribi::newick::GetSimplerNewicksFrequencyPairs(ribi::newick::StringToNewick(s));
    //typedef std::pair<std::vector<int>,int> Pair;
    BOOST_CHECK(n.size() == 4);
    BOOST_CHECK(std::find(n.begin(),n.end(),
      std::make_pair(ribi::newick::StringToNewick("(1,(4,4))"),1))
      != n.end());
    BOOST_CHECK(std::find(n.begin(),n.end(),
      std::make_pair(ribi::newick::StringToNewick("(1,(3,5))"),1))
      != n.end());
    BOOST_CHECK(std::find(n.begin(),n.end(),
      std::make_pair(ribi::newick::StringToNewick("(1,(1,2,4))"),3))
      != n.end());
    BOOST_CHECK(std::find(n.begin(),n.end(),
      std::make_pair(ribi::newick::StringToNewick("(1,(1,3,3))"),4))
      != n.end());
  }
}

BOOST_AUTO_TEST_CASE(ribi_newick_GetSimplerNewicks_2)
{
  const std::string s("((1,1),2)");
  const std::vector<std::vector<int>> n = ribi::newick::GetSimplerNewicks(
    ribi::newick::StringToNewick(s));
  BOOST_CHECK(n.size() == 3);
  BOOST_CHECK(std::find(n.begin(),n.end(),
    ribi::newick::StringToNewick("(2,2)"))
    != n.end());
  BOOST_CHECK(std::find(n.begin(),n.end(),
    ribi::newick::StringToNewick("((1,1),1)"))
    != n.end());
}

BOOST_AUTO_TEST_CASE(ribi_newick_GetSimplerNewicksFrequencyPairs_2)
{
  const std::string s("((1,1),2)");
  typedef std::pair<std::vector<int>,int> Pair;
  const std::vector<Pair> n
    = ribi::newick::GetSimplerNewicksFrequencyPairs(ribi::newick::StringToNewick(s));
  BOOST_CHECK(n.size() == 3);
  BOOST_CHECK(std::find(n.begin(),n.end(),
    std::make_pair(ribi::newick::StringToNewick("(2,2)"),1))
    != n.end());
  BOOST_CHECK(std::find(n.begin(),n.end(),
    std::make_pair(ribi::newick::StringToNewick("((1,1),1)"),2))
    != n.end());
}

BOOST_AUTO_TEST_CASE(ribi_newick_GetSimplerNewicks_3)
{
  const std::string s("((2,1),4)");
  const std::vector<std::vector<int>> n = ribi::newick::GetSimplerNewicks(
    ribi::newick::StringToNewick(s));
  BOOST_CHECK(n.size() == 3);
  BOOST_CHECK(std::find(n.begin(),n.end(),
    ribi::newick::StringToNewick("(3,4)"))
    != n.end());
  BOOST_CHECK(std::find(n.begin(),n.end(),
    ribi::newick::StringToNewick("((1,1),4)"))
    != n.end());
  BOOST_CHECK(std::find(n.begin(),n.end(),
    ribi::newick::StringToNewick("((2,1),3)"))
    != n.end());
}

BOOST_AUTO_TEST_CASE(ribi_newick_GetSimplerNewicksFrequencyPairs_3)
{
  const std::string s("((2,1),4)");
  typedef std::pair<std::vector<int>,int> Pair;
  const std::vector<Pair> n
    = ribi::newick::GetSimplerNewicksFrequencyPairs(ribi::newick::StringToNewick(s));
  BOOST_CHECK(n.size() == 3);
  BOOST_CHECK(std::find(n.begin(),n.end(),
    std::make_pair(ribi::newick::StringToNewick("(3,4)"),1))
    != n.end());
  BOOST_CHECK(std::find(n.begin(),n.end(),
    std::make_pair(ribi::newick::StringToNewick("((1,1),4)"),2))
    != n.end());
  BOOST_CHECK(std::find(n.begin(),n.end(),
    std::make_pair(ribi::newick::StringToNewick("((2,1),3)"),4))
    != n.end());
}


BOOST_AUTO_TEST_CASE(ribi_newick_GetSimplerNewicks_4)
{
  const std::string s("((2,3),4)");
  const std::vector<std::vector<int> > n = ribi::newick::GetSimplerNewicks(
    ribi::newick::StringToNewick(s));
  BOOST_CHECK(n.size() == 3);
  BOOST_CHECK(std::find(n.begin(),n.end(),
    ribi::newick::StringToNewick("((1,3),4)"))
    != n.end());
  BOOST_CHECK(std::find(n.begin(),n.end(),
    ribi::newick::StringToNewick("((2,2),4)"))
    != n.end());
  BOOST_CHECK(std::find(n.begin(),n.end(),
    ribi::newick::StringToNewick("((2,3),3)"))
    != n.end());
}

BOOST_AUTO_TEST_CASE(ribi_newick_GetSimplerNewicksFrequencyPairs_4)
{
  const std::string s("((2,3),4)");
  typedef std::pair<std::vector<int>,int> Pair;
  const std::vector<Pair> n
    = ribi::newick::GetSimplerNewicksFrequencyPairs(ribi::newick::StringToNewick(s));
  BOOST_CHECK(n.size() == 3);
  BOOST_CHECK(std::find(n.begin(),n.end(),
    std::make_pair(ribi::newick::StringToNewick("((1,3),4)"),2))
    != n.end());
  BOOST_CHECK(std::find(n.begin(),n.end(),
    std::make_pair(ribi::newick::StringToNewick("((2,2),4)"),3))
    != n.end());
  BOOST_CHECK(std::find(n.begin(),n.end(),
    std::make_pair(ribi::newick::StringToNewick("((2,3),3)"),4))
    != n.end());
}

BOOST_AUTO_TEST_CASE(ribi_newick_compare_GetSimplerNewicks_and_GetSimplerNewicksFrequencyPairs)
{
  //Compare GetSimplerNewicks and
  //GetSimplerNewicksFrequencyPairs
  const std::vector<std::string> newicks
    = ribi::newick::CreateValidNewicks();
  for(const std::string& newick_str: newicks)
  {
    const std::vector<int> newick
      = ribi::newick::StringToNewick(newick_str);
    const std::vector<std::vector<int> > v1
      = ribi::newick::GetSimplerNewicks(newick);
    const std::vector<std::pair<std::vector<int>,int> > v2
      = ribi::newick::GetSimplerNewicksFrequencyPairs(newick);
    BOOST_CHECK(v1.size() == v2.size());

    //Create sets that should be equal
    std::set<std::vector<int>> s1;
    std::copy(std::begin(v1), std::end(v1),
      std::inserter(s1, std::end(s1))
    );
    std::set<std::vector<int>> s2;
    std::transform(std::begin(v2), std::end(v2),
      std::inserter(s2, std::end(s2)),
      [](const auto& p) { return p.first;}
    );
    BOOST_CHECK(ribi::newick::GetSimplerNewicksFrequencyPairs(newick)
      == NewickCpp98().GetSimplerNewicksFrequencyPairs(newick));
  }
}

#include "newickcpp98.h"

#include <boost/numeric/conversion/cast.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

#include "fuzzy_equal_to.h"
#include "newick.h"

#include <set>
#include <boost/test/unit_test.hpp>

using namespace ribi;

BOOST_AUTO_TEST_CASE(ribi_newickcpp98_all)
{
  //Check difference between C++98 and C++0x
  BOOST_CHECK(newick::CreateValidTrinaryNewicks() == NewickCpp98().CreateValidTrinaryNewicks());
  BOOST_CHECK(newick::GetKnownProbabilities() == NewickCpp98().GetKnownProbabilities());

  //Check conversions from std::string to std::vector #1
  {
    const std::vector<int> v = newick::StringToNewick("(11,(22,33))");
    BOOST_CHECK(v.size() == 7);
    BOOST_CHECK(v[0]==newick::bracket_open);
    BOOST_CHECK(v[1]==11);
    BOOST_CHECK(v[2]==newick::bracket_open);
    BOOST_CHECK(v[3]==22);
    BOOST_CHECK(v[4]==33);
    BOOST_CHECK(v[5]==newick::bracket_close);
    BOOST_CHECK(v[6]==newick::bracket_close);
  }
  //Check if well-formed Newicks are accepted
  {
    const std::vector<std::string> v = newick::CreateValidNewicks();
    for(const std::string& s: v)
    {
      BOOST_CHECK(newick::IsNewick(s));
      const std::vector<int> v = newick::StringToNewick(s);
      BOOST_CHECK(newick::IsNewick(v));
    }
  }

  //Check if ill-formed Newicks are rejected
  {
    #ifndef NDEBUG
    const std::vector<std::string> v = newick::CreateInvalidNewicks();
    for(const std::string& s: v)
    {
      #ifdef TRACE_REJECTED_NEWICKS
      const std::string debug = "I must be rejected: " + s;
      TRACE(debug);
      #endif
      BOOST_CHECK(!newick::IsNewick(s));
      //Cannot test if std::vector<int> versions are rejected,
      //because newick::StringToNewick assumes a valid Newick
      //const std::vector<int> v = newick::StringToNewick(s);
      //BOOST_CHECK(!newick::IsNewick(v));
    }
    #endif
  }
  //Check conversions from std::string to std::vector #2
  {
    const std::vector<int> v = newick::StringToNewick("((11,22),33)");
    BOOST_CHECK(v.size() == 7);
    BOOST_CHECK(v[0]==newick::bracket_open);
    BOOST_CHECK(v[1]==newick::bracket_open);
    BOOST_CHECK(v[2]==11);
    BOOST_CHECK(v[3]==22);
    BOOST_CHECK(v[4]==newick::bracket_close);
    BOOST_CHECK(v[5]==33);
    BOOST_CHECK(v[6]==newick::bracket_close);
  }
  BOOST_CHECK(newick::CalcNumOfSymmetriesBinary(newick::StringToNewick("(1,(3,1))"))==0);
  BOOST_CHECK(newick::CalcNumOfSymmetriesBinary(newick::StringToNewick("(3,(1,1))"))==1);
  BOOST_CHECK(newick::CalcNumOfSymmetriesBinary(newick::StringToNewick("(1,((1,1),(1,1)))"))==3);
  BOOST_CHECK(newick::CalcNumOfSymmetriesBinary(newick::StringToNewick("(1,((1,1),(2,2)))"))==2);
  BOOST_CHECK(newick::CalcNumOfSymmetriesBinary(newick::StringToNewick("(1,(2,3))"))==0);
  BOOST_CHECK(newick::CalcNumOfSymmetriesBinary(newick::StringToNewick("(99,99)"))==1);
  BOOST_CHECK(newick::CalcNumOfSymmetriesBinary(newick::StringToNewick("(3,(2,2))"))==1);
  BOOST_CHECK(newick::CalcNumOfSymmetriesBinary(newick::StringToNewick("(2,(2,2))"))==1);
  BOOST_CHECK(newick::CalcNumOfSymmetriesBinary(newick::StringToNewick("((3,3),(2,2))"))==2);
  BOOST_CHECK(newick::CalcNumOfSymmetriesBinary(newick::StringToNewick("((3,3),(3,3))"))==3);
  BOOST_CHECK(newick::CalcNumOfSymmetriesBinary(newick::StringToNewick("((3,3),(3,4))"))==1);
  BOOST_CHECK(newick::CalcNumOfSymmetriesBinary(newick::StringToNewick("(((3,3),(4,4)),5)"))==2);
  BOOST_CHECK(newick::CalcNumOfSymmetriesBinary(newick::StringToNewick("(((3,3),(5,5)),5)"))==2);
  BOOST_CHECK(newick::CalcNumOfSymmetriesBinary(newick::StringToNewick("(((5,5),(5,5)),5)"))==3);
  BOOST_CHECK(newick::CalcNumOfSymmetriesBinary(newick::StringToNewick("(((5,5),(5,5)),(4,4))"))==4);
  BOOST_CHECK(newick::CalcNumOfSymmetriesBinary(newick::StringToNewick("(((5,5),(4,4)),(4,4))"))==3);
  BOOST_CHECK(newick::CalcNumOfSymmetriesBinary(newick::StringToNewick("(((4,4),(4,4)),(4,4))"))==4);
  BOOST_CHECK(newick::CalcNumOfCombinationsBinary(newick::StringToNewick("(3,(1,1))"))==10);
  BOOST_CHECK(newick::CalcNumOfCombinationsBinary(newick::StringToNewick("(1,(3,1))"))==20);
  BOOST_CHECK(newick::CalcNumOfCombinationsBinary(newick::StringToNewick("(1,(1,(1,(1,1))))"))==60);
  BOOST_CHECK(newick::CalcNumOfCombinationsBinary(newick::StringToNewick("(1,((1,1),(1,1)))"))==15);
  BOOST_CHECK(bigIntegerToString(newick::FactorialBigInt(1))=="1");
  BOOST_CHECK(bigIntegerToString(newick::FactorialBigInt(2))=="2");
  BOOST_CHECK(bigIntegerToString(newick::FactorialBigInt(3))=="6");
  BOOST_CHECK(bigIntegerToString(newick::FactorialBigInt(4))=="24");
  BOOST_CHECK(bigIntegerToString(newick::FactorialBigInt(5))=="120");
  BOOST_CHECK(bigIntegerToString(newick::FactorialBigInt(6))=="720");
  BOOST_CHECK(newick::GetLeafMaxArity(newick::StringToNewick("(1)"))   == 1);
  BOOST_CHECK(newick::GetLeafMaxArity(newick::StringToNewick("(12)"))  == 1);
  BOOST_CHECK(newick::GetLeafMaxArity(newick::StringToNewick("(123)")) == 1);
  BOOST_CHECK(newick::GetLeafMaxArity(newick::StringToNewick("(1,2)"))   == 2);
  BOOST_CHECK(newick::GetLeafMaxArity(newick::StringToNewick("(12,2)"))  == 2);
  BOOST_CHECK(newick::GetLeafMaxArity(newick::StringToNewick("(123,2)")) == 2);
  BOOST_CHECK(newick::GetLeafMaxArity(newick::StringToNewick("(1,(1,2))"))   == 2);
  BOOST_CHECK(newick::GetLeafMaxArity(newick::StringToNewick("(1,(12,2))"))  == 2);
  BOOST_CHECK(newick::GetLeafMaxArity(newick::StringToNewick("(1,(123,2))")) == 2);
  BOOST_CHECK(newick::GetLeafMaxArity(newick::StringToNewick("((1,2),3)"))   == 2);
  BOOST_CHECK(newick::GetLeafMaxArity(newick::StringToNewick("((12,2),3)"))  == 2);
  BOOST_CHECK(newick::GetLeafMaxArity(newick::StringToNewick("((123,2),3)")) == 2);
  BOOST_CHECK(newick::GetLeafMaxArity(newick::StringToNewick("(1,2,3)"))   == 3);
  BOOST_CHECK(newick::GetLeafMaxArity(newick::StringToNewick("(12,2,3)"))  == 3);
  BOOST_CHECK(newick::GetLeafMaxArity(newick::StringToNewick("(123,2,3)")) == 3);
  BOOST_CHECK(newick::GetLeafMaxArity(newick::StringToNewick("(1,(1,2,3))"))   == 3);
  BOOST_CHECK(newick::GetLeafMaxArity(newick::StringToNewick("(1,(12,2,3))"))  == 3);
  BOOST_CHECK(newick::GetLeafMaxArity(newick::StringToNewick("(1,(123,2,3))")) == 3);
  BOOST_CHECK(newick::GetLeafMaxArity(newick::StringToNewick("((1,2,3),4)"))   == 3);
  BOOST_CHECK(newick::GetLeafMaxArity(newick::StringToNewick("((12,2,3),4)"))  == 3);
  BOOST_CHECK(newick::GetLeafMaxArity(newick::StringToNewick("((123,2,3),4)")) == 3);

  BOOST_CHECK(fuzzy_equal_to()(  2.0,newick::CalcDenominator(newick::StringToNewick("(1,1)"),10.0)));
  BOOST_CHECK(fuzzy_equal_to()(  6.0,newick::CalcDenominator(newick::StringToNewick("((1,1),1)"),10.0)));
  BOOST_CHECK(fuzzy_equal_to()( 26.0,newick::CalcDenominator(newick::StringToNewick("(1,2)"),10.0)));
  BOOST_CHECK(fuzzy_equal_to()( 32.0,newick::CalcDenominator(newick::StringToNewick("((1,1),2)"),10.0)));
  BOOST_CHECK(fuzzy_equal_to()( 32.0,newick::CalcDenominator(newick::StringToNewick("(2,(1,1))"),10.0)));
  BOOST_CHECK(fuzzy_equal_to()( 50.0,newick::CalcDenominator(newick::StringToNewick("((1,1),3)"),10.0)));
  BOOST_CHECK(fuzzy_equal_to()( 80.0,newick::CalcDenominator(newick::StringToNewick("((1,2),3)"),10.0)));
  BOOST_CHECK(fuzzy_equal_to()( 80.0,newick::CalcDenominator(newick::StringToNewick("((3,1),2)"),10.0)));
  BOOST_CHECK(fuzzy_equal_to()( 80.0,newick::CalcDenominator(newick::StringToNewick("((2,3),1)"),10.0)));
  BOOST_CHECK(fuzzy_equal_to()(102.0,newick::CalcDenominator(newick::StringToNewick("((2,1),4)"),10.0)));
  BOOST_CHECK(fuzzy_equal_to()(152.0,newick::CalcDenominator(newick::StringToNewick("(2,(1,(3,3)))"),10.0)));
  BOOST_CHECK(fuzzy_equal_to()(162.0,newick::CalcDenominator(newick::StringToNewick("((2,3),4)"),10.0)));
  BOOST_CHECK(fuzzy_equal_to()(180.0,newick::CalcDenominator(newick::StringToNewick("((1,2),(3,4))"),10.0)));
  BOOST_CHECK(fuzzy_equal_to()(180.0,newick::CalcDenominator(newick::StringToNewick("((4,1),(2,3))"),10.0)));
  BOOST_CHECK(fuzzy_equal_to()(180.0,newick::CalcDenominator(newick::StringToNewick("((3,4),(1,2))"),10.0)));
  BOOST_CHECK(fuzzy_equal_to()(180.0,newick::CalcDenominator(newick::StringToNewick("((2,3),(4,1))"),10.0)));
  //Test GetRootBranches
  {
    const std::vector<std::vector<int> > v = newick::GetRootBranches(newick::StringToNewick("(1,2)"));
    BOOST_CHECK(v.size() == 2);
    BOOST_CHECK(std::find(v.begin(),v.end(),
      newick::StringToNewick("(1)")) != v.end());
    BOOST_CHECK(std::find(v.begin(),v.end(),
      newick::StringToNewick("(2)")) != v.end());
  }
  {
    const std::vector<std::vector<int> > v = newick::GetRootBranches(newick::StringToNewick("(1,(2,3))"));
    BOOST_CHECK(v.size() == 2);
    BOOST_CHECK(std::find(v.begin(),v.end(),
      newick::StringToNewick("(1)")) != v.end());
    BOOST_CHECK(std::find(v.begin(),v.end(),
      newick::StringToNewick("(2,3)")) != v.end());
  }
  {
    const std::vector<std::vector<int> > v = newick::GetRootBranches(newick::StringToNewick("(1,2,(3,4))"));
    BOOST_CHECK(v.size() == 3);
    BOOST_CHECK(std::find(v.begin(),v.end(),
      newick::StringToNewick("(1)")) != v.end());
    BOOST_CHECK(std::find(v.begin(),v.end(),
      newick::StringToNewick("(2)")) != v.end());
    BOOST_CHECK(std::find(v.begin(),v.end(),
      newick::StringToNewick("(3,4)")) != v.end());
  }
  //Compare C++98 and C++0x version
  {
    const std::vector<std::string> v = newick::CreateValidBinaryNewicks();
    for(const std::string& s: v)
    {
      const std::vector<int> n = newick::StringToNewick(s);
      BOOST_CHECK(newick::GetRootBranches(n) == NewickCpp98().GetRootBranches(n));
    }
  }

  //Check if binary and trinary Newicks are detected correctly
  {
    const std::vector<std::string> v = newick::CreateValidBinaryNewicks();
    for(const std::string& s: v)
    {
      const std::vector<int> n = newick::StringToNewick(s);
      BOOST_CHECK(newick::IsBinaryNewick(n));
    }
  }
  //Check if unary Newicks are detected correctly
  {
    const std::vector<std::string> v = newick::CreateValidUnaryNewicks();
    for(const std::string& s: v)
    {
      const std::vector<int> n = newick::StringToNewick(s);
      BOOST_CHECK( newick::GetLeafMaxArity(n)<=1);
      BOOST_CHECK( newick::IsUnaryNewick(n));
      BOOST_CHECK(!newick::IsBinaryNewick(n));
      BOOST_CHECK(!newick::IsTrinaryNewick(n));
    }
  }
  //Check if binary Newicks are detected correctly
  {
    const std::vector<std::string> v = newick::CreateValidBinaryNewicks();
    for(const std::string& s: v)
    {
      const std::vector<int> n = newick::StringToNewick(s);
      BOOST_CHECK( newick::GetLeafMaxArity(n)<=2);
      BOOST_CHECK(!newick::IsUnaryNewick(n));
      BOOST_CHECK( newick::IsBinaryNewick(n));
      BOOST_CHECK(!newick::IsTrinaryNewick(n));
    }
  }
  //Check if trinary Newicks are detected correctly
  {
    const std::vector<std::string> v = newick::CreateValidTrinaryNewicks();
    for(const std::string& s: v)
    {
      //TRACE(s);
      const std::vector<int> n = newick::StringToNewick(s);
      BOOST_CHECK( newick::GetLeafMaxArity(n)<=3);
      BOOST_CHECK(!newick::IsUnaryNewick(n));
      BOOST_CHECK(!newick::IsBinaryNewick(n));
      BOOST_CHECK( newick::IsTrinaryNewick(n));
    }
  }
  //Test binary Newick
  {
    const std::string s("(1,(2,3))");
    const std::vector<std::vector<int> > n
      = newick::GetSimplerNewicks(newick::StringToNewick(s));
    //#define DEBUG_1_BO_1_2_3_BC
    #ifdef  DEBUG_1_BO_1_2_3_BC
    for(const auto& t: n)
    {
      TRACE(newick::NewickToString(t));
    }
    #endif
    BOOST_CHECK(n.size() == 2);
    BOOST_CHECK(std::find(n.begin(),n.end(),newick::StringToNewick("(1,(1,3))"))
      != n.end());
    BOOST_CHECK(std::find(n.begin(),n.end(),newick::StringToNewick("(1,(2,2))"))
      != n.end());
  }
  {
    const std::string s("(1,(2,3,4))");
    const std::vector<std::vector<int> > n
      = newick::GetSimplerNewicks(newick::StringToNewick(s));
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
    const std::vector<std::vector<int> > n
     = newick::GetSimplerNewicks(newick::StringToNewick(s));
    //#define DEBUG_1_BO_1_3_4_BC
    #ifdef  DEBUG_1_BO_1_3_4_BC
    TRACE(boost::lexical_cast<std::string>(n.size()));
    for(const auto& t: n)
    {
      TRACE(newick::NewickToString(t));
    }
    #endif
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
  {
    const std::string s("(1,(1,3,4))");
    const std::vector<std::pair<std::vector<int>,int> > n
      = NewickCpp98().GetSimplerNewicksFrequencyPairs(newick::StringToNewick(s));
    #ifdef TRACE_GETSIMPLERNEWICKSFREQUENCYPAIRS_1_134
    typedef std::pair<std::vector<int>,int> Pair;
    for(const Pair& p: n)
    {
      std::cout << newick::NewickToString(p.first) << '\n';
    }
    #endif
    BOOST_CHECK(n.size() == 4);
    BOOST_CHECK(std::find(n.begin(),n.end(),
      std::make_pair(newick::StringToNewick("(1,(4,4))"),1))
      != n.end());
    BOOST_CHECK(std::find(n.begin(),n.end(),
      std::make_pair(newick::StringToNewick("(1,(3,5))"),1))
      != n.end());
    BOOST_CHECK(std::find(n.begin(),n.end(),
      std::make_pair(newick::StringToNewick("(1,(1,2,4))"),3))
      != n.end());
    BOOST_CHECK(std::find(n.begin(),n.end(),
      std::make_pair(newick::StringToNewick("(1,(1,3,3))"),4))
      != n.end());
  }
  {
    const std::string s("((1,1),2)");
    const std::vector<std::vector<int> > n
      = newick::GetSimplerNewicks(
        newick::StringToNewick(s)
      );
    //#define DEBUG_BO_1_1_BC_2
    #ifdef  DEBUG_BO_1_1_BC_2
    for(const auto& t: n)
    {
      TRACE(newick::NewickToString(t));
    }
    #endif
    BOOST_CHECK(n.size() == 3);
    BOOST_CHECK(std::find(n.begin(),n.end(),
      newick::StringToNewick("(2,2)"))
      != n.end());
    BOOST_CHECK(std::find(n.begin(),n.end(),
      newick::StringToNewick("((1,1),1)"))
      != n.end());
  }
  {
    const std::string s("((1,1),2)");
    typedef std::pair<std::vector<int>,int> Pair;
    const std::vector<Pair> n
      = NewickCpp98().GetSimplerNewicksFrequencyPairs(newick::StringToNewick(s));
    #ifdef TRACE_GETSIMPLERNEWICKSFREQUENCYPAIRS_11_2
    for(const Pair& p: n)
    {
      std::clog << newick::NewickToString(p.first) << '\n';
    }
    #endif
    BOOST_CHECK(n.size() == 3);
    BOOST_CHECK(std::find(n.begin(),n.end(),
      std::make_pair(newick::StringToNewick("(2,2)"),1))
      != n.end());
    BOOST_CHECK(std::find(n.begin(),n.end(),
      std::make_pair(newick::StringToNewick("((1,1),1)"),2))
      != n.end());
  }
  {
    const std::string s("((2,1),4)");
    const std::vector<std::vector<int> > n
      = newick::GetSimplerNewicks(
        newick::StringToNewick(s)
      );
    BOOST_CHECK(n.size() == 3);
    BOOST_CHECK(std::find(n.begin(),n.end(),
      newick::StringToNewick("(3,4)"))
      != n.end());
    BOOST_CHECK(std::find(n.begin(),n.end(),
      newick::StringToNewick("((1,1),4)"))
      != n.end());
    BOOST_CHECK(std::find(n.begin(),n.end(),
      newick::StringToNewick("((2,1),3)"))
      != n.end());
  }
  {
    const std::string s("((2,1),4)");
    typedef std::pair<std::vector<int>,int> Pair;
    const std::vector<Pair> n
      = NewickCpp98().GetSimplerNewicksFrequencyPairs(newick::StringToNewick(s));
    #ifdef TRACE_GETSIMPLERNEWICKSFREQUENCYPAIRS_21_2
    for(const Pair& p: n)
    {
      TRACE(newick::NewickToString(p.first));
    }
    #endif
    BOOST_CHECK(n.size() == 3);
    BOOST_CHECK(std::find(n.begin(),n.end(),
      std::make_pair(newick::StringToNewick("(3,4)"),1))
      != n.end());
    BOOST_CHECK(std::find(n.begin(),n.end(),
      std::make_pair(newick::StringToNewick("((1,1),4)"),2))
      != n.end());
    BOOST_CHECK(std::find(n.begin(),n.end(),
      std::make_pair(newick::StringToNewick("((2,1),3)"),4))
      != n.end());
  }
  {
    const std::string s("((2,3),4)");
    const std::vector<std::vector<int> > n
      = newick::GetSimplerNewicks(
        newick::StringToNewick(s)
      );
    BOOST_CHECK(n.size() == 3);
    BOOST_CHECK(std::find(n.begin(),n.end(),
      newick::StringToNewick("((1,3),4)"))
      != n.end());
    BOOST_CHECK(std::find(n.begin(),n.end(),
      newick::StringToNewick("((2,2),4)"))
      != n.end());
    BOOST_CHECK(std::find(n.begin(),n.end(),
      newick::StringToNewick("((2,3),3)"))
      != n.end());
  }
  {
    const std::string s("((2,3),4)");
    typedef std::pair<std::vector<int>,int> Pair;
    const std::vector<Pair> n
      = NewickCpp98().GetSimplerNewicksFrequencyPairs(newick::StringToNewick(s));
    #ifdef TRACE_GETSIMPLERNEWICKSFREQUENCYPAIRS_23_4
    for(const Pair& p: n)
    {
      std::cout << newick::NewickToString(p.first) << '\n';
    }
    #endif
    BOOST_CHECK(n.size() == 3);
    BOOST_CHECK(std::find(n.begin(),n.end(),
      std::make_pair(newick::StringToNewick("((1,3),4)"),2))
      != n.end());
    BOOST_CHECK(std::find(n.begin(),n.end(),
      std::make_pair(newick::StringToNewick("((2,2),4)"),3))
      != n.end());
    BOOST_CHECK(std::find(n.begin(),n.end(),
      std::make_pair(newick::StringToNewick("((2,3),3)"),4))
      != n.end());
  }
  //Compare GetSimplerNewicks and
  //GetSimplerNewicksFrequencyPairs
  {
    const std::vector<std::string> newicks
      = newick::CreateValidNewicks();
    for(const std::string& newick_str: newicks)
    {
      const std::vector<int> newick
        = newick::StringToNewick(newick_str);
      const std::vector<std::vector<int> > v1
        = newick::GetSimplerNewicks(newick);
      const std::vector<std::pair<std::vector<int>,int> > v2
        = newick::GetSimplerNewicksFrequencyPairs(newick);
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
      BOOST_CHECK(s1 == s2);
      BOOST_CHECK(newick::GetSimplerNewicksFrequencyPairs(newick)
        == NewickCpp98().GetSimplerNewicksFrequencyPairs(newick));
    }
  }
}

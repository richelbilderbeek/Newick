


#include "newickcpp98.h"

#include <boost/numeric/conversion/cast.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

#include "fuzzy_equal_to.h"
#include "newick.h"



ribi::NewickCpp98::NewickCpp98()
{

}

///CreateValidTrinaryNewicks creates std::strings
///that can be converted to a TrinaryNewickVector.
///From http://www.richelbilderbeek.nl/CppCreateValidTinaryNewicks.htm
std::vector<std::string> ribi::NewickCpp98::CreateValidTrinaryNewicks() noexcept
{
  std::vector<std::string> v;
  v.push_back("(1,1,1)");
  v.push_back("(1,2,3)");
  v.push_back("((1,1),1,1)");
  v.push_back("(1,(1,1),1)");
  v.push_back("(1,1,(1,1))");
  v.push_back("(1,(2,3,4))");
  v.push_back("(1,2,(3,4))");
  v.push_back("(1,2,(3,4,5))");
  v.push_back("((1,2,3),4,5)");
  v.push_back("(11,22,33)");
  v.push_back("(11,(22,33,44))");
  v.push_back("(11,22,(33,44))");
  v.push_back("(11,22,(33,44,55))");
  v.push_back("((11,22,33),44,55)");
  v.push_back("((1,2),(3,4),(5,6))");
  v.push_back("((1,2,3),(4,5),(6,7))");
  v.push_back("((1,2),(3,4,5),(6,7))");
  v.push_back("((1,2),(3,4),(5,6,7))");
  v.push_back("((1,2,3),(4,5),(6,7))");
  v.push_back("((1,2),(3,4,5),(6,7))");
  v.push_back("((1,2),(3,4),(5,6,7))");
  v.push_back("((1,2,3),(4,5,6),(7,8))");
  v.push_back("((1,2),(3,4,5),(6,7,8))");
  v.push_back("((1,2,3),(4,5),(6,7,8))");
  v.push_back("((1,2,3),(4,5,6),(7,8,9))");
  v.push_back("((11,22,33),(44,55,66),(77,88,99))");
  return v;
}

//Need space between '> >', for C++98 support
std::vector<boost::tuple<std::string,double,double> > //!OCLINT Yes, too long, but I do not care about C++98
ribi::NewickCpp98::GetKnownProbabilities() noexcept
{
  std::vector<boost::tuple<std::string,double,double> > v;
  //Sum equals 1
  v.push_back(boost::make_tuple("(1)"  , 10.0, 1.0000000));
  //Sum equals 2
  v.push_back(boost::make_tuple("(2)"  , 10.0, 0.0909091));
  v.push_back(boost::make_tuple("(1,1)", 10.0, 0.9090909));
  //Sum equals 3
  v.push_back(boost::make_tuple("(3)"  , 10.0, 0.0151515));
  v.push_back(boost::make_tuple("(1,2)", 10.0, 0.0757576));
  v.push_back(boost::make_tuple("(2,1)", 10.0, 0.0757576));
  v.push_back(boost::make_tuple("(1,(1,1))", 10.0, 0.2525253));
  v.push_back(boost::make_tuple("((1,1),1)", 10.0, 0.2525253));
  //Trinary
  v.push_back(boost::make_tuple("(1,1,1)"  , 10.0, 0.7575758));
  //Sum equals 4
  v.push_back(boost::make_tuple("(4)"  , 10.0, 0.0034965));
  v.push_back(boost::make_tuple("(1,3)", 10.0, 0.0116550));
  v.push_back(boost::make_tuple("(2,2)", 10.0, 0.0058275));
  v.push_back(boost::make_tuple("(3,1)", 10.0, 0.0116550));
  v.push_back(boost::make_tuple("(1,(1,2))", 10.0, 0.0194250));
  v.push_back(boost::make_tuple("(1,(2,1))", 10.0, 0.0194250));
  v.push_back(boost::make_tuple("(2,(1,1))", 10.0, 0.0194250));
  v.push_back(boost::make_tuple("((1,2),1)", 10.0, 0.0194250));
  v.push_back(boost::make_tuple("((2,1),1)", 10.0, 0.0194250));
  v.push_back(boost::make_tuple("((1,1),2)", 10.0, 0.0194250));
  //Trinary
  v.push_back(boost::make_tuple("(1,1,2)", 10.0, 0.0582751));
  v.push_back(boost::make_tuple("(1,2,1)", 10.0, 0.0582751));
  v.push_back(boost::make_tuple("(2,1,1)", 10.0, 0.0582751));
  v.push_back(boost::make_tuple("(1,1,(1,1))", 10.0, 0.1295001));   //(1)(confirmed)
  v.push_back(boost::make_tuple("(1,(1,1),1)", 10.0, 0.1295001));
  v.push_back(boost::make_tuple("((1,1),1,1)", 10.0, 0.1295001));
  v.push_back(boost::make_tuple("(1,(1,1,1))", 10.0, 0.0971251));   //(2)(confirmed)
  v.push_back(boost::make_tuple("((1,1,1),1)", 10.0, 0.0971251));
  //Quadrary
  v.push_back(boost::make_tuple("(1,1,1,1)", 10.0, 0.5827505));
  //Sum equals 5
  v.push_back(boost::make_tuple("(1,4)", 10.0, 0.0024975));
  v.push_back(boost::make_tuple("(2,3)", 10.0, 0.0008325));
  v.push_back(boost::make_tuple("(3,2)", 10.0, 0.0008325));
  v.push_back(boost::make_tuple("(4,1)", 10.0, 0.0024975));
  v.push_back(boost::make_tuple("(1,(1,3))", 10.0, 0.0028305));
  v.push_back(boost::make_tuple("(1,(2,2))", 10.0, 0.0012950));
  v.push_back(boost::make_tuple("(1,(3,1))", 10.0, 0.0028305));
  v.push_back(boost::make_tuple("(2,(1,2))", 10.0, 0.0014338));
  v.push_back(boost::make_tuple("(2,(2,1))", 10.0, 0.0014338));
  v.push_back(boost::make_tuple("(3,(1,1))", 10.0, 0.0026640));
  //Trinary
  v.push_back(boost::make_tuple("(1,1,(1,2))"  , 10.0, 0.0092731));   //(3)(confirmed)
  v.push_back(boost::make_tuple("(1,1,(2,1))"  , 10.0, 0.0092731));
  v.push_back(boost::make_tuple("(1,1,(1,1,1))", 10.0, 0.0348263));   //(4)(confirmed)
  v.push_back(boost::make_tuple("(1,(1,1,1),1)", 10.0, 0.0348263));
  v.push_back(boost::make_tuple("((1,1,1),1,1)", 10.0, 0.0348263));
  v.push_back(boost::make_tuple("(2,(1,1,1))"  , 10.0, 0.0070069));   //(5)(confirmed)
  v.push_back(boost::make_tuple("((1,1,1),2)"  , 10.0, 0.0070069));
  v.push_back(boost::make_tuple("(1,1,1,(1,1))", 10.0, 0.0692918));   //(6)(confirmed)
  v.push_back(boost::make_tuple("(1,2,(1,1))"  , 10.0, 0.0092223));   //(7)(confirmed)
  v.push_back(boost::make_tuple("(2,1,(1,1))"  , 10.0, 0.0092223));
  v.push_back(boost::make_tuple("(1,(1,1),2)"  , 10.0, 0.0092223));
  v.push_back(boost::make_tuple("(2,(1,1),1)"  , 10.0, 0.0092223));
  v.push_back(boost::make_tuple("((1,1),1,2)"  , 10.0, 0.0092223));
  v.push_back(boost::make_tuple("((1,1),2,1)"  , 10.0, 0.0092223));
  v.push_back(boost::make_tuple("(1,(1,1,2))"  , 10.0, 0.0069190));   //(9)(confirmed)
  v.push_back(boost::make_tuple("(1,(1,2,1))"  , 10.0, 0.0069190));
  v.push_back(boost::make_tuple("(1,(2,1,1))"  , 10.0, 0.0069190));
  v.push_back(boost::make_tuple("((1,1,2),1)"  , 10.0, 0.0069190));
  v.push_back(boost::make_tuple("((1,2,1),1)"  , 10.0, 0.0069190));
  v.push_back(boost::make_tuple("((2,1,1),1)"  , 10.0, 0.0069190));
  //Quadrary
  v.push_back(boost::make_tuple("(1,(1,1,1,1))", 10.0, 0.0415140));   //(8)(confirmed)
  //Pentary
  v.push_back(boost::make_tuple("(1,1,1,1,1)"  , 10.0, 0.4162504));
  //Sum equals 6
  v.push_back(boost::make_tuple("(1,5)", 10.0, 0.0006660));
  v.push_back(boost::make_tuple("(2,4)", 10.0, 0.0001665));
  v.push_back(boost::make_tuple("(3,3)", 10.0, 0.0001110));
  v.push_back(boost::make_tuple("(1,(1,4))", 10.0, 0.0005804));
  v.push_back(boost::make_tuple("(1,(2,3))", 10.0, 0.0001679));
  v.push_back(boost::make_tuple("(1,(3,2))", 10.0, 0.0001679));
  v.push_back(boost::make_tuple("(1,(4,1))", 10.0, 0.0005804));
  v.push_back(boost::make_tuple("(2,(1,3))", 10.0, 0.0001991));
  v.push_back(boost::make_tuple("(2,(2,2))", 10.0, 0.0000925));
  v.push_back(boost::make_tuple("(2,(3,1))", 10.0, 0.0001991));
  v.push_back(boost::make_tuple("(3,(1,2))", 10.0, 0.0001880));
  v.push_back(boost::make_tuple("(3,(2,1))", 10.0, 0.0001880));
  v.push_back(boost::make_tuple("(4,(1,1))", 10.0, 0.0005043));
  //Trinary
  v.push_back(boost::make_tuple("(1,1,(1,3))"  , 10.0, 0.0012712));
  v.push_back(boost::make_tuple("(1,1,(2,2))"  , 10.0, 0.0005563));
  v.push_back(boost::make_tuple("(1,1,(3,1))"  , 10.0, 0.0012712));
  v.push_back(boost::make_tuple("(1,(1,3),1)"  , 10.0, 0.0012712));
  v.push_back(boost::make_tuple("(1,(2,2),1)"  , 10.0, 0.0005563));
  v.push_back(boost::make_tuple("(1,(3,1),1)"  , 10.0, 0.0012712));
  v.push_back(boost::make_tuple("((1,3),1,1)"  , 10.0, 0.0012712));
  v.push_back(boost::make_tuple("((2,2),1,1)"  , 10.0, 0.0005563));
  v.push_back(boost::make_tuple("((3,1),1,1)"  , 10.0, 0.0012712));
  v.push_back(boost::make_tuple("(1,2,(1,2))"  , 10.0, 0.0006346));
  v.push_back(boost::make_tuple("(2,1,(1,2))"  , 10.0, 0.0006346));
  v.push_back(boost::make_tuple("(1,2,(2,1))"  , 10.0, 0.0006346));
  v.push_back(boost::make_tuple("(2,1,(2,1))"  , 10.0, 0.0006346));
  v.push_back(boost::make_tuple("(1,(2,1),2)"  , 10.0, 0.0006346));
  v.push_back(boost::make_tuple("(1,(1,2),2)"  , 10.0, 0.0006346));
  v.push_back(boost::make_tuple("(2,(2,1),1)"  , 10.0, 0.0006346));
  v.push_back(boost::make_tuple("(2,(1,2),1)"  , 10.0, 0.0006346));
  v.push_back(boost::make_tuple("(1,3,(1,1))"  , 10.0, 0.0011913));
  v.push_back(boost::make_tuple("(1,(1,1),3)"  , 10.0, 0.0011913));
  v.push_back(boost::make_tuple("((1,1),1,3)"  , 10.0, 0.0011913));
  v.push_back(boost::make_tuple("(3,(1,1),1)"  , 10.0, 0.0011913));
  v.push_back(boost::make_tuple("((1,1),3,1)"  , 10.0, 0.0011913));
  v.push_back(boost::make_tuple("(1,1,(1,1,2))", 10.0, 0.0023165));
  v.push_back(boost::make_tuple("(1,1,(1,2,1))", 10.0, 0.0023165));
  v.push_back(boost::make_tuple("(1,1,(2,1,1))", 10.0, 0.0023165));
  v.push_back(boost::make_tuple("(1,(1,1,2),1)", 10.0, 0.0023165));
  v.push_back(boost::make_tuple("(1,(1,2,1),1)", 10.0, 0.0023165));
  v.push_back(boost::make_tuple("(1,(2,1,1),1)", 10.0, 0.0023165));
  v.push_back(boost::make_tuple("((1,1,2),1,1)", 10.0, 0.0023165));
  v.push_back(boost::make_tuple("((1,2,1),1,1)", 10.0, 0.0023165));
  v.push_back(boost::make_tuple("((2,1,1),1,1)", 10.0, 0.0023165));
  v.push_back(boost::make_tuple("(1,2,(1,1,1))", 10.0, 0.0023323));
  v.push_back(boost::make_tuple("(2,1,(1,1,1))", 10.0, 0.0023323));
  v.push_back(boost::make_tuple("(1,(1,1,1),2)", 10.0, 0.0023323));
  v.push_back(boost::make_tuple("(2,(1,1,1),1)", 10.0, 0.0023323));
  v.push_back(boost::make_tuple("((1,1,1),1,2)", 10.0, 0.0023323));
  v.push_back(boost::make_tuple("((1,1,1),2,1)", 10.0, 0.0023323));
  //Quadrary
  v.push_back(boost::make_tuple("(1,(1,1,1,2))"  , 10.0, 0.0027574));
  v.push_back(boost::make_tuple("(1,(1,1,2,1))"  , 10.0, 0.0027574));
  v.push_back(boost::make_tuple("(1,(1,2,1,1))"  , 10.0, 0.0027574));
  v.push_back(boost::make_tuple("(1,(2,1,1,1))"  , 10.0, 0.0027574));
  v.push_back(boost::make_tuple("((1,1,1,2),1)"  , 10.0, 0.0027574));
  v.push_back(boost::make_tuple("((1,1,2,1),1)"  , 10.0, 0.0027574));
  v.push_back(boost::make_tuple("((1,2,1,1),1)"  , 10.0, 0.0027574));
  v.push_back(boost::make_tuple("((2,1,1,1),1)"  , 10.0, 0.0027574));
  v.push_back(boost::make_tuple("(2,(1,1,1,1))"  , 10.0, 0.0028154));
  v.push_back(boost::make_tuple("((1,1,1,1),2)"  , 10.0, 0.0028154));
  //Pentary
  v.push_back(boost::make_tuple("(1,(1,1,1,1,1))", 10.0, 0.0183824));   //(7)
  v.push_back(boost::make_tuple("((1,1,1,1,1),1)", 10.0, 0.0183824));
  //Hexary
  v.push_back(boost::make_tuple("(1,1,1,1,1,1)"  , 10.0, 0.2775003));
  return v;
}

///GetRootBranches obtains the root branches from a non-unary Newick.
///Examples:
///(1,2)               -> { 1     , 2             }
///(1,2,3)             -> { 1     , 2     , 3     }
///((1,1),(2,2),(3,3)) -> { (1,1) , (2,2) , (3,3) }
///From http://www.richelbilderbeek.nl/CppGetRootBranchesBinary.htm
std::vector<std::vector<int> > //!OCLINT Yes, too long, but I do not care about C++98
  ribi::NewickCpp98::GetRootBranches(const std::vector<int>& n)
{
  assert(newick::IsNewick(n));
  assert(!newick::IsUnaryNewick(n));

  const int size = boost::numeric_cast<int>(n.size());
  std::vector<std::vector<int> > v;

  if (newick::IsSimple(n))
  {
    for (int i=1; i!=size-1; ++i) //Skip brackets
    {
      v.push_back(
        newick::CreateVector(
          static_cast<int>(newick::bracket_open),
          n[i],
          static_cast<int>(newick::bracket_close)
        )
      );
    }
    assert(newick::IsNewick(v.back()));
    assert(v.size() > 1);
    return v;
  }
  //Complex newick
  assert(!newick::IsSimple(n));
  const std::vector<int> depth = newick::GetDepth(n);

  assert(depth.size() == n.size());
  //Search for open and closing brackets in depth 1
  for (int i=0; i!=size; ++i)
  {
    if (depth[i] == 0 && n[i] > 0)
    {
      //C++0x initialization list

      std::vector<int> tmp = newick::Surround(n[i]);
      v.push_back(tmp);
      assert(newick::IsNewick(v.back()));
      continue;
    }
    if (depth[i] != 1 || n[i]!=newick::bracket_open) continue;
    for (int j=i+1; j!=size; ++j)
    {
      if (depth[j] != 1 || n[j]!=newick::bracket_close) continue;
      std::vector<int> w;
      w.push_back(newick::bracket_open);
      std::copy(n.begin() + i + 1,n.begin() + j,std::back_inserter(w));
      w.push_back(newick::bracket_close);
      assert(newick::IsNewick(w));
      v.push_back(w);
      //Set from index i after current end
      i = j;
      break; //!OCLINT Sure, this can be avoided, but I do not care about C++98
    }
  }
  assert(v.size() > 1);
  return v;
}

///GetSimplerNewicksFrequencyPairs creates simpler, derived Newicks from a Newick.
///Its simpler Newicks are identical to those created by GetSimplerNewicks.
///From http://www.richelbilderbeek.nl/CppGetSimplerNewicksFrequencyPairs.htm
std::vector<std::pair<std::vector<int>,int> > //!OCLINT Yes, too long, but I do not care about C++98
  ribi::NewickCpp98::GetSimplerNewicksFrequencyPairs(const std::vector<int>& n)
{
  assert(newick::IsNewick(n));

  std::vector<std::pair<std::vector<int>,int> > newicks;
  const std::vector<int> depths = newick::GetDepth(n);


  const int size = boost::numeric_cast<int>(n.size());
  for (int i = 0; i!=size; ++i)
  {
    assert(i >= 0);
    assert(i < size);
    if (n[i] < 1) continue;
    assert(n[i] > 0);
    if (n[i] > 1)
    {
      std::vector<int> new_newick(n);
      --new_newick[i];
      newicks.push_back( std::make_pair(new_newick,n[i]) );
      continue;
    }
    assert(n[i] == 1); //Most difficult...
    const int depth = depths[i];
    //j must first decrement, later increment with the same code
    int j_end  = -1;
    int j_step = -1;
    for (int j=i-1; ; j+=j_step)
    {
      //j must first decrement, later increment with the same code
      if (j == j_end
        //|| depths[j] < depth
        || (depths[j] == depth && n[j] < 0))
      {
        if (j_step == -1)
        {
          j = i + 1;
          j_end = size;
          j_step = 1;
        }
        else
        {
          break;
        }
      }
      assert(i!=j);
      assert(j >= 0);
      assert(j < size);
      //Only take frequencies of the same depth into account
      if (n[j] < 1 || depths[j] != depth) continue;
      std::vector<int> new_newick_with_zero(n);
      --new_newick_with_zero[i];
      assert(new_newick_with_zero[i] == 0);
      ++new_newick_with_zero[j];
      //If there is only one or two values between
      //the brackets, and one of these values was a
      //1 becoming added to the other, nullify the
      //1 and both brackets:
      //'((1,1),2)' -> '(00102)' -> '(1,2)'
      if (std::abs(i - j) == 1)
        //|| (std::abs(i - j) == 2 && n[i] == 1))
      {
        const int index_bracket_open  = std::min(i,j) - 1;
        const int index_bracket_close = std::max(i,j) + 1;
        if ( new_newick_with_zero[index_bracket_open]  == newick::bracket_open
          && new_newick_with_zero[index_bracket_close] == newick::bracket_close)
        {
          new_newick_with_zero[index_bracket_open]  = 0;
          new_newick_with_zero[index_bracket_close] = 0;
        }
      }
      //Remove decremented i and possibly nullified brackets
      std::vector<int> new_newick;
      std::remove_copy(
        new_newick_with_zero.begin(),
        new_newick_with_zero.end(),
        std::back_inserter(new_newick),
        0);
      //Add brackets if these are removed
      if (new_newick.front() != newick::bracket_open
        || new_newick.back() != newick::bracket_close)
      {
        new_newick = newick::Surround(new_newick);
      }
      assert(newick::IsNewick(new_newick));
      newicks.push_back(std::make_pair(new_newick, 1));
      continue;  //!OCLINT Sure, this can be avoided, but I do not care about C++98
    }
  }
  return newicks;
}

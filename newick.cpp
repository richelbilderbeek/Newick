//---------------------------------------------------------------------------
/*
Newick, Newick functions
Copyright (C) 2010-2015 Richel Bilderbeek

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
//From http://www.richelbilderbeek.nl/CppNewick.htm
//---------------------------------------------------------------------------
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#include "newick.h"

#include <algorithm>
#include <cassert>
#include <deque>
#include <iostream>
#include <functional>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <sstream>
#include <vector>


#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include "BigIntegerLibrary.hh"

#include "newickcpp98.h"
#pragma GCC diagnostic pop

//From http://www.richelbilderbeek.nl/CppAccumulate_if.htm
template
  <
  typename InputIterator,
  typename ElementType,
  typename Predicate
  >
const ElementType Accumulate_if(
  InputIterator first,
  const InputIterator last,
  ElementType init,
  const Predicate predicate)
{
  for (; first != last; ++first)
    if (predicate(*first)) init += *first;
  return init;
}

BigInteger ribi::newick::CalcComplexity(const std::vector<int>& v)
{
  if (v.empty()) return 0;
  //assert(IsNewick(v));
  BigInteger complexity = 1;
  int n_frequencies = 0;
  const int sz = v.size();
  for (int i=0; i!=sz; ++i)
  {
    const int x = v[i];
    if (x < 0) continue; //Ignore if x is not a number
    ++n_frequencies;
    complexity*=x;
  }
  complexity*=n_frequencies;
  return complexity;
}

double ribi::newick::CalcDenominator(const std::vector<int>& v,const double theta)
{
  int sum_above_zero = 0;
  int sum_above_one  = 0;
  for(const int& i: v)
  {
    if (i > 0) sum_above_zero+= i;
    if (i > 1) sum_above_one += i;
  }
  return boost::numeric_cast<double>(
      sum_above_zero * (sum_above_zero - 1))
    + (boost::numeric_cast<double>(sum_above_one) * theta)
  ;
}

BigInteger ribi::newick::CalcNumOfCombinationsBinary(const std::vector<int>& v)
{
  assert(IsNewick(v));

  //Get all positives
  std::vector<BigInteger> positives;
  std::copy_if(
    std::begin(v),std::end(v),
    std::back_inserter(positives),
    std::bind2nd(std::greater<BigInteger>(),0)
  );

  //Obtain numerator = (SUM(x))!
  const int sum_values = Accumulate_if(v.begin(),v.end(),0,std::bind2nd(std::greater<int>(),0));

  BigInteger numerator = FactorialBigInt(sum_values);

  //Obtain factorialated positives
  BigInteger denominator = 1;
  for(const int& i: v)
  {
    if (i<=0) continue;
    const BigInteger i_temp = FactorialBigInt(i);
    denominator*=i_temp;
  }

  //Obtain number_of_symmetries
  const BigInteger number_of_symmetries = CalcNumOfSymmetriesBinary(v);

  //Add number_of_symmetries times a 2 to denominator terms
  for(BigInteger i=0; i!=number_of_symmetries; ++i)
  {
    denominator*=2;
  }

  //Return the division
  numerator/=denominator;
  return numerator;
}

BigInteger ribi::newick::CalcNumOfSymmetriesBinary(std::vector<int> v)
{
  assert(IsNewick(v));
  assert(IsBinaryNewick(v));
  if (v.size() == 3) return 0;
  if (v.size() == 4) return (v[1] > 0 && v[1]==v[2] ? 1 : 0);

  const int n_reserved
    = *std::max_element(std::begin(v),std::end(v))
    + std::count_if(v.begin(), v.end(), std::bind2nd(std::greater<int>(),0));

  BigInteger n_symmetries = 0;
  int id = n_reserved + 1;

  std::map<std::pair<int,int>,int> ids;

  while (1)
  {
    //std::copy(v.begin(),v.end(),std::ostream_iterator<int>(std::clog," ")); std::clog << '\n';
    //Count number of symmetries
    assert(!v.empty());
    const std::size_t sz = v.size();
    assert(sz >= 2);
    const std::size_t j = sz - 1;
    for (std::size_t i = 0; i!=j; ++i)
    {
      if (v[i] > 0 && v[i]==v[i+1]) ++n_symmetries;
    }
    //Collect all leaves and store new leaves
    //std::vector<std::pair<int,int>> leaves;
    for (std::size_t i = 0; i!=j; ++i)
    {
      if (v[i] > 0 && v[i+1] > 0)
      {
        //Keep pair sorted
        const std::pair<int,int> p
          = (v[i] <= v[i+1]
          ? std::make_pair(v[i+0],v[i+1])
          : std::make_pair(v[i+1],v[i+0]));
        //If this leaf is new, store it
        if (ids.find(p)==ids.end())
        {
          ids[p] = id;
          ++id;
        }
      }
    }
    //Generalize all leaves
    for (std::size_t i = 0; i < v.size()-1; ++i)
    {
      assert(v.size()>2);
      if (v[i] > 0 && v[i+1] > 0)
      {
        //Keep pair sorted
        const std::pair<int,int> p
          = (v[i] <= v[i+1]
          ? std::make_pair(v[i+0],v[i+1])
          : std::make_pair(v[i+1],v[i+0]));
        //If this leaf is new, store it
        assert(ids.find(p)!=ids.end() && "Leaf should have been stored already");
        assert(i > 0);
        std::vector<int> v_new;
        std::copy(v.begin(),v.begin() + i - 1,std::back_inserter(v_new));
        const int id = ids[p];
        v_new.push_back(id);
        std::copy(v.begin() + i + 3,v.end(),std::back_inserter(v_new));
        v = v_new;
        i = (i-1 > 0 ? i-1 : 0);
        //Check if there are more leaves to be generalized
        if (v.size()<=4)
        {
          //Check if the last (X,Y) is symmetrical...
          return n_symmetries + (v[1] > 0 && v[1]==v[2] ? 1 : 0);
        }
      }
    }
  }
}

double ribi::newick::CalcProbabilitySimpleNewick(
  const std::vector<int>& v,
  const double theta)
{
  assert(newick::IsNewick(v));
  assert(IsSimple(v));
  const int sz = v.size();

  int n=0;
  int k=0;

  double probability = 1.0;

  for (int i=0; i!=sz; ++i)
  {
    if (v[i]>0)
    {
      const int ni = v[i];
      ++k;
      ++n;
      for (int p=1; p!=ni; ++p, ++n)
      {
        probability *= (static_cast<double>(p)
          / ( static_cast<double>(n) + theta));
      }
      probability /= ( static_cast<double>(n) + theta);
    }
  }
  probability *= (static_cast<double>(n)+theta)
    * std::pow(theta,static_cast<double>(k-1));
  return probability;
}

void ribi::newick::CheckNewickForMinimalSize(const std::string& s)
{
  if (s.size()<3)
  {
    throw std::invalid_argument(
      "The Newick std::string must have a size of "
      "at least three characters"
    );
  }
}

void ribi::newick::CheckNewickForOpeningBracket(const std::string& s)
{
  if (s[0]!='(')
  {
    throw std::invalid_argument(
      "The Newick std::string must start with "
      "an opening bracket ('(')."
    );
  }
}

void ribi::newick::CheckNewickForClosingBracket(const std::string& s)
{
  if (s[s.size()-1]!=')')
  {
    throw std::invalid_argument(
      "The Newick std::string must end with "
      "a closing bracket (')')."
    );
  }
}

void ribi::newick::CheckNewickForMatchingBrackets(const std::string& s)
{
  if (std::count(std::begin(s),std::end(s),'(')
    !=std::count(std::begin(s),std::end(s),')'))
  {
    throw std::invalid_argument(
       "The Newick std::string must have as much opening "
       "as closing brackets"
    );
  }
}

void ribi::newick::CheckNewickForZero(const std::string& s)
{
  if (s.find("(0")!=std::string::npos)
  {
    throw std::invalid_argument(
      "A std::string Newick frequency cannot be or "
      "start with a zero (#1)"
    );
  }
  if (s.find(",0")!=std::string::npos)
  {
    throw std::invalid_argument(
      "A std::string Newick frequency cannot be or "
      "start with a zero (#2)"
    );
  }
}

void ribi::newick::CheckNewickForBracketDistance(const std::string& s)
{
  if (s.find("()")!=std::string::npos)
  {
    throw std::invalid_argument(
      "The Newick std::string cannot have "
      "a consecutive opening and closing bracket"
    );
  }
}

void ribi::newick::CheckNewickForConsecutiveCommas(const std::string& s)
{
  if (s.find(",,") != std::string::npos)
  {
    throw std::invalid_argument(
      "A Newick std::string can have no consecutive comma's"
    );
  }
}

void ribi::newick::CheckNewickForCommaAfterBracketOpen(const std::string& s)
{
  if (s.find("(,")!=std::string::npos)
  {
    throw std::invalid_argument(
      "A Newick std::string cannot have comma "
      "directly after an opening bracket"
    );
  }
}


void ribi::newick::CheckNewickForCommaBeforeBracketClose(const std::string& s)
{
  if (s.find(",)")!=std::string::npos)
  {
    throw std::invalid_argument(
      "A Newick std::string cannot have comma "
      "directly before a closing bracket"
    );
  }
}

void ribi::newick::CheckNewick(const std::string& s)
{
  CheckNewickForMinimalSize(s);
  CheckNewickForOpeningBracket(s);
  CheckNewickForClosingBracket(s);
  CheckNewickForMatchingBrackets(s);
  CheckNewickForZero(s);
  CheckNewickForBracketDistance(s);
  CheckNewickForConsecutiveCommas(s);
  CheckNewickForCommaAfterBracketOpen(s);
  CheckNewickForCommaBeforeBracketClose(s);


  std::string s_copy = s;
  while(s_copy.size()>2) //Find a leaf and cut it until the string is empty
  {
    //Find a leaf
    //Find index i (starting opening bracket) and j (closing bracket)
    const std::size_t sz = s_copy.size();
    const auto p = FindOpeningAndClosingBracketIndices(s_copy);
    std::size_t i = p.first;
    std::size_t j = p.second;
    assert(s_copy[i]=='(');
    assert(s_copy[j]==')');
    //Check the range
    for (size_t k=i+1; k!=j; ++k)
    {
      if (!IsNumberOrComma(s_copy[k]))
      {
        std::stringstream err_msg;
        err_msg << "Invalid non-number character in input: '" << s_copy[k] << "'";
        throw std::invalid_argument(err_msg.str().c_str());
      }
    }
    if (i > 0 && s_copy[i-1]=='(' &&
      j +1 < sz && s_copy[j + 1] == ')')
    {
      throw std::invalid_argument(
        "Newicks must not have the form ((X))");
    }
    //Check if there is a comma somewhere between brackets
    if (i > 0 //< (1) is valid, (1,(2)) not, ((1),2) not
      && std::find(
        s_copy.begin()+i,s_copy.begin()+j,',')
          == s_copy.begin()+j)
    {
      throw std::invalid_argument(
        "The Newick std::string cannot have the sequence "
        "of an opening bracket, a value and a closing bracket "
        "as a \'complex\' leaf");
    }
    //Range is assumed valid
    //Cut the leaf (turns '(1,2)' to (9) )
    assert(s_copy[i]=='(');
    assert(s_copy[j]==')');
    const std::string s_new_1 = s_copy.substr(0,i);
    const std::string s_new_2 = s_copy.substr(j+1,sz-j-1); //OK
    const std::string s_new =  s_new_1 + "9" + s_new_2;
    s_copy = s_new;
  }
}

void ribi::newick::CheckNewickForMinimalSize(const std::vector<int>& v)
{
  if (v.size()<3)
  {
    throw std::invalid_argument(
      "The Newick std::vector<int> must have "
      "a size of at least three characters"
    );
  }
}

void ribi::newick::CheckNewickForOpeningBracket(const std::vector<int>& v)
{
  CheckNewickForMinimalSize(v);
  if (v[0]!=bracket_open)
  {
    throw std::invalid_argument(
      "The Newick std::vector<int> must start with "
      "an opening bracket ('(')."
    );
  }
}

void ribi::newick::CheckNewickForClosingBracket(const std::vector<int>& v)
{
  CheckNewickForMinimalSize(v);
  if (v.back() != bracket_close)
  {
    throw std::invalid_argument(
      "The Newick std::vector<int> must end with "
      "a closing bracket (')')."
    );
  }
}

void ribi::newick::CheckNewickForMatchingBrackets(const std::vector<int>& v)
{
  if (std::count(std::begin(v),std::end(v),static_cast<int>(bracket_open))
   != std::count(std::begin(v),std::end(v),static_cast<int>(bracket_close))
  )
  {
    throw std::invalid_argument(
      "The Newick std::string must have as much opening "
      "as closing brackets"
    );
  }
}

void ribi::newick::CheckNewickForZero(const std::vector<int>& v)
{
  if (std::count(std::begin(v),std::end(v),0))
  {
    throw std::invalid_argument(
      "A std::vector<int> Newick frequency cannot be zero"
    );
  }
}

void ribi::newick::CheckNewickForBracketDistance(const std::vector<int>& v)
{
  CheckNewickForMinimalSize(v);
  auto left = std::begin(v);
  auto right = std::begin(v);
  ++right;
  for ( ; right != std::end(v); ++left, ++right)
  {
    const auto a = *left;
    const auto b = *right;
    if (a == bracket_open && b == bracket_close)
    {
      throw std::invalid_argument(
        "The Newick std::vector<int> cannot have "
        "a consecutive opening and closing bracket"
      );
    }
  }
}

void ribi::newick::CheckNewick(const std::vector<int>& v)
{
  CheckNewickForMinimalSize(v);
  CheckNewickForOpeningBracket(v);
  CheckNewickForClosingBracket(v);
  CheckNewickForMatchingBrackets(v);
  CheckNewickForZero(v);
  CheckNewickForBracketDistance(v);

  std::vector<int> v_copy = v;
  while(v_copy.size()>2) //Find a leaf and cut it until the string is empty
  {
    //Find a leaf
    //Find index i (starting opening bracket) and j (closing bracket)
    const std::size_t sz = v_copy.size();
    std::size_t i = 0;
    std::size_t j = 0;
    for (i=0 ; i!=sz; ++i) //Index of opening bracket
    {
      if (v_copy[i]!=bracket_open) continue;

      assert(v_copy[i]==bracket_open);
      assert(i+1 < v_copy.size());
      assert(v_copy[i+1]!=bracket_close); //Already checked

      for (j=i+1; j!=sz; ++j)
      {
        if (v_copy[j]==bracket_open) { j = 0; break; }
        if (v_copy[j]!=bracket_close) continue;
        break;
      }
      if (i + 2 == j && j < sz - 1) //< (1) is valid, (1,(2)) not, ((1),2) not
      {
        throw std::invalid_argument(
          "The Newick std::vector<int> cannot have the sequence"
          "of an opening bracket, a value and a closing bracket"
          "as a \'complex\' leaf");
      }

      if (j ==  0) continue; //j cannot be 0 after previous for loop
      assert(j != sz); //CheckNewickForMatchingBrackets
      break;
    }
    assert(v_copy[i] == bracket_open); //CheckNewickForMatchingBrackets
    //Indices i and j found
    //Is range between i and j valid?
    assert(v_copy[i]==bracket_open);
    assert(v_copy[j]==bracket_close);
    //Check the range
    for (size_t k=i+1; k!=j; ++k)
    {
      if (v_copy[k] < 0)
      {
        std::ostringstream err_msg;
        err_msg << "Invalid non-number in input: '" << v_copy[k] << "'";
        throw std::invalid_argument(err_msg.str().c_str());
      }
    }
    //Range is assumed valid
    //Cut the leaf
    //Changes '(1,2)' to '(999)'
    assert(v_copy[i]==bracket_open);
    assert(v_copy[j]==bracket_close);
    std::vector<int> v_new(v_copy.begin(),v_copy.begin() + i);
    v_new.push_back(999);
    std::copy(v_copy.begin() + j + 1, v_copy.end(),std::back_inserter(v_new));
    v_copy = v_new;
  }
}

std::vector<std::string> ribi::newick::CreateInvalidNewicks() noexcept
{
  std::vector<std::string> v;
  v.push_back("");
  v.push_back("(");
  v.push_back(")");
  v.push_back("1");
  v.push_back("1234");
  v.push_back(")1234(");
  v.push_back("()1234()");
  v.push_back("(1234,)");
  v.push_back("(,1234,)");
  v.push_back("()");
  v.push_back("(0)");
  v.push_back("(-)");
  v.push_back("(-1)");
  v.push_back("(0,0)");
  v.push_back("(1,0)");
  v.push_back("(0,1)");
  v.push_back("(0,(1,1))");
  v.push_back("(1,(0,1))");
  v.push_back("(1,(1,0))");
  v.push_back("((0,1),1)");
  v.push_back("((1,0),1)");
  v.push_back("((1,1),0)");
  v.push_back("((2))");
  v.push_back("(1,(2,3)");
  v.push_back("(1,(2))");
  v.push_back("(1,((3)))");
  v.push_back("(11,(22,33)");
  v.push_back("(22,33),33)");
  v.push_back("1,2");
  v.push_back("(1,1),");
  v.push_back("(2,2),");
  v.push_back("((2,2),2),");
  v.push_back(",(1,1)");
  v.push_back(",(2,2)");
  v.push_back(",((2,2),2)");
  v.push_back(",(1,1),");
  v.push_back(",(2,2),");
  v.push_back("(2,(1,1)),");
  v.push_back(",((2,2),2),");
  v.push_back("(1,2");
  v.push_back("(-1,2");
  v.push_back("(1,-2");
  v.push_back("(0,-2");
  v.push_back("(-0,2");
  v.push_back("(1,,2)");
  v.push_back("(1,2))");
  v.push_back("(1,(2),3)");
  v.push_back("((1,2),(3,4))()");
  return v;
}

std::string ribi::newick::CreateRandomNewick(const int n,const int max) noexcept
{
  const std::vector<int> v = CreateRandomBinaryNewickVector(n,max);
  return NewickToString(v);
}

std::vector<int> ribi::newick::CreateRandomBinaryNewickVector(const int n,const int max) noexcept
{
  assert(n>0);
  assert(max>1);

  std::vector<int> v;
  v.reserve(2 + (n*2));

  v.push_back(bracket_open);

  v.push_back(1 + (std::rand() % (max-1) ));
  if (n==1)
  {
    v.push_back(bracket_close);
    return v;
  }
  v.push_back(1 + (std::rand() % (max-1) ));

  v.push_back(bracket_close); //??? IntVector format has no trailing bracket

  if (n==2)
  {
    return v;
  }

  for (int i=2; i!=n; ++i)
  {
    if ((std::rand() >> 4) % 2)
    {
      //Append
      std::vector<int> new_v;
      new_v.reserve(2 + v.size());
      new_v.push_back(bracket_open);

      std::copy(v.begin(),v.end(),std::back_inserter(new_v));

      new_v.push_back(1 + (std::rand() % (max-1)));
      new_v.push_back(bracket_close);
      std::swap(v,new_v);
    }
    else
    {
      //Prepend
      std::vector<int> new_v;
      new_v.reserve(2 + v.size());
      new_v.push_back(bracket_open);
      new_v.push_back(1 + (std::rand() % (max-1)));

      std::copy(v.begin(),v.end(),std::back_inserter(new_v));

      new_v.push_back(bracket_close);
      std::swap(v,new_v);
    }
    assert(std::count(v.begin(),v.end(),static_cast<int>(bracket_open ))
        == std::count(v.begin(),v.end(),static_cast<int>(bracket_close)));
  }
  assert(newick::IsNewick(v));
  return v;
}

std::vector<std::string> ribi::newick::CreateValidBinaryNewicks() noexcept
{
  return {
    "(1,2)",
    "(11,22)",
    "(1,(1,1))",
    "(1,(1,2))",
    "(1,(2,1))",
    "(1,(2,2))",
    "(1,(2,3))",
    "(2,(1,1))",
    "(2,(1,2))",
    "(2,(2,1))",
    "(2,(2,2))",
    "(4,(2,3))",
    "((2,3),4)",
    "(2,((2,3),4))",
    "(11,(22,33))",
    "((22,33),33)",
    "((1,2),(3,4))",
    "(((1,2),(3,4)),5)",
    "(1,((2,3),(4,5)))",
    "((11,2),(3,44))",
    "(((1,22),(33,4)),(55,6))"
  };
}

std::vector<std::string> ribi::newick::CreateValidTrinaryNewicks() noexcept
{
  return NewickCpp98().CreateValidTrinaryNewicks();
}

std::vector<std::string> ribi::newick::CreateValidNewicks() noexcept
{
  std::vector<std::string> v;
  {
    std::vector<std::string> w = CreateValidUnaryNewicks();
    std::copy(w.begin(),w.end(),std::back_inserter(v));
  }
  {
    std::vector<std::string> w = CreateValidBinaryNewicks();
    std::copy(w.begin(),w.end(),std::back_inserter(v));
  }
  {
    std::vector<std::string> w = CreateValidTrinaryNewicks();
    std::copy(w.begin(),w.end(),std::back_inserter(v));
  }
  {
    const std::vector<std::string> w = {
      "(1,2,3,4)",
      "(1,2,3,4,5)",
      "(1,2,3,4,5,6)",
      "(1,2,3,4,5,6,7)",
      "(1,2,3,4,5,6,7,8)",
      "((1,2,3,4,5,6,7,8),(1,2,3,4,5,6,7,8),(1,2,3,4,5,6,7,8))"
    };
    std::copy(w.begin(),w.end(),std::back_inserter(v));
  }
  return v;
}

std::vector<std::string> ribi::newick::CreateValidUnaryNewicks() noexcept
{
  return {
    "(1)",
    "(9)",
    "(123)"
  };
}

std::string ribi::newick::DumbNewickToString(const std::vector<int>& v) noexcept
{
  std::string s;
  s.reserve(2 * v.size()); //Just a guess
  const int sz = v.size();
  for (int i=0; i!=sz; ++i)
  {
    const int x = v[i];
    if (x >= 0)
    {
      s+=boost::lexical_cast<std::string>(x);
      const int next = v[i+1];
      if (next > 0 || next == bracket_open)
      {
        s+=",";
      }
    }
    else if (x==bracket_open)
    {
      s+="(";
    }
    else if (x==bracket_close)
    {
      s+=")";
      //Final closing bracket?
      if (i+1==sz) break;
      const int next = v[i+1];
      if (next > 0 || next == bracket_open)
      {
        s+=",";
      }
    }
    else
    {
      s+="x"; //Unknown character
    }
  }
  return s;
}

std::vector<int> ribi::newick::Factorial(const std::vector<int>& v_original) noexcept
{
  std::vector<int> v(v_original);
  std::transform(v.begin(),v.end(),v.begin(),
    [](const int i) { return Factorial(i); }
  );
  return v;
}

int ribi::newick::Factorial(const int n) noexcept
{
  assert(n>=0);
  int result = 1;
  for (int i=1; i<=n; ++i)
  {
    result*=i;
  }
  return result;
}

BigInteger ribi::newick::FactorialBigInt(const int n) noexcept
{
  assert(n>=0);
  BigInteger result = 1;
  for (int i=1; i<=n; ++i)
  {
    result*=i;
  }
  return result;
}

int ribi::newick::FindPosAfter(const std::vector<int>& v,const int x, const int index) noexcept
{
  assert(!v.empty());
  const int sz = v.size();
  for (int i=index+1; i!=sz; ++i)
  {
    assert(i >= 0);
    assert(i < sz);
    if (v[i]==x) return i;
  }
  return sz;
}

int ribi::newick::FindPosBefore(const std::vector<int>& v,const int x, const int index) noexcept
{
  for (int i=index-1; i!=-1; --i)
  {
    if (v[i]==x) return i;
  }
  return -1;
}

std::pair<std::size_t, std::size_t>
ribi::newick::FindOpeningAndClosingBracketIndices(
  const std::string& s
)
{
  CheckNewickForMatchingBrackets(s);
  const std::size_t sz = s.size();
  std::size_t i = 0;
  std::size_t j = 0;
  for (i=0; i!=sz; ++i) //Index of opening bracket
  {
    if (s[i]!='(') continue;

    assert(s[i]=='(');
    assert(i+1 < s.size());
    for (j=i+1; j!=sz; ++j)
    {
      if (s[j]=='(')
      {
        //j = 0;
        break;
      }
      if (s[j]!=')') continue;
      if (s[i]=='(' && s[j] == ')') return std::make_pair(i, j);
    }
    //if (j == 0) continue; //j cannot be 0 after previous for loop
    //assert(j != sz); // Must have as much opening as closing brackets #2"
    //break;
  }
  assert(!"Should not get here"); //!OCLINT valid idiom
  throw std::logic_error("could not find inner pair");
}

std::vector<int> ribi::newick::GetDepth(const std::vector<int>& n) noexcept
{
  assert(IsNewick(n));
  std::vector<int> v;
  int depth = -1;
  for(const int& x: n)
  {
    if (x == newick::bracket_open) ++depth;
    v.push_back(depth);
    if (x == newick::bracket_close) --depth;
  }
  assert(n.size() == v.size());
  return v;
}

std::vector<int> ribi::newick::GetFactorialTerms(const int n) noexcept
{
  assert(n > 0);
  std::vector<int> v(n);
  std::iota(std::begin(v), std::end(v), 1);
  assert(std::count(v.begin(),v.end(),0) == 0);
  return v;
}

std::vector<boost::tuple<std::string,double,double> > ribi::newick::GetKnownProbabilities() noexcept
{
  return NewickCpp98().GetKnownProbabilities();
}

int ribi::newick::GetLeafMaxArity(const std::vector<int>& n) noexcept
{
  assert(IsNewick(n));
  const int size = boost::numeric_cast<int>(n.size());
  if (IsSimple(n)) return size - 2;

  int max = 0;
  for (int from = 0; from!=size; ++from)
  {
    if (n[from] != newick::bracket_open) continue;
    for (int to = from+1; to!=size; ++to)
    {
      if (n[to] == newick::bracket_open) break;
      if (n[to]  > 0) continue;
      if (n[to] == newick::bracket_close)
      {
        assert(from < to);
        max = to - from - 1;
        break;
      }
    }
  }
  return max;
}

std::vector<std::vector<int> >
  ribi::newick::GetRootBranches(const std::vector<int>& n) noexcept
{
  return NewickCpp98().GetRootBranches(n);
}

std::pair<std::vector<int>,std::vector<int> >
  ribi::newick::GetRootBranchesBinary(const std::vector<int>& n) noexcept
{
  assert(IsNewick(n));
  assert(IsBinaryNewick(n) && "Only binary Newicks can have two roots");

  assert(n[0] == bracket_open);
  assert(n[n.size()-1] == bracket_close);
  const int sz = boost::numeric_cast<int>(n.size());
  //Return the answer directly is Newick consists
  //out of two values only
  if (sz==4)
  {
    assert(n[1] > 0);
    assert(n[2] > 0);
    return std::make_pair(
      newick::CreateVector(
        static_cast<int>(bracket_open),
        n[1],
        static_cast<int>(bracket_close)),
      newick::CreateVector(
        static_cast<int>(bracket_open),
        n[2],
        static_cast<int>(bracket_close)));
  }
  typedef std::vector<int>::const_iterator Iterator;
  const Iterator j = n.end() - 1;
  for (Iterator i = n.begin() + 2; i!=j; ++i)
  {
    std::vector<int> lhs(n.begin() + 1,i);
    std::vector<int> rhs(i,n.end() - 1);
    if ( lhs.front() != newick::bracket_open
      || lhs.back()  != newick::bracket_close)
    {
      lhs = newick::Surround(lhs);
    }
    if ( rhs.front() != newick::bracket_open
      || rhs.back()  != newick::bracket_close)
    {
      rhs = newick::Surround(rhs);
    }
    if (newick::IsNewick(lhs) && newick::IsNewick(rhs))
    {
      return std::make_pair(lhs,rhs);
    }
  }
  assert(!"Should not get here"); //!OCLINT
  throw std::logic_error("Should not get here");
}

std::vector<std::vector<int> >
  ribi::newick::GetSimplerBinaryNewicks(const std::vector<int>& n) noexcept
{
  assert(IsNewick(n));
  assert(IsUnaryNewick(n) || IsBinaryNewick(n));

  std::vector<std::vector<int> > newicks;

  //If newick is simple (by counting the number of opening brackets),
  //there are no simpler Newicks
  if (std::count( n.begin(),n.end(),
    static_cast<int>(bracket_open))==1)
  {
    //Simple newicks do not need to be simplified
    assert(n.size()==3 || n.size() == 4);
    assert(n[0]==bracket_open);
    assert(n[n.size()-1]==bracket_close);
    if (n.size() == 3)
    {
      if (n[1]>1)
      {
        std::vector<int> next(n);
        --next[1];
        assert(newick::IsNewick(next));
        newicks.push_back(next);
      }
      return newicks;
    }
    assert(n.size()==4);
    if (n[1] == 1)
    {
      std::vector<int> next
        = newick::CreateVector(
            static_cast<int>(bracket_open),
            n[2]+1,
            static_cast<int>(bracket_close)
          );
      assert(newick::IsNewick(next));
      newicks.push_back(next);
    }
    else
    {
      std::vector<int> next(n);
      --next[1];
      assert(newick::IsNewick(next));
      newicks.push_back(next);
    }
    if (n[2] == 1)
    {
      std::vector<int> next
        = newick::CreateVector(
            static_cast<int>(bracket_open),
            n[1]+1,
            static_cast<int>(bracket_close)
          );
      assert(newick::IsNewick(next));
      newicks.push_back(next);
    }
    else
    {
      std::vector<int> next(n);
      --next[2];
      assert(newick::IsNewick(next));
      newicks.push_back(next);
    }
    return newicks;
  }

  //newick is complex
  //Generate other Newicks and their coefficients
  const int sz = n.size();
  for (int i=0; i!=sz; ++i)
  {
    const int x = n[i];
    if (x < 0) //x is not a digit
    {
      continue;
    }
    if (x == 1)
    {
      //If x is not next to a digit, there is no simplification
      if (n[i-1]<0 && n[i+1]<1) continue;
      //If next to the x in a digit, merge these and remove their brackets
      //Is the 1 left of a value?
      if (n[i-1]==bracket_open)
      {
        assert(n[i+1] > 0);
        const int new_value = n[i+1] + 1;
        std::vector<int> next(n.begin(),n.begin() + i - 1);
        next.push_back(new_value);
        assert(n[i+2] < 0);
        std::copy(n.begin() + i + 3,n.end(),std::back_inserter(next));
        assert(newick::IsNewick(next));
        newicks.push_back(next);
      }
      else
      {
        //Is the 1 to the right of a value?
        assert(n[i-1] > 0); //< The other value
        const int new_value = n[i-1] + 1;
        std::vector<int> next(n.begin(),n.begin()+i-2);
        next.push_back(new_value);
        assert(n[i+1] < 0);
        std::copy(n.begin() + i + 2,n.end(),std::back_inserter(next));
        assert(newick::IsNewick(next));
        newicks.push_back(next);
      }
    }
    else
    {
      std::vector<int> next = n;
      --next[i];
      assert(newick::IsNewick(next));
      newicks.push_back(next);
    }
  }
  return newicks;
}

std::vector<std::vector<int>>
  ribi::newick::GetSimplerNewicksEasy(const std::vector<int>& n) noexcept
{
  assert(IsNewick(n));
  std::vector<std::vector<int>> newicks;
  const int size = boost::numeric_cast<int>(n.size());

  for (int i = 0; i!=size; ++i)
  {
    assert(i >= 0);
    assert(i < size);
    if (n[i] < 1) continue;
    assert(n[i] > 0);
    //If a frequency is above one, it is easy to create the simpler newicks
    if (n[i] > 1)
    {
      std::vector<int> new_newick(n);
      --new_newick[i];
      newicks.push_back(new_newick);
    }
  }
  return newicks;
}

std::vector<std::vector<int>>
  ribi::newick::GetSimplerNewicksHard(const std::vector<int>& n) noexcept
{
  assert(IsNewick(n));
  std::vector<std::vector<int>> newicks;
  const int size = boost::numeric_cast<int>(n.size());
  const std::vector<int> depths = GetDepth(n);

  //Go through all positions
  for (int i = 0; i!=size; ++i)
  {
    assert(i >= 0);
    assert(i < size);
    if (n[i] != 1) continue;
    //If a frequency is one, the Newick needs to be simplified
    assert(n[i] == 1); //Most difficult...
    const int depth = depths[i];
    //j must first decrement, later increment with the same code
    int j_end  = -1;
    int j_step = -1;
    for (int j=i-1; ; j+=j_step)
    {
      //j must first decrement, later increment with the same code
      if (j == j_end || (depths[j] == depth && n[j] < 0))
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
      //Remove brackets after possibly lonely value
      //If there is only one or two values between
      //the brackets, and one of these values was a
      //1 becoming added to the other, nullify the
      //1 and both brackets:
      //'((1,1),2)' -> '(00102)' -> '(1,2)'
      if (std::abs(i - j) == 1)
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
      assert(IsNewick(new_newick));
      newicks.push_back(new_newick);
      continue;
    }
  }
  return newicks;
}

std::vector<std::vector<int>>
  ribi::newick::GetSimplerNewicks(const std::vector<int>& n) noexcept
{
  assert(IsNewick(n));
  std::vector<std::vector<int>> newicks = GetSimplerNewicksEasy(n);
  const std::vector<std::vector<int>> hard_newicks = GetSimplerNewicksHard(n);
  std::copy(
    std::begin(hard_newicks),
    std::end(hard_newicks),
    std::back_inserter(newicks)
  );
  return newicks;
}

std::vector<std::pair<std::vector<int>,int> >
  ribi::newick::GetSimplerNewicksFrequencyPairs(
  const std::vector<int>& n
)
{
  return NewickCpp98().GetSimplerNewicksFrequencyPairs(n);
}


std::vector<std::pair<std::vector<int>,int> >
  ribi::newick::GetSimplerBinaryNewicksFrequencyPairs(
  const std::vector<int>& n) noexcept
{
  assert(IsNewick(n));
  assert(IsBinaryNewick(n));

  //If newick is simple (by counting the number of opening brackets),
  //there are no simpler Newicks
  if (IsSimple(n))
  {
    return GetSimplerBinaryNewicksFrequencyPairsSimple(n);
  }
  return GetSimplerBinaryNewicksFrequencyPairsComplex(n);
}

std::vector<std::pair<std::vector<int>,int> >
  ribi::newick::GetSimplerBinaryNewicksFrequencyPairsSimple(
  const std::vector<int>& n) noexcept
{
  //Simple newicks do not need to be simplified
  assert(IsSimple(n));
  assert(n.size()==3 || n.size() == 4);
  assert(n[0]==bracket_open);
  assert(n[n.size()-1]==bracket_close);
  std::vector<std::pair<std::vector<int>,int>> v;
  if (n.size() == 3)
  {
    if (n[1]>1)
    {
      std::vector<int> next(n);
      --next[1];
      assert(newick::IsNewick(next));
      v.push_back(std::make_pair(next,n[1]));
    }
    return v;
  }
  assert(n.size()==4);
  if (n[1] == 1)
  {
    std::vector<int> next
      = newick::CreateVector(
          static_cast<int>(bracket_open),
          n[2]+1,
          static_cast<int>(bracket_close)
        );
    assert(newick::IsNewick(next));
    v.push_back(std::make_pair(next,1));
  }
  else
  {
    std::vector<int> next(n);
    --next[1];
    assert(newick::IsNewick(next));
    v.push_back(std::make_pair(next,n[1]));
  }
  if (n[2] == 1)
  {
    std::vector<int> next
      = newick::CreateVector(
          static_cast<int>(bracket_open),
          n[1]+1,
          static_cast<int>(bracket_close)
        );
    assert(newick::IsNewick(next));
    v.push_back(std::make_pair(next,1));
  }
  else
  {
    std::vector<int> next(n);
    --next[2];
    assert(newick::IsNewick(next));
    v.push_back(std::make_pair(next,n[2]));
  }
  return v;
}

std::vector<std::pair<std::vector<int>,int> >
  ribi::newick::GetSimplerBinaryNewicksFrequencyPairsComplex(
  const std::vector<int>& n) noexcept
{
  assert(!IsSimple(n));
  std::vector<std::pair<std::vector<int>,int>> v;
  //newick is complex
  //Generate other Newicks and their coefficients
  const int sz = n.size();
  for (int i=0; i!=sz; ++i)
  {
    const int x = n[i];
    if (x < 0) //x is not a digit
    {
      continue;
    }
    if (x == 1)
    {
      //If x is not next to a digit, there is no simplification
      if (n[i-1]<0 && n[i+1]<1) continue;
      //If next to the x in a digit, merge these and remove their brackets
      //Is the 1 left of a value?
      if (n[i-1]==bracket_open)
      {
        assert(n[i+1] > 0);
        const int new_value = n[i+1] + 1;
        std::vector<int> next(n.begin(),n.begin() + i - 1);
        next.push_back(new_value);
        assert(n[i+2] < 0);
        std::copy(n.begin() + i + 3,n.end(),std::back_inserter(next));
        assert(newick::IsNewick(next));
        v.push_back(std::make_pair(next,x));
      }
      else
      {
        //Is the 1 to the right of a value?
        assert(n[i-1] > 0); //< The other value
        const int new_value = n[i-1] + 1;
        std::vector<int> next(n.begin(),n.begin()+i-2);
        next.push_back(new_value);
        assert(n[i+1] < 0);
        std::copy(n.begin() + i + 2,n.end(),std::back_inserter(next));
        assert(newick::IsNewick(next));
        v.push_back(std::make_pair(next,x));
      }
    }
    else
    {
      std::vector<int> next = n;
      --next[i];
      assert(newick::IsNewick(next));
      v.push_back(std::make_pair(next,x));
    }
  }
  return v;
}

std::string ribi::newick::GetNewickVersion() noexcept
{
  return "2.0";
}

std::vector<std::string> ribi::newick::GetNewickVersionHistory() noexcept
{
  return {
    "20xx-xx-xx: Version 1.0: initial version",
    "2013-05-29: Version 1.1: version history",
    "2016-10-19: Version 2.0: changed to free functions"
  };
}

void ribi::newick::InspectInvalidNewick(std::ostream& os, const std::vector<int>& v) noexcept
{
  os << "InspectInvalidNewick on: "
    << DumbNewickToString(v) << '\n';
  try
  {
    newick::CheckNewick(v);
  }
  catch (std::exception& e)
  {
    os << "Invalidity caused by: " << e.what() << '\n';
  }
}

bool ribi::newick::IsNewick(const std::string& s) noexcept
{
  try
  {
    CheckNewick(s);
  }
  catch (const std::exception&)
  {
    return false;
  }
  return true;
}

bool ribi::newick::IsNumberOrComma(const char c) noexcept
{
  const std::string s{"0123456789,"};
  return s.find(c) != std::string::npos;
}

bool ribi::newick::IsSimple(const std::vector<int>& v) noexcept
{
  assert(newick::IsNewick(v));
  //A Newick is simple if it contains no '(' after the initial one
  return std::count(
    v.begin()+1,v.end(),
    static_cast<int>(bracket_open)
  ) == 0;
}

bool ribi::newick::IsBinaryNewick(std::vector<int> v) noexcept
{
  assert(newick::IsNewick(v));
  if (IsUnaryNewick(v)) return false;

  while (1)
  {
    const int sz = boost::numeric_cast<int>(v.size());
    //Only numbers?
    if (IsSimple(v))
    {
      //Binary Newick has size 4, for example '(1,2)'
      return sz == 4;
    }
    if (GetLeafMaxArity(v) > 2) return false;
    v = newick::ReplaceLeave(v,42);
  }
}

bool ribi::newick::IsNewick(const std::vector<int>& v) noexcept
{
  try
  {
    CheckNewick(v);
  }
  catch (...)
  {
    return false;
  }
  return true;
}

///IsTrinaryNewick checks if a Newick is a trinary tree,
///that is: each node splits in three or less branches
///From http://www.richelbilderbeek.nl/CppIsTrinaryNewick.htm
bool ribi::newick::IsTrinaryNewick(std::vector<int> v) noexcept
{
  assert(newick::IsNewick(v));
  if (IsUnaryNewick(v)) return false;
  if (IsBinaryNewick(v)) return false;

  bool trinarity_found = false;

  while (1)
  {
    const int sz = boost::numeric_cast<int>(v.size());
    //Only numbers?
    if (IsSimple(v))
    {
      //Ternary Newick has size 5, for example '(1,2,3)'
      return trinarity_found || sz == 5;
    }
    const int leaf_max_arity = GetLeafMaxArity(v);
    if (leaf_max_arity > 3) return false;
    if (leaf_max_arity == 3) trinarity_found = true;

    v = newick::ReplaceLeave(v,42);
  }
}

bool ribi::newick::IsUnaryNewick(const std::vector<int>& v) noexcept
{
  assert(newick::IsNewick(v));
  return v.size() == 3
    && v[0] == newick::bracket_open
    && v[1] >  0
    && v[2] == newick::bracket_close;
}

std::string ribi::newick::NewickToString(const std::vector<int>& v)
{
  assert(v.size() > 2 && "A Newick must at least have one single value");
  assert(v[0] == bracket_open
    && "A std::vector<int> Newick must start with a bracket_open");
  assert(v[v.size() - 1] == bracket_close
    && "A std::vector<int> Newick must end with a bracket_close");
  std::string s;
  s.reserve(2 * v.size()); //Just a guess
  const int sz = v.size();
  for (int i=0; i!=sz; ++i)
  {
    const int x = v[i];
    if (x >= 0)
    {
      s+=boost::lexical_cast<std::string>(x);
      assert(i+1<sz && "Must not end with number");
      const int next = v[i+1];
      if (next > 0 || next == bracket_open)
      {
        s+=",";
      }
    }
    else if (x==bracket_open)
    {
      s+="(";
    }
    else if (x==bracket_close)
    {
      s+=")";
      //Final closing bracket?
      if (i+1==sz) break;
      const int next = v[i+1];
      if (next > 0 || next == bracket_open)
      {
        s+=",";
      }
    }
    else
    {
      assert(!"A std::vector<int> Newick must consist of brackets and values only"); //!OCLINT
    }
  }
  return s;
}

std::vector<int> ribi::newick::ReplaceLeave(
  const std::vector<int>& newick,
  const int value
)
{
  assert(IsNewick(newick) && "Only a valid Newick can have its leaves replaced");
  assert(!IsSimple(newick) && "There must a leaf to simplify");
  typedef std::vector<int>::const_iterator Iterator;
  const Iterator end = newick.end();
  for (Iterator from = newick.begin(); from!=end; ++from)
  {
    if (*from != newick::bracket_open) continue;

    for (Iterator to = from + 1; to!=end; ++to)
    {
      if (*to > 0) continue;
      if (*to == newick::bracket_open) break;
      if (*to == newick::bracket_close)
      {
        //Found
        std::vector<int> new_newick(newick.begin(),from);
        new_newick.push_back(value);
        std::copy(to + 1,newick.end(),std::back_inserter(new_newick));
        assert(newick::IsNewick(new_newick));
        return new_newick;
      }
    }
  }
  assert(!"Should not get here"); //!OCLINT
  throw std::logic_error("Should not get here");
}

std::vector<int> ribi::newick::StringToNewick(const std::string& newick)
{
  assert(IsNewick(newick));
  assert(!newick.empty()
    && "s must not be empty");
  assert(newick[              0]=='('
    && "s must begin with a '('");
  assert(newick[newick.size()-1]==')'
    && "s must end with a ')'");

  std::vector<int> v;
  int value = 0;

  for(const char i: newick)
  {
    if (i == '(')
    {
      if (value!=0) v.push_back(value);
      value = 0;
      v.push_back(bracket_open);
      continue;
    }
    if (i == ')')
    {
      if (value!=0) v.push_back(value);
      value = 0;
      v.push_back(bracket_close);
      continue;
    }
    if (i == ',')
    {
      if (value!=0) v.push_back(value);
      value = 0;
      continue;
    }
    assert(i >= '0' && i <= '9'); //Should be a number
    value*=10;
    value+=boost::lexical_cast<int>(i);
  }
  assert(value == 0 && "Final bracket close must set value to zero");
  return v;
}

std::vector<int> ribi::newick::Surround(const std::vector<int>& newick) noexcept
{
  std::vector<int> new_newick;
  new_newick.push_back(newick::bracket_open);
  std::copy(newick.begin(),newick.end(),std::back_inserter(new_newick));
  new_newick.push_back(newick::bracket_close);
  return new_newick;
}

std::vector<int> ribi::newick::Surround(const int f) noexcept
{
  return CreateVector(
    static_cast<int>(newick::bracket_open),
    f,
    static_cast<int>(newick::bracket_close)
  );
}

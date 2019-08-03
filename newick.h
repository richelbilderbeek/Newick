#ifndef NEWICK_H
#define NEWICK_H

#include <cmath>
#include <string>
#include <vector>

#pragma GCC diagnostic push



#include <boost/lexical_cast.hpp>
#include <boost/tuple/tuple.hpp>

#include "BigIntegerLibrary.hh"
#include "newickcpp98.h"
#include "newickstorage.h"
#pragma GCC diagnostic pop

namespace ribi {

namespace newick {

///CalcComplexity calculates the complexity of a Newick.
///From http://www.richelbilderbeek.nl/CppCalcComplexity.htm
BigInteger CalcComplexity(const std::vector<int>& v);

///CalcNumOfCombinations returns the number of combinations a Newick can have.
///
///The number of possible combinations equals
///     !(v0 + v1 + v2 + etc)
/// N = -------------------------- / 2^number_of_symmetries
///     !v0 * !v1 * !v2 * etc
///
///      n
///   = --- / 2^number_of_symmetries
///      d
///
/// where v denotes an element in vector v in range [1,-> >
/// where v0 denotes the first element in vector v
/// and where !v0 denotes the factorial of v0
///     {factorial(!SUM(v)) product terms}
/// N = --------------------------------------------------
///     {product terms} + { number_symmetries times a '2'}
///
///     numerator_terms
/// N = --------------------------------------------------
///     denominator_terms with appended number_symmetries times a '2'
///
///From http://www.richelbilderbeek.nl/CppCalcNumOfCombinationsBinary.htm
BigInteger CalcNumOfCombinationsBinary(const std::vector<int>& v);

///CalcNumOfSymmetries calculates the number of symmetries in a Newick.
///From http://www.richelbilderbeek.nl/CppCalcNumOfSymmetriesBinary.htm
BigInteger CalcNumOfSymmetriesBinary(std::vector<int> v);

double CalcDenominator(const std::vector<int>& v,const double theta);

///CalcProbabilitySimpleNewick returns the probability of
///a Newick for a value of theta
///using the Ewens formula
double CalcProbabilitySimpleNewick(const std::vector<int>& v,const double theta);

///Count the number of adjacent non-zero positive values
int CountAdjacentNonZeroPositives(const std::vector<int>& v);

///CreateInvalidNewicks creates std::strings
///that cannot and must not be converted to a Newick
///From http://www.richelbilderbeek.nl/CppCreateInvalidNewicks.htm
std::vector<std::string> CreateInvalidNewicks() noexcept;

///CreateRandomNewick creates an unsorted Newick string,
///with n values, with each value e [0,max>.
///From http://www.richelbilderbeek.nl/CppCreateRandomNewick.htm
std::string CreateRandomNewick(const int n,const int max) noexcept;

///Put the integers in a pair ordered
std::pair<int,int> CreateSortedPair(const int a, const int b) noexcept;

///CreateRandomBinaryNewickVector creates an unsorted Newick
///std::vector<int>, with n values, with each value e [0,max>.
///From http://www.richelbilderbeek.nl/CppCreateRandomBinaryNewickVector.htm
std::vector<int> CreateRandomBinaryNewickVector(const int n,const int max) noexcept;

///DumbNewickToString converts a Newick std::vector<int> to a
///standard-format std::string without error checking.
///From http://www.richelbilderbeek.nl/CppDumbNewickToString.htm
std::string DumbNewickToString(const std::vector<int>& v) noexcept;

///Factorial calculates the factorial of all std::vector elements.
///From http://www.richelbilderbeek.nl/CppFactorial.htm
std::vector<int> Factorial(const std::vector<int>& v_original) noexcept;

///FactorialBigInt returns the factorial of an integer
///as a BigInteger.
///From http://www.richelbilderbeek.nl/CppFactorialBigInt.htm
BigInteger FactorialBigInt(const int n) noexcept;

///Factorial calculates the factorial of a value.
///From http://www.richelbilderbeek.nl/CppFactorial.htm
int Factorial(const int n) noexcept;

int FindPosAfter(const std::vector<int>& v,const int index,const int value) noexcept;
int FindPosBefore(const std::vector<int>& v,const int index,const int value) noexcept;

///GetDepth returns the depth of each Newick element
///Example #1
///(1,2,3)
///01 1 10
///Example #2
///((1,2),(3,(4,5)))
///000 00 00 00 0000 <- depth layer 0 (comma's are skipped)
///.11 11 11 11 111. <- depth layer 1
///... .. .. .2 22.. <- depth layer 2
///011 11 11 22 2210 <- result of GetDepth
std::vector<int> GetDepth(const std::vector<int>& n) noexcept;


///GetFactorialTerms returns all terms from a factorial.
///For example, 4! return {1,2,3,4}
///From http://www.richelbilderbeek.nl/CppGetFactorialTerms.htm
std::vector<int> GetFactorialTerms(const int n) noexcept;

std::vector<boost::tuple<std::string,double,double> > GetKnownProbabilities() noexcept;
int GetLeafMaxArity(const std::vector<int>& n) noexcept;


///GetRootBranches obtains the root branches from a non-unary Newick.
///Examples:
///(1,2)               -> { 1     , 2             }
///(1,2,3)             -> { 1     , 2     , 3     }
///((1,1),(2,2),(3,3)) -> { (1,1) , (2,2) , (3,3) }
///From http://www.richelbilderbeek.nl/CppGetRootBranchesBinary.htm
std::vector<std::vector<int> >
  GetRootBranches(const std::vector<int>& n) noexcept;

///GetRootBranchesBinary obtains the two root branches from a binary Newick.
///Examples:
///(1,2)                 -> { 1             , 2     }
///(1,(2,3))             -> { 1             , (2,3) }
///((1,2),(3,4))         -> { (1,2)         , (3,4) }
///(((1,2),(3,4)),(5,6)) -> { ((1,2),(3,4)) , (5,6) }
///From http://www.richelbilderbeek.nl/CppGetRootBranchesBinary.htm
std::pair<std::vector<int>,std::vector<int> >
  GetRootBranchesBinary(const std::vector<int>& n) noexcept;

///GetSimplerBinaryNewicks creates simpler, derived Newicks from a binary Newick.
///From http://www.richelbilderbeek.nl/CppGetSimplerBinaryNewicks.htm
std::vector<std::vector<int>> GetSimplerBinaryNewicks(
  const std::vector<int>& n
) noexcept;

///Used by GetSimplerBinaryNewicks
std::vector<std::vector<int>> GetSimplerBinaryNewicksSimple(
  const std::vector<int>& n
) noexcept;

///Used by GetSimplerBinaryNewicks
std::vector<std::vector<int>> GetSimplerBinaryNewicksComplex(
  const std::vector<int>& n
) noexcept;

///GetSimplerBinaryNewicksFrequencyPairs creates simpler, derived Newicks from a
///binary Newick as well as the frequency that is simplified.
std::vector<std::pair<std::vector<int>,int> >
  GetSimplerBinaryNewicksFrequencyPairs(
  const std::vector<int>& n
) noexcept;

///GetSimplerBinaryNewicksFrequencyPairs for a simple newick.
std::vector<std::pair<std::vector<int>,int> >
  GetSimplerBinaryNewicksFrequencyPairsSimple(
  const std::vector<int>& n
) noexcept;

///GetSimplerBinaryNewicksFrequencyPairs for a complex newick.
std::vector<std::pair<std::vector<int>,int> >
  GetSimplerBinaryNewicksFrequencyPairsComplex(
  const std::vector<int>& n
) noexcept;

std::string GetNewickVersion() noexcept;
std::vector<std::string> GetNewickVersionHistory() noexcept;

///InspectInvalidNewick writes the cause of the Newick invalidity
///to the std::ostream.
///From http://www.richelbilderbeek.nl/CppInspectInvalidNewick.htm
void InspectInvalidNewick(std::ostream& os, const std::vector<int>& v) noexcept;

///IsBinaryNewick checks if a Newick is a binary tree,
///that is: each node splits in two (not more) branches
///From http://www.richelbilderbeek.nl/CppIsBinaryNewick.htm
bool IsBinaryNewick(std::vector<int> v) noexcept;

bool IsTrinaryNewick(std::vector<int> v) noexcept;

///IsUnaryNewick checks if a Newick is a unary tree,
///that is: there is only one node.
///From http://www.richelbilderbeek.nl/CppIsUnaryNewick.htm
bool IsUnaryNewick(const std::vector<int>& v) noexcept;

///IsSimple returns true if the Newick std::vector contains
///leaves only. For example, the Newick '(1,2,3)' is simple,
///the Newick '((1,2),3)' is not simple
///From http://www.richelbilderbeek.nl/CppIsNewick.htm
bool IsSimple(const std::vector<int>& v) noexcept;

///NewickToString converts a Newick std::vector<int> to a
///standard-format std::string.
///From http://www.richelbilderbeek.nl/CppNewickToString.htm
std::string NewickToString(const std::vector<int>& v);

///ReplaceLeave replaces the first leaf that it finds by a value.
///For example, using ReplaceLeave on '((1,2),(3,4))' with a value
///of 42 results in '(42,(3,4))'.
std::vector<int> ReplaceLeave(const std::vector<int>& newick, const int value);

///StringToNewick converts a std::string to a Newick std::vector<int>
///StringToNewick assumes that the input is well-formed and
///has both trailing and tailing brackets.
///From http://www.richelbilderbeek.nl/CppNewickToVector.htm
std::vector<int> StringToNewick(const std::string& newick);

enum { bracket_open  = -1 };
enum { bracket_close = -2 };
enum { comma         = -3 };
enum { new_line      = -4 };
enum { null          = -5 };

///CheckNewick checks if a std::vector<int> is a valid Newick.
///If this std::vector<int> is not a valid Newick,
///CheckNewick throws an exception with a detailed description
///From http://www.richelbilderbeek.nl/CppCheckNewick.htm
void CheckNewick(const std::vector<int>& v);

///CheckNewick checks if a std::string is a valid Newick.
///If this std::string is not a valid Newick,
///CheckNewick throws an exception with a detailed description
///From http://www.richelbilderbeek.nl/CppCheckNewick.htm
void CheckNewick(const std::string& s);

///Throws if Newick is too short to be valid
void CheckNewickForMinimalSize(const std::vector<int>& v);

///Throws if Newick is too short to be valid
void CheckNewickForMinimalSize(const std::string& s);

///Throws if Newick has no opening bracket
void CheckNewickForOpeningBracket(const std::vector<int>& v);

///Throws if Newick has no opening bracket
void CheckNewickForOpeningBracket(const std::string& s);

///Throws if Newick has no closing bracket
void CheckNewickForClosingBracket(const std::vector<int>& v);

///Throws if Newick has no closing bracket
void CheckNewickForClosingBracket(const std::string& s);

///Throws if Newick has an equal number of opening and closing brackets
void CheckNewickForMatchingBrackets(const std::vector<int>& v);

///Throws if Newick has an equal number of opening and closing brackets
void CheckNewickForMatchingBrackets(const std::string& s);

///Throws if Newick has a frequency of zero
void CheckNewickForZero(const std::vector<int>& v);

///Throws if Newick has a frequency of zero
void CheckNewickForZero(const std::string& s);

///Throws if Newick has no value between brackets, e.g '(1,())'
void CheckNewickForBracketDistance(const std::vector<int>& v);

///Throws if Newick has no value between brackets, e.g '(1,())'
void CheckNewickForBracketDistance(const std::string& s);

///Throws if Newick has two consecutive comma's, e.g '(1,,1)'
void CheckNewickForConsecutiveCommas(const std::string& s);

///Throws if Newick has a comma after a bracket open, e.g '(,1)'
void CheckNewickForCommaAfterBracketOpen(const std::string& s);

///Throws if Newick has a comma before a bracket close, e.g '(1,)'
void CheckNewickForCommaBeforeBracketClose(const std::string& s);

///CreateValidBinaryNewicks creates std::strings
///that can be converted to a BinaryNewickVector.
///From http://www.richelbilderbeek.nl/CppCreateValidBinaryNewicks.htm
std::vector<std::string> CreateValidBinaryNewicks() noexcept;

///CreateValidNewicks creates std::strings
///that are valid newicks.
///From http://www.richelbilderbeek.nl/CppCreateValidNewicks.htm
std::vector<std::string> CreateValidNewicks() noexcept;

///CreateValidTrinaryNewicks creates std::strings
///that can be converted to a TrinaryNewickVector.
///From http://www.richelbilderbeek.nl/CppCreateValidTinaryNewicks.htm
std::vector<std::string> CreateValidTrinaryNewicks() noexcept;

///CreateValidUnaryNewicks creates unary Newick std::strings
std::vector<std::string> CreateValidUnaryNewicks() noexcept;


///CreateVector creates a std::vector from three arguments
///From http://www.richelbilderbeek.nl/CppCreateVector.htm
template <class T>
std::vector<T> CreateVector(const T& a, const T& b, const T& c)
{
  std::vector<T> v;
  v.reserve(3);
  v.push_back(a);
  v.push_back(b);
  v.push_back(c);
  return v;
}

///Finds the indices of the innermost brackets
///Examples:
///(1,(1,1))
///   ^   ^
std::pair<std::size_t, std::size_t> FindOpeningAndClosingBracketIndices(
  const std::string& s
);

///Finds the indices of the innermost brackets
///Examples:
///(1,(1,1))
///   ^   ^
std::pair<std::size_t, std::size_t> FindOpeningAndClosingBracketIndices(
  const std::vector<int>& v
);


///GetSimplerNewicks creates simpler, derived Newicks from a Newick.
///From http://www.richelbilderbeek.nl/CppGetSimplerNewicks.htm
std::vector<std::vector<int>> GetSimplerNewicks(
  const std::vector<int>& n
) noexcept;

///Used by GetSimplerNewicks
std::vector<std::vector<int>> GetSimplerNewicksEasy(
  const std::vector<int>& n
) noexcept;

///GetSimplerNewicksFrequencyPairs creates simpler, derived Newicks from a Newick.
///Its simpler Newicks are identical to those created by GetSimplerNewicks.
///From http://www.richelbilderbeek.nl/CppGetSimplerNewicksFrequencyPairs.htm
std::vector<std::pair<std::vector<int>,int> >
  GetSimplerNewicksFrequencyPairs(
  const std::vector<int>& n
);

///Used by GetSimplerNewicks
std::vector<std::vector<int>> GetSimplerNewicksHard(
  const std::vector<int>& n
) noexcept;

///Used by GetSimplerNewicksHard
std::vector<std::vector<int>> GetSimplerNewicksHardFromHere(
  const std::vector<int>& n,
  const int i
) noexcept;


///IsNewick returns true if a std::string is a valid Newick
///and false otherwise.
///From http://www.richelbilderbeek.nl/CppIsNewick.htm
bool IsNewick(const std::string& s) noexcept;

///IsNewick returns true if a std::vector<int> is a valid Newick
///and false otherwise.
///From http://www.richelbilderbeek.nl/CppIsNewick.htm
bool IsNewick(const std::vector<int>& v) noexcept;

///Returns true if c is a number or comma
bool IsNumberOrComma(const char c) noexcept;

///Used by CalcNumOfSymmetriesBinary
void StoreAllNewLeafs(
  const std::vector<int>& v,
  std::map<std::pair<int,int>,int>& ids,
  int& id
);

///Surround surrounds the Newick with brackets
std::vector<int> Surround(const std::vector<int>& newick) noexcept;

///Surround surrounds the frequency with brackets
std::vector<int> Surround(const int f) noexcept;

template <class NewickType>
double CalculateProbability(
  const NewickType& n,
  const double theta,
  NewickStorage<NewickType>& storage
)
{
  while(1)
  {
    try
    {
      //Is n already known?
      {
        const double p = storage.Find(n);
        if (p!=0.0)
        {
          return p;
        }
      }

      //Check for simple phylogeny
      if (n.IsSimple())
      {
        const double p = n.CalcProbabilitySimpleNewick(theta);
        storage.Store(n,p);
        return p;
      }
      //Complex
      //Generate other Newicks and their coefficients
      std::vector<double> coefficients;
      std::vector<NewickType> newicks;
      {
        const double d = n.CalcDenominator(theta);
        typedef std::pair<std::vector<int>,int> NewickFrequencyPair;
        const std::vector<NewickFrequencyPair> newick_freqs
          = GetSimplerNewicksFrequencyPairs(n.Peek());
        for(const NewickFrequencyPair& p: newick_freqs)
        {
          const int frequency = p.second;
          assert(frequency > 0);
          if (frequency == 1)
          {
            newicks.push_back(p.first);
            coefficients.push_back(theta / d);
          }
          else
          {
            const double f_d = static_cast<double>(frequency);
            newicks.push_back(p.first);
            coefficients.push_back( (f_d*(f_d-1.0)) / d);
          }
        }
      }
      //Ask help about these new Newicks
      {
        const int sz = newicks.size();
        assert(newicks.size() == coefficients.size() );
        double p = 0.0;
        for (int i=0; i!=sz; ++i)
        {
          //Recursive function call
          p+=(coefficients[i] * CalculateProbability(newicks[i],theta,storage));
        }
        storage.Store(n,p);
        return p;
      }
    }
    catch (std::bad_alloc& e)
    {
      storage.CleanUp();
      std::cerr << "std::bad_alloc\n";
    }
    catch (std::exception& e)
    {
      storage.CleanUp();
      std::cerr << "std::exception";
    }
    catch (...)
    {
      storage.CleanUp();
      std::cerr << "Unknown exception\n";
    }
  }
}

} //~namespace newick
} //~namespace ribi

#endif // NEWICK_H


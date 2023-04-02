#ifndef COUNTLE_INCLUDED
#define COUNTLE_INCLUDED

#include "flint/fmpz_polyxx.h"
#include "flint/fmpzxx.h"
#include "htd/main.hpp"
#include "set.hpp"

#include "FeatureCalculation.cpp"
#include "Globals.cpp"
#include "Init.cpp"

extern std::vector<fixed_set<setSize>> subtreeVertices;
extern unsigned graphSize;
extern htd::ITreeDecomposition* ntd;
extern std::vector<std::unordered_map<int, int>> variableIndices;
extern std::vector<std::unordered_map<int, htd::vertex_t>> indexToVertex;
extern bool run;

extern int mat[][graphSizeLimit];

std::vector<std::vector<flint::fmpz_polyxx>> DPTable;
std::vector<std::unordered_map<htd::vertex_t, int>> bagExp;
std::vector<std::vector<std::unordered_map<htd::vertex_t, int>>>
    bagExpByPermutation;
std::vector<std::vector<slong>> exponentInBag;

// Helpers
inline void uniVariateToExpVector(htd::vertex_t node, slong uniVariateExp,
                                  int vec[]) {
  for (size_t ind = ntd->maximumBagSize() - 1; ind > 0; ind--) {
    if (indexToVertex[node][ind]) {
      vec[indexToVertex[node][ind]] = uniVariateExp / exponentInBag[node][ind];
      uniVariateExp = uniVariateExp % exponentInBag[node][ind];
    }
  }
  vec[indexToVertex[node][0]] = uniVariateExp;
}

inline slong calculateMonomialDegreeAndExpVector(htd::vertex_t node,
                                                 slong uniVariateExp,
                                                 int vec[]) {
  slong monomialDegree = 0;
  for (size_t ind = ntd->maximumBagSize() - 1; ind > 0; ind--) {
    if (indexToVertex[node][ind]) {
      vec[indexToVertex[node][ind]] = uniVariateExp / exponentInBag[node][ind];
      monomialDegree += vec[indexToVertex[node][ind]];
      uniVariateExp = uniVariateExp % exponentInBag[node][ind];
    }
  }
  vec[indexToVertex[node][0]] = uniVariateExp;
  monomialDegree += uniVariateExp;
  return monomialDegree;
}

inline slong expVectorTouniVariate(htd::vertex_t node, int vec[]) {
  slong ret = 0;
  for (htd::vertex_t bagVertex : ntd->bagContent(node)) {
    ret +=
        exponentInBag[node][variableIndices[node][bagVertex]] * vec[bagVertex];
  }
  return ret;
}

inline void formPermutationOfNumber(int x, int len,
                                    const flint::fmpzxx factorial[],
                                    int* perm) {
  std::vector<int> unused(len);
  for (int i = 0; i < len; i++) unused[i] = i;

  for (int i = 0; i < len; i++) {
    int count = x / factorial[len - i - 1].to<ulong>();
    perm[i] = unused[count];
    unused.erase(unused.begin() + count);
    x -= count * factorial[len - i - 1].to<ulong>();
  }
}

inline void formPermutationOfNumber(int x, int len, const int factorial[],
                                    int* perm) {
  std::vector<int> unused(len);
  for (int i = 0; i < len; i++) unused[i] = i;

  for (int i = 0; i < len; i++) {
    int count = x / factorial[len - i - 1];
    perm[i] = unused[count];
    unused.erase(unused.begin() + count);
    x -= count * factorial[len - i - 1];
  }
}

inline int permutationNumber(int perm[], const int len,
                             const flint::fmpzxx factorial[]) {
  int n = 0, ind = 0;
  std::vector<int> unused(len);
  for (int i = 0; i < len; i++) unused[i] = i;

  for (int i = 0; i < len; i++) {
    for (int j = 0; j < len; j++) {
      if (unused[j] == perm[i]) {
        ind = j;
        break;
      }
    }
    unused.erase(unused.begin() + ind);
    n += factorial[len - i - 1].to<ulong>() * ind;
  }
  return n;
}

inline int permutationNumber(int perm[], const int len, const int factorial[]) {
  int n = 0, ind = 0;
  std::vector<int> unused;
  unused.reserve(len);
  for (int i = 0; i < len; i++) {
    unused.push_back(i);
  }

  for (int i = 0; i < len; i++) {
    for (int j = 0; j < len; j++) {
      if (unused[j] == perm[i]) {
        ind = j;
        break;
      }
    }
    unused.erase(unused.begin() + ind);
    n += factorial[len - i - 1] * ind;
  }
  return n;
}

/*
Finds index of given vertex toFind in given bag
*/
inline int findBagIndex(htd::vertex_t vertex, htd::vertex_t toFind) {
  for (size_t i = 0; i < ntd->bagSize(vertex); i++) {
    if (ntd->bagContent(vertex)[i] == toFind) return i;
  }
  return indexNotFound;
}

inline void formBagIndices(htd::vertex_t vertex,
                           std::unordered_map<htd::vertex_t, int>& indices) {
  for (size_t i = 0; i < ntd->bagSize(vertex); i++) {
    indices[ntd->bagContent(vertex)[i]] = i;
  }
}

// assumes lower bound set to 0 and ub set to bagsize+1 beforehand
inline void formValidInterval(int* lowerBound, int* upperBound,
                              std::vector<htd::vertex_t>& before,
                              std::vector<htd::vertex_t>& after,
                              htd::vertex_t vertex, int* perm) {
  int index;
  for (htd::vertex_t adj : after) {
    index = findBagIndex(vertex, adj);
    if (index != indexNotFound)
      *upperBound = std::min(*upperBound, perm[index] + 1);
  }
  for (htd::vertex_t adj : before) {
    index = findBagIndex(vertex, adj);
    if (index != indexNotFound)
      *lowerBound = std::max(*lowerBound, perm[index] + 1);
  }
}

/*
Finds vertex at position index in permutation perm in bag of node vertex
*/
inline htd::vertex_t findFollower(htd::vertex_t vertex, int perm[], int index) {
  for (size_t i = 0; i < ntd->bagSize(vertex); i++) {
    if (perm[i] == index) return ntd->bagContent(vertex)[i];
  }
  return indexNotFound;
}

inline void calculateExponentsInBag(htd::vertex_t node) {
  slong prod = 1;
  int var;
  exponentInBag[node][0] = 1;
  for (std::size_t i = 1; i < ntd->maximumBagSize(); i++) {
    var = indexToVertex[node][i - 1];
    prod *= (bagExp[node][var] + 1);
    exponentInBag[node][i] = prod;
  }
}

/**
 * Converts polynomials encoded with Kronecker-substitution using exponents of
 * [fromNode] to polynomials with Kronecker-substitution exponents from bag
 * [toNode].
 */

inline void reEncodeExponents(htd::vertex_t fromNode, htd::vertex_t toNode) {
  int vec[graphSize + 1];
  flint::fmpz_polyxx temp;
  for (unsigned encodedMonomial = 0; encodedMonomial < DPTable[toNode].size();
       encodedMonomial++) {
    temp = flint::fmpz_polyxx();
    for (slong j = 0; j < DPTable[toNode][encodedMonomial].length(); j++) {
      if (DPTable[toNode][encodedMonomial].get_coeff(j) == 0) continue;
      uniVariateToExpVector(toNode, j, vec);
      slong newExp = expVectorTouniVariate(fromNode, vec);
      temp.set_coeff(newExp, DPTable[toNode][encodedMonomial].get_coeff(j));
    }
    DPTable[toNode][encodedMonomial] = flint::fmpz_polyxx(temp);
  }
}

// Countle main
// algorithm---------------------------------------------------------------------

// Leaf node
//---------------------------------------

inline void handleLeaf(htd::vertex_t node, flint::fmpzxx factorial[],
                       std::vector<htd::vertex_t> from[],
                       std::vector<htd::vertex_t> to[]) {
  subtreeVertices[node] = fixed_set<setSize>::empty();
  for (htd::vertex_t bagVertex : ntd->bagContent(node)) {
    subtreeVertices[node].set(bagVertex);
  }
  calculateExponentsInBag(node);
  if (leavesEmpty) {
    DPTable[node][0] = factorial[graphSize];
  } else {
    int bagSize = ntd->bagSize(node);
    int numPerms = factorial[bagSize].to<slong>();
    int perm[bagSize];
    int index, compIndex;
    bool valid;

    for (int i = 0; i < numPerms; i++) {
      valid = 1;
      formPermutationOfNumber(i, bagSize, factorial, perm);
      for (htd::vertex_t vertex : ntd->bagContent(node)) {
        index = findBagIndex(node, vertex);
        for (htd::vertex_t before : from[vertex]) {
          if (!valid) break;
          compIndex = findBagIndex(node, before);
          if (compIndex == indexNotFound) continue;
          if (perm[index] < perm[compIndex]) valid = 0;
        }
        for (htd::vertex_t after : to[vertex]) {
          if (!valid) break;
          compIndex = findBagIndex(node, after);
          if (compIndex == indexNotFound) continue;
          if (perm[index] > perm[compIndex]) valid = 0;
        }
      }
      if (valid) DPTable[node][i] = factorial[graphSize];
    }
  }
}

// Forget node
//------------------------------------------

/*
Removes vertex at bag index forgottenIndex from permtuation and moves everything
after it one step towards beginning of array to preserve bijectivity.
*/

inline void formPermutationWithoutForgotten(int perm[], int len,
                                            int forgottenIndex, int newPerm[]) {
  for (int i = 0, j = 0; i < len; i++) {
    if (i != forgottenIndex) {
      if (perm[i] > perm[forgottenIndex]) {
        newPerm[j++] = perm[i] - 1;
      } else {
        newPerm[j++] = perm[i];
      }
    }
  }
}

inline void updateBagExponentsForForgetNode(htd::vertex_t vertex,
                                            htd::vertex_t child,
                                            htd::vertex_t forgotten,
                                            flint::fmpzxx factorial[]) {
  int vec[graphSize + 1];
  std::size_t bagSize = ntd->bagSize(child);
  int perm[bagSize];
  int newPerm[bagSize - 1];
  std::unordered_map<htd::vertex_t, int> indices;
  formBagIndices(child, indices);

  for (std::size_t i = 0; i < DPTable[child].size(); i++) {
    if (DPTable[child][i].is_zero()) continue;
    formPermutationOfNumber(i, bagSize, factorial, perm);
    formPermutationWithoutForgotten(perm, bagSize, indices[forgotten], newPerm);
    int permNro = permutationNumber(newPerm, bagSize - 1, factorial);

    for (slong j = 0; j < DPTable[child][i].length(); j++) {
      if (DPTable[child][i].get_coeff(j) == 0) continue;
      uniVariateToExpVector(child, j, vec);

      for (htd::vertex_t bagVertex : ntd->bagContent(vertex)) {
        if (perm[indices[bagVertex]] == perm[indices[forgotten]] + 1) {
          bagExpByPermutation[vertex][permNro][bagVertex] =
              std::max(bagExpByPermutation[vertex][permNro][bagVertex],
                       vec[bagVertex] + vec[forgotten] + 1);
        } else {
          bagExpByPermutation[vertex][permNro][bagVertex] = std::max(
              bagExpByPermutation[vertex][permNro][bagVertex], vec[bagVertex]);
        }
      }
    }

    for (htd::vertex_t bagVertex : ntd->bagContent(vertex)) {
      bagExp[vertex][bagVertex] =
          std::max(bagExp[vertex][bagVertex],
                   bagExpByPermutation[vertex][permNro][bagVertex]);
    }
  }
}

inline void forgetBeforeLast(htd::vertex_t node, htd::vertex_t child,
                             int permNro, slong uniVariateExp,
                             const flint::fmpzxx& coeff,
                             htd::vertex_t successor, htd::vertex_t forgotten,
                             flint::fmpzxx factorial[]) {
  int expVector[graphSize + 1];
  uniVariateToExpVector(child, uniVariateExp, expVector);
  int successorExponent = expVector[successor];
  int forgottenExponent = expVector[forgotten];

  flint::fmpzxx coef =
      (coeff * factorial[forgottenExponent] * factorial[successorExponent] /
       factorial[successorExponent + forgottenExponent + 1])
          .evaluate();
  expVector[successor] = successorExponent + forgottenExponent + 1;
  slong newExp = expVectorTouniVariate(node, expVector);
  DPTable[node][permNro].set_coeff(
      newExp, coef + DPTable[node][permNro].get_coeff(newExp));
}

inline void forgetLast(htd::vertex_t node, htd::vertex_t child,
                       htd::vertex_t forgotten, int bagSize, int permNro,
                       slong uniVariateExp, const flint::fmpzxx& coeff,
                       flint::fmpzxx factorial[]) {
  int expVector[graphSize + 1];
  int monomialDegree =
      calculateMonomialDegreeAndExpVector(child, uniVariateExp, expVector);

  int hiddenAfterLast =
      subtreeVertices[child].popcount() - monomialDegree - bagSize;
  int forgottenExponent = expVector[forgotten];

  flint::fmpzxx coef =
      (coeff * factorial[forgottenExponent] * factorial[hiddenAfterLast] /
       factorial[hiddenAfterLast + forgottenExponent + 1])
          .evaluate();
  slong newExp = expVectorTouniVariate(node, expVector);
  DPTable[node][permNro].set_coeff(
      newExp, coef + DPTable[node][permNro].get_coeff(newExp));
}

inline void handleForget(htd::vertex_t vertex, flint::fmpzxx factorial[]) {
  htd::vertex_t child = ntd->childAtPosition(vertex, 0);
  htd::vertex_t forgotten = ntd->forgottenVertexAtPosition(vertex, 0);
  int forgottenIndex = findBagIndex(child, forgotten);
  int bagSize = ntd->bagSize(child);
  int perm[bagSize];
  int newPerm[bagSize - 1];

  updateBagExponentsForForgetNode(vertex, child, forgotten, factorial);
  calculateExponentsInBag(vertex);

  for (size_t i = 0; i < DPTable[child].size(); i++) {
    if (DPTable[child][i].is_zero()) continue;
    formPermutationOfNumber(i, bagSize, factorial, perm);
    formPermutationWithoutForgotten(perm, bagSize, forgottenIndex, newPerm);
    int permNro = permutationNumber(newPerm, bagSize - 1, factorial);
    int successor = findFollower(child, perm, perm[forgottenIndex] + 1);

    for (slong j = 0; j < DPTable[child][i].length(); j++) {
      if (DPTable[child][i].get_coeff(j) == 0) continue;
      if (successor != indexNotFound) {
        forgetBeforeLast(vertex, child, permNro, j,
                         DPTable[child][i].get_coeff(j).evaluate(), successor,
                         forgotten, factorial);
      } else {
        forgetLast(vertex, child, forgotten, bagSize, permNro, j,
                   DPTable[child][i].get_coeff(j).evaluate(), factorial);
      }
    }
  }
}

// Introduce node
//----------------------------------------------

/*
Form new permutation in array newPerm from old permutation perm with length len
by inserting vertex at index introducedIndex to rank where.
*/
inline void formPermutationAfterIntroduction(int perm[], int len,
                                             int introducedIndex, int where,
                                             int newPerm[]) {
  for (int i = 0, j = 0; i < len; i++) {
    if (i == introducedIndex) {
      newPerm[i] = where;
    } else if (perm[j] < where) {
      newPerm[i] = perm[j++];
    } else if (perm[j] >= where) {
      newPerm[i] = perm[j++] + 1;
    }
  }
}

inline void updateBagExponentsForIntroduction(htd::vertex_t cur,
                                              htd::vertex_t child,
                                              flint::fmpzxx factorial[],
                                              std::vector<htd::vertex_t>& from,
                                              std::vector<htd::vertex_t>& to,
                                              htd::vertex_t introduced) {
  int vec[graphSize + 1];
  int bagSize = ntd->bagSize(child);
  int childPerm[bagSize], newPerm[bagSize + 1];
  std::unordered_map<htd::vertex_t, int> indices;
  formBagIndices(cur, indices);
  int introducedIndex = indices[introduced];

  for (std::size_t childPermutation = 0;
       childPermutation < DPTable[child].size(); childPermutation++) {
    int lowerBound = 0;
    int upperBound = bagSize + 1;
    formPermutationOfNumber(childPermutation, bagSize, factorial, childPerm);
    formValidInterval(&lowerBound, &upperBound, from, to, child, childPerm);
    for (int introducedAt = lowerBound; introducedAt < upperBound;
         introducedAt++) {
      formPermutationAfterIntroduction(childPerm, bagSize + 1, introducedIndex,
                                       introducedAt, newPerm);
      int permNro = permutationNumber(newPerm, bagSize + 1, factorial);
      for (slong j = 0; j < DPTable[child][childPermutation].length(); j++) {
        if (DPTable[child][childPermutation].get_coeff(j) == 0) continue;

        int monomialDegree = calculateMonomialDegreeAndExpVector(child, j, vec);
        for (htd::vertex_t bagVertex : ntd->bagContent(child)) {
          bagExpByPermutation[cur][permNro][bagVertex] = std::max(
              bagExpByPermutation[cur][permNro][bagVertex], vec[bagVertex]);
          bagExp[cur][bagVertex] =
              std::max(bagExp[cur][bagVertex],
                       bagExpByPermutation[cur][permNro][bagVertex]);
          if (newPerm[indices[bagVertex]] == introducedAt + 1) {
            bagExpByPermutation[cur][permNro][introduced] = std::max(
                bagExpByPermutation[cur][permNro][introduced], vec[bagVertex]);
          }
        }
        if (introducedAt == bagSize) {
          bagExpByPermutation[cur][permNro][introduced] =
              std::max(bagExpByPermutation[cur][permNro][introduced],
                       (int)subtreeVertices[child].popcount() - monomialDegree -
                           bagSize);
        }
      }
      bagExp[cur][introduced] =
          std::max(bagExp[cur][introduced],
                   bagExpByPermutation[cur][permNro][introduced]);
    }
  }
}

inline void introduceVertexBeforeLast(int vertex, int child,
                                      slong uniVariateExp, flint::fmpzxx coeff,
                                      int introduced, int permNro,
                                      htd::vertex_t successor,
                                      flint::fmpzxx factorial[]) {
  int expVector[graphSize + 1];
  uniVariateToExpVector(child, uniVariateExp, expVector);
  int successorExponent = expVector[successor];
  for (int k = 0; k <= successorExponent; k++) {
    flint::fmpzxx coef =
        ((coeff * factorial[successorExponent] / factorial[k]) /
         factorial[successorExponent - k])
            .evaluate();
    expVector[introduced] = k;
    expVector[successor] = successorExponent - k;
    slong newExp = expVectorTouniVariate(vertex, expVector);
    DPTable[vertex][permNro].set_coeff(
        newExp, DPTable[vertex][permNro].get_coeff(newExp) + coef);
  }
}

inline void introduceVertexLast(htd::vertex_t vertex, htd::vertex_t child,
                                int bagSize, int introduced, int permNro,
                                slong uniVariateExp, flint::fmpzxx coeff,
                                flint::fmpzxx factorial[]) {
  int expVector[graphSize + 1];
  int monomialDegree =
      calculateMonomialDegreeAndExpVector(child, uniVariateExp, expVector);
  int hiddenAfterLast =
      subtreeVertices[child].popcount() - monomialDegree - bagSize;

  for (int k = 0; k <= hiddenAfterLast; k++) {
    flint::fmpzxx coef = ((coeff * factorial[hiddenAfterLast] / factorial[k]) /
                          factorial[hiddenAfterLast - k])
                             .evaluate();
    expVector[introduced] = k;
    slong newExp = expVectorTouniVariate(vertex, expVector);
    DPTable[vertex][permNro].set_coeff(
        newExp, DPTable[vertex][permNro].get_coeff(newExp) + coef);
  }
}

inline void handleIntroduction(htd::vertex_t vertex, flint::fmpzxx factorial[],
                               std::vector<htd::vertex_t>& from,
                               std::vector<htd::vertex_t>& to) {
  htd::vertex_t child = ntd->childAtPosition(vertex, 0);
  htd::vertex_t introduced = ntd->introducedVertices(vertex)[0];
  int bagSize = ntd->bagSize(child);
  int childPerm[bagSize];
  int newPerm[bagSize + 1];
  int introducedIndex = findBagIndex(vertex, introduced);
  subtreeVertices[vertex].set(introduced);
  int lb, ub;

  updateBagExponentsForIntroduction(vertex, child, factorial, from, to,
                                    introduced);
  calculateExponentsInBag(vertex);

  for (size_t i = 0; i < DPTable[child].size(); i++) {
    if (DPTable[child][i].is_zero()) continue;
    lb = 0;
    ub = bagSize + 1;
    formPermutationOfNumber(i, bagSize, factorial, childPerm);
    formValidInterval(&lb, &ub, from, to, child, childPerm);

    for (int j = lb; j < ub; j++) {
      formPermutationAfterIntroduction(childPerm, bagSize + 1, introducedIndex,
                                       j, newPerm);

      int successor = findFollower(child, childPerm, j);
      int permNro = permutationNumber(newPerm, bagSize + 1, factorial);

      if (successor != indexNotFound) {
        for (slong k = 0; k < DPTable[child][i].length(); k++) {
          if (DPTable[child][i].get_coeff(k).is_zero()) continue;
          introduceVertexBeforeLast(vertex, child, k,
                                    DPTable[child][i].get_coeff(k).evaluate(),
                                    introduced, permNro, successor, factorial);
        }
      } else {
        for (slong k = 0; k < DPTable[child][i].length(); k++) {
          if (DPTable[child][i].get_coeff(k).is_zero()) continue;
          introduceVertexLast(vertex, child, bagSize, introduced, permNro, k,
                              DPTable[child][i].get_coeff(k).evaluate(),
                              factorial);
        }
      }
    }
  }
}

// Join node
//---------------------------------------------

inline void updateBagExponentsForJoin(htd::vertex_t cur, htd::vertex_t y,
                                      htd::vertex_t z) {
  for (std::size_t i = 0; i < DPTable[cur].size(); i++) {
    if (DPTable[y][i].is_zero() || DPTable[z][i].is_zero()) continue;
    for (htd::vertex_t bagVertex : ntd->bagContent(cur)) {
      bagExpByPermutation[cur][i][bagVertex] =
          bagExpByPermutation[y][i][bagVertex] +
          bagExpByPermutation[z][i][bagVertex];
      bagExp[cur][bagVertex] = std::max(bagExp[cur][bagVertex],
                                        bagExpByPermutation[cur][i][bagVertex]);
    }
  }
}

inline void handleJoin(htd::vertex_t node, flint::fmpzxx factorial[]) {
  htd::vertex_t y = ntd->childAtPosition(node, 0);
  htd::vertex_t z = ntd->childAtPosition(node, 1);
  subtreeVertices[node] |= subtreeVertices[z];
  updateBagExponentsForJoin(node, y, z);
  calculateExponentsInBag(node);
  reEncodeExponents(node, y);
  reEncodeExponents(node, z);
  for (std::size_t i = 0; i < DPTable[node].size(); i++) {
    if (DPTable[y][i].is_zero() || DPTable[z][i].is_zero()) continue;
    DPTable[node][i] = DPTable[y][i] * DPTable[z][i] / factorial[graphSize];
  }
}

// Countle main

void runAlgorithm(flint::fmpzxx factorial[], std::vector<htd::vertex_t> from[],
                  std::vector<htd::vertex_t> to[]) {
  int numPerms;
  htd::PostOrderTreeTraversal traversal;
  traversal.traverse(*ntd, [&](htd::vertex_t vertex, htd::vertex_t parent,
                               size_t distanceToRoot) {
    numPerms = factorial[ntd->bagSize(vertex)].to<slong>();
    bagExpByPermutation[vertex].resize(numPerms);
    DPTable[vertex].resize(numPerms);

    if (!ntd->isLeaf(vertex)) {
      subtreeVertices[vertex] =
          fixed_set<setSize>(subtreeVertices[ntd->childAtPosition(vertex, 0)]);
    }

    if (ntd->isLeaf(vertex)) {
      handleLeaf(vertex, factorial, from, to);
    } else if (ntd->isForgetNode(vertex)) {
      handleForget(vertex, factorial);
    } else if (ntd->isIntroduceNode(vertex)) {
      htd::vertex_t introduced = ntd->introducedVertexAtPosition(vertex, 0);
      handleIntroduction(vertex, factorial, from[introduced], to[introduced]);
    } else if (ntd->isJoinNode(vertex)) {
      handleJoin(vertex, factorial);
    } else if (ntd->bagContent(vertex) ==
               ntd->bagContent(ntd->childAtPosition(vertex, 0))) {
      // Unnecessary extra node
      htd::vertex_t child = ntd->childAtPosition(vertex, 0);
      for (size_t i = 0; i < DPTable[child].size(); i++) {
        DPTable[vertex][i] = DPTable[child][i];
        bagExpByPermutation[vertex][i] = bagExpByPermutation[child][i];
        for (htd::vertex_t bv : ntd->bagContent(vertex)) {
          bagExpByPermutation[vertex][i][bv] =
              bagExpByPermutation[child][i][bv];
        }
      }
      bagExp[vertex] = bagExp[child];
      for (htd::vertex_t bv : ntd->bagContent(vertex)) {
        bagExp[vertex][bv] = bagExp[child][bv];
      }
      calculateExponentsInBag(vertex);
    } else {
      std::cout
          << "Error, given node is not one of leaf, forget, introduce or join!"
          << std::endl;
      std::cout << "Children's bags:" << std::endl;
      for (htd::vertex_t child : ntd->children(vertex)) {
        std::cout << "Child: " << child << std::endl;
        for (htd::vertex_t bv : ntd->bagContent(child)) {
          std::cout << bv << " ";
        }
        std::cout << std::endl;
      }
      return;
    }
  });
  flint::print(DPTable[ntd->root()][0].lead().evaluate());
  std::cout << std::endl;
  delete ntd;
}

int main(int argc, char* argv[]) {
  std::srand(0);
  htd::IMutableDirectedGraph* graph =
      manager->directedGraphFactory().createInstance();
  init(argc, argv, graph);

  if (!run) {
    delete graph;
    return 0;
  }
  int ntdSize = (int)ntd->vertexCount();

  setresourcelimit();

  subtreeVertices = std::vector<fixed_set<setSize>>(ntdSize + 1);
  DPTable = std::vector<std::vector<flint::fmpz_polyxx>>(
      ntdSize + 1, std::vector<flint::fmpz_polyxx>());
  exponentInBag = std::vector<std::vector<slong>>(
      ntdSize + 1, std::vector<slong>(ntd->maximumBagSize()));

  flint::fmpzxx flintorial[graphSize + 1];
  initializeFactorials(flintorial, graphSize);
  initializeExponents(graphSize);

  bagExp = std::vector<std::unordered_map<htd::vertex_t, int>>(
      ntdSize + 1, std::unordered_map<htd::vertex_t, int>());
  bagExpByPermutation =
      std::vector<std::vector<std::unordered_map<htd::vertex_t, int>>>(
          ntdSize + 1, std::vector<std::unordered_map<htd::vertex_t, int>>());
  std::vector<htd::vertex_t> from[graphSize + 1];
  std::vector<htd::vertex_t> to[graphSize + 1];

  for (htd::vertex_t i = 1; i <= graphSize; i++) {
    for (htd::vertex_t j = 1; j <= graphSize; j++) {
      if (mat[i][j]) {
        from[i].push_back(j);
        to[j].push_back(i);
      }
    }
  }

  if (run && ntd != nullptr) {
    runAlgorithm(flintorial, from, to);
  }
  delete graph;
  return 0;
}

#endif  // !COUNTLE_INCLUDED
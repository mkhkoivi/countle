#ifndef COUNTLEFEATCALC_INCLUDED
#define COUNTLEFEATCALC_INCLUDED

#include <algorithm>
#include <cmath>
#include <fstream>
#include <numeric>
#include <vector>

#include <htd/main.hpp>

#include "DefaultNormalizationParameters.cpp"
#include "Globals.cpp"

extern std::vector<double> featureMean;
extern std::vector<double> featureSTD;

extern std::vector<htd::ITreeDecomposition *> tdPool;
extern std::vector<double (*)(std::vector<double> &)> aggregates;
extern std::vector<double> results[poolSize];
extern std::string avgFile, stdFile;
extern bool fixedNormalizationParameters;
// Aggregates

inline double count(std::vector<double> &vec) {
  if (vec.empty()) return 0;
  std::sort(vec.begin(), vec.end());
  int count = vec.size();
  for (unsigned i = 0; i < vec.size() - 1; i++) {
    if (fabs(vec[i] - vec[i + 1]) < 1e-9) count--;
  }
  return (double)count;
}

inline double max(std::vector<double> &vec) {
  if (vec.empty()) return 0;
  double r = vec[0];
  for (size_t i = 0; i < vec.size(); i++) r = std::max(r, vec[i]);
  return r;
}

inline double median(std::vector<double> &vec) {
  if (vec.empty()) return 0;
  sort(vec.begin(), vec.end());
  return vec[vec.size() / 2];
}

inline double min(std::vector<double> &vec) {
  if (vec.empty()) return 0;
  double r = vec[0];
  for (size_t i = 0; i < vec.size(); i++) r = std::min(r, vec[i]);
  return r;
}

inline double sum(std::vector<double> &vec) {
  if (vec.empty()) return 0;
  double r = 0;
  for (size_t i = 0; i < vec.size(); i++) r += vec[i];
  return r;
}

inline double average(std::vector<double> &vec) {
  if (vec.empty()) return 0;
  return sum(vec) / vec.size();
}

inline double var(std::vector<double> &vec) {
  if (vec.empty()) return 0;
  double a = average(vec);
  double r = 0;
  for (size_t i = 0; i < vec.size(); i++) r += (vec[i] - a) * (vec[i] - a);
  return r / vec.size();
}

inline double sd(std::vector<double> &vec) { return std::sqrt(var(vec)); }

inline double separateVar(std::vector<double> &vec, double mean) {
  double r = 0;
  for (size_t i = 0; i < vec.size(); i++)
    r += (vec[i] - mean) * (vec[i] - mean);
  return r / vec.size();
}

inline double separateSD(std::vector<double> &vec, double mean) {
  return std::sqrt(separateVar(vec, mean));
}

// Features

/*
 *Computes features related to bag sizes in tree decomposition
 */

void computeAllBagSizeFeatures(const htd::ITreeDecomposition &decomposition,
                               unsigned tdIndex) {
  std::vector<double> nodes[7];
  for (htd::vertex_t node : decomposition.leaves())
    nodes[0].push_back(decomposition.bagSize(node));
  for (htd::vertex_t node : decomposition.forgetNodes())
    nodes[1].push_back(decomposition.bagSize(node));
  for (htd::vertex_t node : decomposition.introduceNodes())
    nodes[2].push_back(decomposition.bagSize(node));
  for (htd::vertex_t node : decomposition.joinNodes())
    nodes[3].push_back(decomposition.bagSize(node));
  for (htd::vertex_t node : decomposition.vertices()) {
    if (!decomposition.isLeaf(node))
      nodes[4].push_back(decomposition.bagSize(node));
  }
  for (htd::vertex_t node : decomposition.vertices()) {
    if (decomposition.bagSize(node) > 0)
      nodes[5].push_back(decomposition.bagSize(node));
  }
  for (htd::vertex_t node : decomposition.vertices())
    nodes[6].push_back(decomposition.bagSize(node));

  for (int i = 0; i < 6; i++) {
    results[tdIndex].push_back(aggregates[i](nodes[6]));
  }

  for (int j = 0; j < 6; j++) {
    for (int i = 0; i < 6; i++)
      results[tdIndex].push_back(aggregates[i](nodes[j]));
  }

  // Total bagsize sum counted from non-empty bags
  results[tdIndex].push_back(sum(nodes[5]));
  for (int i = 0; i < 4; i++) results[tdIndex].push_back(sum(nodes[i]));
}

void computeAllDepthFeatures(const htd::ITreeDecomposition &decomposition,
                             unsigned tdIndex) {
  std::vector<double> nodes[7];
  for (htd::vertex_t node : decomposition.leaves())
    nodes[0].push_back(decomposition.depth(node));
  for (htd::vertex_t node : decomposition.forgetNodes())
    nodes[1].push_back(decomposition.depth(node));
  for (htd::vertex_t node : decomposition.introduceNodes())
    nodes[2].push_back(decomposition.depth(node));
  for (htd::vertex_t node : decomposition.joinNodes())
    nodes[3].push_back(decomposition.depth(node));
  for (htd::vertex_t node : decomposition.vertices()) {
    if (!decomposition.isLeaf(node))
      nodes[4].push_back(decomposition.depth(node));
  }
  for (htd::vertex_t node : decomposition.vertices()) {
    if (decomposition.bagSize(node) > 0)
      nodes[5].push_back(decomposition.depth(node));
  }

  for (htd::vertex_t node : decomposition.vertices())
    nodes[6].push_back(decomposition.bagSize(node));

  for (int i = 0; i < 6; i++) {
    results[tdIndex].push_back(aggregates[i](nodes[6]));
  }

  for (int j = 0; j < 6; j++) {
    for (int i = 0; i < 6; i++)
      results[tdIndex].push_back(aggregates[i](nodes[j]));
  }
}

/*
Same as containerCount with aggregate average, count average amount of times a
vertex is met in tree decomposition
*/

double decompositionOverheadRatio(const htd::ITreeDecomposition &decomposition,
                                  const htd::IMultiHypergraph &graph) {
  std::vector<int> nodes(decomposition.vertexCount());
  for (htd::vertex_t vertex : decomposition.vertices())
    nodes.push_back(decomposition.bagSize(vertex));
  return std::accumulate(nodes.begin(), nodes.end(), 0) / graph.vertexCount();
}

/*
 *Computes features related to how many bags each vertex is present.
 */
void computeContainerCountFeatures(const htd::ITreeDecomposition &decomposition,
                                   unsigned tdIndex,
                                   const htd::IMultiHypergraph &graph) {
  std::vector<double> timesMet(graph.vertexCount(), 0);
  for (htd::vertex_t node : decomposition.vertices()) {
    for (htd::vertex_t vertex : decomposition.bagContent(node))
      timesMet[vertex - 1]++;
  }
  for (int i = 1; i < 6; i++)
    results[tdIndex].push_back(aggregates[i](timesMet));
}

/*
 *Computes features related to how many levels (depths) have bags containing
 *each vertex
 */

void computeItemLifeTimeFeatures(const htd::ITreeDecomposition &decomposition,
                                 unsigned tdIndex,
                                 const htd::IMultiHypergraph &graph) {
  std::vector<std::unordered_set<int>> vec(graph.vertexCount());
  for (htd::vertex_t node : decomposition.vertices()) {
    for (htd::vertex_t vertex : decomposition.bagContent(node))
      vec[vertex - 1].insert(decomposition.depth(node));
  }
  std::vector<double> vect(graph.vertexCount());
  for (unsigned i = 0; i < graph.vertexCount(); i++) vect[i] = vec[i].size();
  for (int i = 1; i < 6; i++) results[tdIndex].push_back(aggregates[i](vect));
}

/*
 *Computes percentage of given node types in tree decomposition.
 */

void computePercentages(const htd::ITreeDecomposition &decomposition,
                        unsigned tdIndex) {
  results[tdIndex].push_back(decomposition.leafCount() /
                             decomposition.vertexCount());
  results[tdIndex].push_back(decomposition.forgetNodeCount() /
                             decomposition.vertexCount());
  results[tdIndex].push_back(decomposition.introduceNodeCount() /
                             decomposition.vertexCount());
  results[tdIndex].push_back(decomposition.joinNodeCount() /
                             decomposition.vertexCount());
}

void computeSkipFunction(const htd::ITreeDecomposition &decomposition,
                         std::vector<std::vector<int>> &skipFun, int height,
                         int vertices) {
  // std::cout << "Height: " << decomposition.height() << std::endl;
  for (htd::vertex_t vertex : decomposition.vertices()) {
    // if (decomposition.isRoot(vertex))
    //   std::cout << vertex << " : " << decomposition.parent(vertex) << " : "
    //             << decomposition.depth(vertex) << std::endl;
    skipFun[0][vertex] = decomposition.parent(vertex);
  }

  for (int i = 0; i < height; i++) {
    for (htd::vertex_t vertex : decomposition.vertices()) {
      skipFun[i + 1][vertex] = skipFun[i][skipFun[i][vertex]];
    }
  }
}

void computeDistancesWithSkipFunction(
    const htd::ITreeDecomposition &decomposition, std::vector<double> &result,
    std::vector<std::vector<int>> &skipFun, int logHeight) {
  if (decomposition.joinNodeCount() < 2) return;
  int dist, exp, real;

  for (htd::vertex_t node1 : decomposition.joinNodes()) {
    for (htd::vertex_t node2 : decomposition.joinNodes()) {
      if (node1 >= node2) continue;
      htd::vertex_t a = node1, b = node2;
      if (decomposition.depth(a) < decomposition.depth(b)) std::swap(a, b);
      dist = decomposition.depth(a) - decomposition.depth(b);
      exp = logHeight;
      real = 1 << logHeight;
      while (decomposition.depth(a) > decomposition.depth(b)) {
        if (decomposition.depth(a) >= decomposition.depth(b) + real)
          a = skipFun[exp][a];
        exp--;
        real >>= 1;
      }

      exp = logHeight;
      real = 1 << logHeight;
      while (real > 0) {
        if (skipFun[exp][a] != skipFun[exp][b]) {
          a = skipFun[exp][a];
          b = skipFun[exp][b];
          dist += 2 * real;
        }
        real >>= 1;
        exp--;
      }

      if (a != b) dist += 2;
      result.push_back(dist);
    }
  }
}

void computeDistancesNaive(const htd::ITreeDecomposition &decomposition,
                           std::vector<double> &result) {
  for (htd::vertex_t node1 : decomposition.joinNodes()) {
    for (htd::vertex_t node2 : decomposition.joinNodes()) {
      if (node1 >= node2) continue;
      htd::vertex_t a = node1, b = node2;
      if (decomposition.depth(a) < decomposition.depth(b)) std::swap(a, b);
      int r = decomposition.depth(a) - decomposition.depth(b);
      while (decomposition.depth(a) > decomposition.depth(b))
        a = decomposition.parent(a);
      while (a != b) {
        r += 2;
        a = decomposition.parent(a);
        b = decomposition.parent(b);
      }
      result.push_back(r);
    }
  }
}

/*
O(n^3) solution, could be easily improved.
*/

void computeJoinNodeDistanceFeatures(
    const htd::ITreeDecomposition &decomposition, unsigned tdIndex) {
  std::vector<double> distances;

  if (decomposition.joinNodeCount() > 1) {
    if (decomposition.height() > 20) {
      int logHeight = std::log2(decomposition.height());
      int vertices = decomposition.vertexCount();
      std::vector<std::vector<int>> skip(logHeight + 1,
                                         std::vector<int>(vertices + 1));
      computeSkipFunction(decomposition, skip, logHeight, vertices);
      computeDistancesWithSkipFunction(decomposition, distances, skip,
                                       logHeight);
    } else
      computeDistancesNaive(decomposition, distances);
  }

  for (int i = 1; i < 6; i++)
    results[tdIndex].push_back(aggregates[i](distances));
}

void computeBranchingFactorFeatures(
    const htd::ITreeDecomposition &decomposition, unsigned tdIndex) {
  std::vector<double> vec;
  for (htd::vertex_t node : decomposition.vertices()) {
    vec.push_back(decomposition.childCount(node));
  }
  for (int i = 1; i < 6; i++) results[tdIndex].push_back(aggregates[i](vec));
}

/*
 *Naive O(t^3n) implementation, should be fast enough anyway. Computes ratio of
 *adjacent vertex pairs to total pairs in bag.
 */

void computeBagAdjacencyFeatures(const htd::ITreeDecomposition &decomposition,
                                 unsigned tdIndex,
                                 const htd::IMultiHypergraph &graph) {
  std::vector<double> vec;
  for (htd::vertex_t node : decomposition.vertices()) {
    int r = 0;
    for (htd::vertex_t vertex : decomposition.bagContent(node)) {
      for (htd::vertex_t vertex2 : decomposition.bagContent(node)) {
        if (vertex < vertex2) {
          r += (graph.isEdge(vertex, vertex2) || graph.isEdge(vertex2, vertex));
        }
      }
    }
    r *= 2;
    vec.push_back((double)r /
                  std::max(1ul, (decomposition.bagSize(node) *
                                 (decomposition.bagSize(node) - 1))));
  }
  for (int i = 1; i < 6; i++) results[tdIndex].push_back(aggregates[i](vec));
}

/*
 *Bag connectedness factor.
 */

void computeBagConnectednessFeatures(
    const htd::ITreeDecomposition &decomposition, unsigned tdIndex,
    const htd::IMultiHypergraph &graph) {
  std::vector<double> vec;

  for (htd::vertex_t node : decomposition.vertices()) {
    int r = 0;
    for (htd::vertex_t vertex : decomposition.bagContent(node)) {
      for (htd::vertex_t vertex2 : decomposition.bagContent(node)) {
        r += (vertex != vertex2 && graph.isConnected(vertex, vertex2));
      }
    }
    vec.push_back((double)r /
                  std::max(1ul, (decomposition.bagSize(node) *
                                 (decomposition.bagSize(node) - 1))));
  }
  for (int i = 1; i < 6; i++) results[tdIndex].push_back(aggregates[i](vec));
}

/*
 *Bag neighborhood coverage, ratio of neighbors in bag to neighbors in graph
 *averaged over vertices in bag.
 */

void computeBagNeighborhoodCoverageFeatures(
    const htd::ITreeDecomposition &decomposition, unsigned tdIndex,
    const htd::IMultiHypergraph &graph) {
  std::vector<double> vec;

  for (htd::vertex_t node : decomposition.vertices()) {
    std::vector<double> bagVec;
    for (htd::vertex_t vertex : decomposition.bagContent(node)) {
      int r = 0;
      for (htd::vertex_t vertex2 : decomposition.bagContent(node))
        r += (vertex != vertex2 && graph.isNeighbor(vertex, vertex2));
      bagVec.push_back((double)r / graph.neighborCount(vertex));
    }
    vec.push_back(average(bagVec));
  }
  for (int i = 1; i < 6; i++) results[tdIndex].push_back(aggregates[i](vec));
}

/*
 *For each bag and forgotten/introduced vertex counts how many neighbors are in
 *same bag before aggregate.
 */

void computeVertexNeighborCountFeatures(
    const htd::ITreeDecomposition &decomposition, unsigned tdIndex,
    const htd::IMultiHypergraph &graph) {
  std::vector<double> vec;
  htd::vertex_t vertex;

  for (htd::vertex_t node : decomposition.forgetNodes()) {
    vertex = decomposition.forgottenVertexAtPosition(node, 0);
    int r = 0;
    for (htd::vertex_t vertex2 : decomposition.bagContent(node))
      r += (vertex != vertex2 && graph.isNeighbor(vertex, vertex2));
    vec.push_back((double)r);
  }
  for (int i = 1; i < 6; i++) results[tdIndex].push_back(aggregates[i](vec));
  vec.clear();
  for (htd::vertex_t node : decomposition.introduceNodes()) {
    vertex = decomposition.introducedVertexAtPosition(node, 0);
    int r = 0;
    for (htd::vertex_t vertex2 : decomposition.bagContent(node))
      r += (vertex != vertex2 && graph.isNeighbor(vertex, vertex2));
    vec.push_back((double)r);
  }
  for (int i = 1; i < 6; i++) results[tdIndex].push_back(aggregates[i](vec));
}

/*
 *vertexConnectednessFactor measures the ratio of vertices connected to
 *introduced/forgotten node to all possible connections (pairs) of bag vertices.
 */

void computeVertexConnectednessFeatures(
    const htd::ITreeDecomposition &decomposition, unsigned tdIndex,
    const htd::IMultiHypergraph &graph) {
  std::vector<double> vec;
  double r;
  htd::vertex_t vertex;

  for (htd::vertex_t node : decomposition.forgetNodes()) {
    vertex = decomposition.forgottenVertexAtPosition(node, 0);
    r = 0;
    for (htd::vertex_t vertex2 : decomposition.bagContent(node))
      r += (vertex != vertex2 && graph.isConnected(vertex, vertex2));
    vec.push_back(r / std::max(1ul, (decomposition.bagSize(node) *
                                     (decomposition.bagSize(node) - 1))));
  }
  for (int i = 1; i < 6; i++) results[tdIndex].push_back(aggregates[i](vec));

  for (htd::vertex_t node : decomposition.introduceNodes()) {
    vertex = decomposition.introducedVertexAtPosition(node, 0);
    r = 0;
    for (htd::vertex_t vertex2 : decomposition.bagContent(node))
      r += (vertex != vertex2 && graph.isConnected(vertex, vertex2));
    vec.push_back(r / std::max(1ul, (decomposition.bagSize(node) *
                                     (decomposition.bagSize(node) - 1))));
  }
  for (int i = 1; i < 6; i++) results[tdIndex].push_back(aggregates[i](vec));
}

void featureCalculation(const unsigned tdIndex,
                        const htd::IMultiHypergraph &graph) {
  const htd::ITreeDecomposition *decomposition = tdPool[tdIndex];
  results[tdIndex].clear();

  computeAllDepthFeatures(*decomposition, tdIndex);
  computeAllBagSizeFeatures(*decomposition, tdIndex);

  computePercentages(*decomposition, tdIndex);

  computeContainerCountFeatures(*decomposition, tdIndex, graph);
  computeItemLifeTimeFeatures(*decomposition, tdIndex, graph);
  computeJoinNodeDistanceFeatures(*decomposition, tdIndex);
  computeBranchingFactorFeatures(*decomposition, tdIndex);
  computeBagAdjacencyFeatures(*decomposition, tdIndex, graph);
  computeBagConnectednessFeatures(*decomposition, tdIndex, graph);
  computeBagNeighborhoodCoverageFeatures(*decomposition, tdIndex, graph);
  computeVertexNeighborCountFeatures(*decomposition, tdIndex, graph);
  computeVertexConnectednessFeatures(*decomposition, tdIndex, graph);
}

// Evaluating feature-based selection accuracy and normalizing features.
inline bool normalizationParametersSet() {
  return avgFile != "" && stdFile != "";
}

// Implemented in IO
void readDoubleVecFromCSVFile(std::vector<double> &, std::string);

// If std is below 1e-15 correponding features set to 0 and ignored in
// calculation
inline void normalizeFeatures(unsigned tdCount) {
  std::vector<double> v;
  if (fixedNormalizationParameters && normalizationParametersSet) {
    readDoubleVecFromCSVFile(featureMean, avgFile);
    readDoubleVecFromCSVFile(featureSTD, stdFile);
  }

  // Last result is assumed to be actual runtime for comparision.
  for (unsigned i = 0; i < results[0].size() - 1; i++) {
    if (!fixedNormalizationParameters) {
      v.clear();
      for (unsigned j = 0; j < tdCount; j++) v.push_back(results[j][i]);
      featureMean[i] = average(v);
      featureSTD[i] = fabs(separateSD(v, featureMean[i]));
    }

    for (unsigned j = 0; j < tdCount; j++) {
      results[j][i] -= featureMean[i];
      if (featureSTD[i] > 1e-15)
        results[j][i] /= featureSTD[i];
      else
        results[j][i] = 0;
    }
  }
}

inline double denormalizedRuntime(double normalizedRuntime,
                                  unsigned runtimeIndex) {
  return (normalizedRuntime * featureSTD[runtimeIndex]) +
         featureMean[runtimeIndex];
}

// Note assumes run times are saved as extra feature with coefficient 0 given by
// evaluation!
void reportError(unsigned selected) {
  if (quiet) return;
  unsigned runtimeIndex = results[0].size() - 1;
  double selectedRuntime =
      denormalizedRuntime(results[selected].back(), runtimeIndex);
  int numBetter = 0;

  std::vector<double> runtimes(featuresToSelectFrom, 0);
  for (unsigned i = 0; i < featuresToSelectFrom; i++) {
    runtimes[i] = denormalizedRuntime(results[i].back(), runtimeIndex);
    if (runtimes[i] < selectedRuntime - 1e-9) {
      numBetter++;
    }
  }

  std::sort(runtimes.begin(), runtimes.end());

  // Order: selected   Best  Median  Worst   Average   #BetterChoices  #Choices
  std::cout << selectedRuntime << " " << runtimes[0] << ' '
            << runtimes[featuresToSelectFrom / 2] << ' ' << runtimes.back()
            << " " << featureMean[runtimeIndex] << ' ' << numBetter << " "
            << featuresToSelectFrom << "\n";
}

#endif  // !COUNTLEFEATCALC_INCLUDED

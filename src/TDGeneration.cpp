#ifndef COUNTLE_TDGENERATION_INCLUDED
#define COUNTLE_TDGENERATION_INCLUDED

#include <vector>

#include <htd/main.hpp>

#include "FeatureCalculation.cpp"
#include "Globals.cpp"
#include "IO.cpp"
#include "set.hpp"

extern std::vector<fixed_set<setSize>> subtreeVertices;
extern std::vector<htd::ITreeDecomposition *> tdPool;
extern bool readSeedsFromFile, sampleTD;
extern unsigned selected;
extern long long randSeed;
extern std::unique_ptr<htd::LibraryInstance> manager;
// Seeds could be made local when selected seed not required to be shown
extern int seeds[];

// Helpers

inline std::size_t subTreeSize(htd::vertex_t vertex) {
  return subtreeVertices[vertex].popcount();
}

inline void computeSubTrees(const htd::ITreeDecomposition &decomposition) {
  htd::PostOrderTreeTraversal traversal;
  subtreeVertices =
      std::vector<fixed_set<setSize>>(decomposition.vertexCount() + 1);

  traversal.traverse(
      decomposition, [&](htd::vertex_t vertex, htd::vertex_t parent,
                         std::size_t distanceToRoot) {
        if (decomposition.isLeaf(vertex)) {
          subtreeVertices[vertex] = fixed_set<setSize>::empty();
        } else {
          subtreeVertices[vertex] =
              subtreeVertices[decomposition.childAtPosition(vertex, 0)];
          if (decomposition.isIntroduceNode(vertex)) {
            subtreeVertices[vertex].set(
                decomposition.introducedVertexAtPosition(vertex, 0));
          } else if (decomposition.isJoinNode(vertex)) {
            subtreeVertices[vertex] |=
                subtreeVertices[decomposition.childAtPosition(vertex, 1)];
          }
        }
      });
}

// Heuristics

inline double joinNodeAvgDepth(const htd::ITreeDecomposition &decomposition) {
  std::vector<double> depths;

  for (htd::vertex_t node : decomposition.joinNodes()) {
    depths.push_back(decomposition.depth(node));
  }
  return average(depths);
}

inline double introduceNodeAvgDepth(
    const htd::ITreeDecomposition &decomposition) {
  std::vector<double> depths;

  for (htd::vertex_t node : decomposition.introduceNodes()) {
    depths.push_back(decomposition.depth(node));
  }
  return average(depths);
}

// Memory-based time approximations

inline double forgetTimeApprox(const htd::ITreeDecomposition &decomposition) {
  double s = 0;
  for (htd::vertex_t vertex : decomposition.forgetNodes()) {
    s += std::pow(subTreeSize(vertex), decomposition.bagSize(vertex));
  }
  return s;
}

inline double introduceTimeApprox(
    const htd::ITreeDecomposition &decomposition) {
  double s = 0;
  for (htd::vertex_t vertex : decomposition.introduceNodes()) {
    s += std::pow(subTreeSize(vertex), decomposition.bagSize(vertex));
  }
  return s;
}

inline double joinTimeApprox(const htd::ITreeDecomposition &decomposition) {
  double s = 0;
  int bagSize, popcount;
  for (htd::vertex_t vertex : decomposition.joinNodes()) {
    bagSize = decomposition.bagSize(vertex);
    popcount = subTreeSize(vertex);
    s += std::pow(popcount, bagSize) * bagSize * std::log(popcount);
  }
  return s;
}

inline double overallCostApproximation(
    const htd::ITreeDecomposition &decomposition) {
  computeSubTrees(decomposition);
  return forgetTimeApprox(decomposition) + introduceTimeApprox(decomposition) +
         joinTimeApprox(decomposition);
}

// Forming tree-decompositions

class FitnessFunction : public htd::ITreeDecompositionFitnessFunction {
 public:
  FitnessFunction(void) {}

  ~FitnessFunction() {}

  htd::FitnessEvaluation *fitness(
      const htd::IMultiHypergraph &graph,
      const htd::ITreeDecomposition &decomposition) const {
    HTD_UNUSED(graph)

    return new htd::FitnessEvaluation(2,
                                      -(double)(decomposition.maximumBagSize()),
                                      -overallCostApproximation(decomposition));
  }

  FitnessFunction *clone(void) const { return new FitnessFunction(); }
};

// Generates tree-decompositions with HTD.

void formTDPoolWithoutReset(htd::IMutableDirectedGraph *graph) {
  htd::TreeDecompositionOptimizationOperation *operation =
      new htd::TreeDecompositionOptimizationOperation(manager.get());
  operation->setManagementInstance(manager.get());
  operation->addManipulationOperation(new htd::NormalizationOperation(
      manager.get(), true, leavesEmpty, true, false));
  if (sampleTD)
    operation->setVertexSelectionStrategy(
        new htd::RandomVertexSelectionStrategy(10));

  htd::ITreeDecompositionAlgorithm *baseAlgorithm =
      manager->treeDecompositionAlgorithmFactory().createInstance();
  baseAlgorithm->addManipulationOperation(operation);

  FitnessFunction fitnessFunction;
  htd::IterativeImprovementTreeDecompositionAlgorithm algorithm(
      manager.get(), baseAlgorithm, fitnessFunction);
  algorithm.setIterationCount(poolSize);
  algorithm.setNonImprovementLimit(poolSize);

  htd::ITreeDecomposition *decomposition = algorithm.computeDecomposition(
      *graph, [&](const htd::IMultiHypergraph &graph,
                  const htd::ITreeDecomposition &decomposition,
                  const htd::FitnessEvaluation &fitness) {
        // Disable warnings concerning unused variables.
        HTD_UNUSED(graph)
        HTD_UNUSED(decomposition)

        tdPool.push_back(decomposition.clone());
      });
  HTD_UNUSED(decomposition)
}

void readSeeds(std::string path, int seeds[poolSize]) {
  std::ifstream in(path);
  if (!in.is_open()) {
    std::cout << "Couldn't open seeds file!" << std::endl;
    return;
  }
  for (int i = 0; i < poolSize; i++) {
    in >> seeds[i];
  }
  in.close();
}

inline void addTDToPool(int seed, htd::IMutableDirectedGraph *graph) {
  htd::TreeDecompositionOptimizationOperation *operation =
      new htd::TreeDecompositionOptimizationOperation(manager.get());
  operation->setManagementInstance(manager.get());
  operation->addManipulationOperation(new htd::NormalizationOperation(
      manager.get(), true, leavesEmpty, true, false));

  htd::ITreeDecompositionAlgorithm *baseAlgorithm =
      manager->treeDecompositionAlgorithmFactory().createInstance();
  baseAlgorithm->addManipulationOperation(operation);

  FitnessFunction fitnessFunction;
  htd::IterativeImprovementTreeDecompositionAlgorithm algorithm(
      manager.get(), baseAlgorithm, fitnessFunction);
  algorithm.setIterationCount(1);
  algorithm.setNonImprovementLimit(1);

  srand(seed);

  htd::ITreeDecomposition *decomposition = algorithm.computeDecomposition(
      *graph, [&](const htd::IMultiHypergraph &graph,
                  const htd::ITreeDecomposition &decomposition,
                  const htd::FitnessEvaluation &fitness) {
        // Disable warnings concerning unused variables.
        HTD_UNUSED(graph)
        HTD_UNUSED(decomposition)
        tdPool.push_back(decomposition.clone());
      });
  HTD_UNUSED(decomposition);
}

void formTDPoolByReset(htd::IMutableDirectedGraph *graph) {
  srand(randSeed);

  if (readSeedsFromFile)
    readSeeds("seeds.txt", seeds);
  else
    for (int i = 0; i < poolSize; i++) seeds[i] = rand();

  for (int i = 0; i < poolSize; i++) addTDToPool(seeds[i], graph);
}

// Select nice tree-decomposition for algorithm.

void selectTDFromPoolWithMemoryHeuristic() {
  double bestSoFar = 1e300;
  bool first = 1;
  for (unsigned i = 0; i < tdPool.size(); i++) {
    computeSubTrees(*tdPool[i]);
    double comp = overallCostApproximation(*tdPool[i]);
    if (std::isnan(comp)) {
      std::cout << "Comp: " << comp << " was skipped as it was nan"
                << std::endl;
      continue;
    }

    if (first || comp < bestSoFar) {
      selected = i;
      bestSoFar = comp;
      first = 0;
    }
  }
}

void selectTDFromPool(const htd::IMutableDirectedGraph &graph) {
  for (unsigned i = 0; i < tdPool.size(); i++) featureCalculation(i, graph);

  double bestVal = 1e300;
  bool first = 1;

  normalizeFeatures(poolSize);

  for (unsigned i = 0; i < tdPool.size(); i++) {
    double val = 0;
    for (unsigned j = 0; j < results[i].size(); j++)
      val += results[i][j] * multipliers[j];
    if (first || val < bestVal) {
      bestVal = val;
      selected = i;
      first = 0;
    }
  }
}

void selectBestFeaturesByFeatEval() {
  std::vector<double> realRuntimes, predictions, errors;
  for (unsigned i = 0; i < featuresToSelectFrom; i++) {
    realRuntimes.push_back(results[i].back());
  }
  normalizeFeatures(featuresToSelectFrom);
  double bestVal = 1e300;
  bool first = 1;

  for (unsigned i = 0; i < featuresToSelectFrom; i++) {
    double val = 0;
    for (unsigned j = 0; j < results[i].size(); j++)
      val += results[i][j] * multipliers[j];

    if (first || val < bestVal) {
      first = 0;
      bestVal = val;
      selected = i;
    }
    if (writePredictionsToFile) {
      predictions.push_back(multipliers[multipliersInUse - 1] + val);
      errors.push_back(fabs(predictions.back() - std::log(realRuntimes[i])));
    }
  }

  if (writePredictionsToFile) {
    std::ofstream out(predictionsToFile, std::ofstream::app);
    if (!out.is_open()) {
      std::cout << "Couldn't open output file!" << std::endl;
      return;
    }

    for (unsigned i = 0; i < featuresToSelectFrom; i++) {
      out << std::setprecision(3) << std::fixed << std::log(realRuntimes[i])
          << ' ' << predictions[i] << ' ' << errors[i] << std::endl;
    }
    out.close();
  }
}

htd::ITreeDecomposition *formNiceTreeDecomposition(
    htd::IMutableDirectedGraph *graph) {
  htd::TreeDecompositionOptimizationOperation *operation =
      new htd::TreeDecompositionOptimizationOperation(manager.get());
  operation->setManagementInstance(manager.get());
  operation->addManipulationOperation(new htd::NormalizationOperation(
      manager.get(), true, leavesEmpty, true, false));
  if (sampleTD)
    operation->setVertexSelectionStrategy(
        new htd::RandomVertexSelectionStrategy(10));

  htd::ITreeDecompositionAlgorithm *baseAlgorithm =
      manager->treeDecompositionAlgorithmFactory().createInstance();
  baseAlgorithm->addManipulationOperation(operation);

  FitnessFunction fitnessFunction;
  htd::IterativeImprovementTreeDecompositionAlgorithm algorithm(
      manager.get(), baseAlgorithm, fitnessFunction);
  algorithm.setIterationCount(200);
  algorithm.setNonImprovementLimit(100);
  if (sampleTD) algorithm.setIterationCount(1);
  srand(randSeed);

  htd::ITreeDecomposition *decomposition = algorithm.computeDecomposition(
      *graph, [&](const htd::IMultiHypergraph &graph,
                  const htd::ITreeDecomposition &decomposition,
                  const htd::FitnessEvaluation &fitness) {});
  tdPool.push_back(decomposition);
  return decomposition;
}

#endif  // !COUNTLE_TDGENERATION_INCLUDED
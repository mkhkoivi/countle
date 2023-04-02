#ifndef COUNTLEINIT_INCLUDED
#define COUNTLEINIT_INCLUDED

#include <stdlib.h>
#include <sys/resource.h>
#include <unistd.h>
#include <cstdlib>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include <flint/fmpzxx.h>
#include <htd/main.hpp>

#include "DefaultNormalizationParameters.cpp"
#include "IO.cpp"
#include "TDGeneration.cpp"

#define MEM_ENV_VAR "MAXMEM_MB"
#define TIME_ENV_VAR "MAXTIME"

extern std::vector<htd::ITreeDecomposition *> tdPool;
// Flags
extern bool sampleTD, fixedNormalizationParameters,
    featEval, features2file, header2file, pool, quiet, run,
    readSeedsFromFile, resetBetweenTDSelection,
    writePredictionsToFile;

extern std::string featurePath, evalCoefPath, predictionsToFile,
      avgFile, stdFile;
extern long long randSeed;
extern unsigned graphSize, selected;
extern htd::ITreeDecomposition *ntd;

extern slong exponent[maxTreeWidth + 1];

extern std::vector<std::unordered_map<int, int>> variableIndices;
extern std::vector<std::unordered_map<int, htd::vertex_t>> indexToVertex;

extern std::unique_ptr<htd::LibraryInstance> manager;

// Resource limits
void setresourcelimit() {
  struct rlimit memlimit, timelimit;
  long bytes;

  if (getenv(MEM_ENV_VAR) != NULL) {
    bytes = atol(getenv(MEM_ENV_VAR)) * (1024 * 1024);
    memlimit.rlim_cur = bytes;
    memlimit.rlim_max = bytes;
    setrlimit(RLIMIT_AS, &memlimit);
  }

  if (getenv(TIME_ENV_VAR) != NULL) {
    bytes = atol(getenv(TIME_ENV_VAR));
    timelimit.rlim_cur = bytes;
    timelimit.rlim_max = bytes;
    setrlimit(RLIMIT_CPU, &timelimit);
  }
}

// Note: Doesn't check if TD is valid for given graph! Usually seg fault follows
// from incorrect use though.
void parseOptions(int argc, char *argv[], std::string &path,
                  std::string &featurePath, std::string &tdReadPath,
                  std::string &tdWritePath) {
  int opt;
  opterr = 0;
  while ((opt = getopt(argc, argv, "a:A::B::e::f:F:g:hnpqRsS:x")) != -1) {
    switch (opt) {
      case 'a':
        writePredictionsToFile = 1;
        predictionsToFile = optarg;
        break;
      case 'A':
        fixedNormalizationParameters = 1;
        if (optarg) avgFile = optarg;
        break;
      case 'B':
        fixedNormalizationParameters = 1;
        if (optarg) stdFile = optarg;
        break;
      case 'e':
        featEval = 1;
        if (optarg) evalCoefPath = optarg;
        break;
      case 'f':
        features2file = 1;
        featurePath = optarg;
        break;
      case 'F':
        featurePath = optarg;
        break;
      case 'g':
        path = optarg;
        break;
      case 'h':
        header2file = 1;
        break;
      case 'n':
        run = 0;
        break;
      case 'p':
        pool = 1;
        break;
      case 'q':
        quiet = 1;
        break;
      case 'R':
        resetBetweenTDSelection = 1;
        break;
      case 's':
        sampleTD = 1;
        break;
      case 'S':
        std::srand(std::stoll(optarg));
        randSeed = std::stoll(optarg);
        break;
      case 'x':
        readSeedsFromFile = 1;
        break;
      case '?':
        std::cout << "Invalid argument!: " << std::endl;
        break;
    }
  }
}

void initializeExponents(size_t graphSize) {
  exponent[0] = 1;
  for (std::size_t i = 1; i < ntd->maximumBagSize(); i++) {
    exponent[i] = exponent[i - 1] * graphSize + 1;
  }
}

void initializeFactorials(flint::fmpzxx factorial[], size_t graphSize) {
  factorial[0] = flint::fmpzxx(1);
  if (graphSize < 1) return;
  factorial[1] = flint::fmpzxx(1);
  for (unsigned i = 2; i < graphSize + 1; i++) {
    factorial[i] = factorial[i - 1] * flint::fmpzxx(i);
  }
}

void chooseVariableIndices() {
  htd::PreOrderTreeTraversal traversal;
  for (size_t i = 0; i <= ntd->vertexCount(); i++) {
    variableIndices.push_back(std::unordered_map<int, int>());
    indexToVertex.push_back(std::unordered_map<int, htd::vertex_t>());
  }
  traversal.traverse(*ntd, [&](htd::vertex_t vertex, htd::vertex_t parent,
                               std::size_t distanceToRoot) {
    if (ntd->isRoot(vertex)) return;
    variableIndices[vertex] = variableIndices[parent];
    indexToVertex[vertex] = indexToVertex[parent];
    if (ntd->isIntroduceNode(parent)) {
      indexToVertex[vertex].erase(
          variableIndices[vertex][ntd->introducedVertexAtPosition(parent, 0)]);
      variableIndices[vertex].erase(ntd->introducedVertexAtPosition(parent, 0));
    } else if (ntd->isForgetNode(parent)) {
      std::vector<int> taken;
      for (auto p : variableIndices[vertex]) taken.push_back(p.second);
      std::sort(taken.begin(), taken.end());
      int x = 0;
      for (size_t i = 0; i < taken.size() && taken[i] == x; i++) x++;
      variableIndices[vertex][ntd->forgottenVertexAtPosition(parent, 0)] = x;
      indexToVertex[vertex][x] = ntd->forgottenVertexAtPosition(parent, 0);
    }
  });
}

void selectTreeDecomposition(std::string tdReadPath,
                             htd::IMutableDirectedGraph *graph) {
  std::chrono::high_resolution_clock::time_point before =
      std::chrono::high_resolution_clock::now();

  if (!quiet) std::cout << "Selecting tree decomposition." << std::endl;
  if (pool) {
    if (resetBetweenTDSelection)
      formTDPoolByReset(graph);
    else
      formTDPoolWithoutReset(graph);

    if (featEval)
      selectTDFromPool(*graph);
    else
      selectTDFromPoolWithMemoryHeuristic();
    ntd = tdPool[selected];
    computeSubTrees(*tdPool[selected]);
  } else
    ntd = formNiceTreeDecomposition(graph);
  if (!quiet) std::cout << "Tree decomposition has been selected!" << std::endl;
  std::chrono::high_resolution_clock::time_point after =
      std::chrono::high_resolution_clock::now();
  if (!quiet)
    std::cout << "Time spent computing tree decomposition was: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(after -
                                                                       before)
                     .count()
              << std::endl;
}

void init(int argc, char *argv[], htd::IMutableDirectedGraph *graph) {
  std::string path = "graphs2/gridtree_2_11.txt", tdReadPath = "td.txt",
              tdWritePath = "td.txt";
  aggregates = {count, max, average, median, min, sd};

  parseOptions(argc, argv, path, featurePath, tdReadPath, tdWritePath);
  if (featEval) {
    if (evalCoefPath != "") readFeatureEvaluationMultipliers();
  }
  if (header2file) writeHeader(featurePath);
  readGraph(graph, path);
  if (!quiet) std::cout << "Graph has been read." << std::endl;
  graphSize = (int)graph->vertexCount();
  selectTreeDecomposition(tdReadPath, graph);
  chooseVariableIndices();

  if (features2file) writeFeatures(selected, *graph, featurePath);
}

#endif  // !COUNTLEINIT_INCLUDED
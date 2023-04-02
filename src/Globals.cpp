#ifndef COUNTLE_GLOBALS_INCLUDED
#define COUNTLE_GLOBALS_INCLUDED

#include <flint/fmpzxx.h>
#include <htd/main.hpp>
#include "set.hpp"

// Determines whether leaf bag size is limited or not.
#define leavesEmpty 0
#define indexNotFound -1

#define graphSizeLimit 300
// SetSize has to be graphSizeLimit/64+1 or higher
#define setSize 6
#define maxTreeWidth 4
#define poolSize 100

// Flags
bool featEval = 0, features2file = 0, fixedNormalizationParameters = 0,
     header2file = 0, pool = 0, quiet = 0, run = 1, readSeedsFromFile = 0,
     resetBetweenTDSelection = 0, sampleTD = 0,
     writePredictionsToFile = 0;

int mat[graphSizeLimit][graphSizeLimit], seeds[poolSize];
unsigned featuresToSelectFrom = 0, multipliersInUse = 0, graphSize,
         selected = 0;
long long randSeed = 1;

std::string avgFile = "", stdFile = "", featurePath = "features.csv",
            evalCoefPath = "", path,
            predictionsToFile = "predictionErrors.txt";
htd::ITreeDecomposition *ntd;
slong exponent[maxTreeWidth + 1];

std::vector<double (*)(std::vector<double> &)> aggregates(6);
std::vector<double> results[poolSize];

std::vector<std::unordered_map<int, int>> variableIndices;
std::vector<std::unordered_map<int, htd::vertex_t>> indexToVertex;
std::vector<fixed_set<setSize>> subtreeVertices;
std::vector<htd::ITreeDecomposition *> tdPool;

std::unique_ptr<htd::LibraryInstance> manager(
    htd::createManagementInstance(htd::Id::FIRST));

#endif // !COUNTLE_GLOBALS_INCLUDED
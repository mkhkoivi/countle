#ifndef CountleIO_INCLUDED
#define CountleIO_INCLUDED

#include <fstream>
#include <iostream>
#include <istream>
#include <memory>
#include <stack>

#include <htd/main.hpp>

#include "DefaultNormalizationParameters.cpp"
#include "FeatureCalculation.cpp"
#include "Globals.cpp"
#include "TDGeneration.cpp"

extern std::vector<double> multipliers;

extern int mat[graphSizeLimit][graphSizeLimit];
extern unsigned featuresToSelectFrom, multipliersInUse;

/*
 *Assumes graphs is written in directed 0-1 adjacency matrix format with
 *integers separated by single space.
 */
void readGraph(htd::IMutableDirectedGraph *graph, std::string path) {
  std::ifstream infile(path);
  if (!infile.is_open()) {
    std::cout << "Couldn't open the file!" << std::endl;
    return;
  }
  std::string row;
  unsigned graphSize;

  if (std::getline(infile, row)) {
    unsigned n = row.length();
    graphSize = (n + 1) / 2;
    graph->addVertices(graphSize);
    if (row.compare("\n")) {
      for (unsigned i = 0; i < n; i += 2) {
        if (row[i] == '1') {
          graph->addEdge(1, i / 2 + 1);
          mat[1][i / 2 + 1] = 1;
        }
      }

      for (unsigned i = 2; i <= graphSize; i++) {
        std::getline(infile, row);
        for (unsigned j = 0; j < n; j += 2) {
          if (row[j] == '1') {
            graph->addEdge(i, j / 2 + 1);
            mat[i][j / 2 + 1] = 1;
          }
        }
      }
    }
  }
  infile.close();
}

class CSV : public std::string {};

std::istream &operator>>(std::istream &is, CSV &output) {
  std::getline(is, output, ',');
  return is;
}

void readFeatureEvaluationMultipliers() {
  std::ifstream in(evalCoefPath);
  if (!in.is_open()) {
    std::cout << "Couldn't open coefficient file!" << std::endl;
    return;
  }
  double d;
  while (in >> d) {
    multipliers[multipliersInUse++] = d;
  }
  in.close();
}

void writeHeader(std::string path) {
  std::ofstream out(path);
  if (!out.is_open()) {
    std::cout << "Couldn't open output file!" << std::endl;
    return;
  }
  std::string aggregates[] = {"Count\"",  "Maximum\"", "Average\"",
                              "Median\"", "Minimum\"", "StandardDeviation\""};
  std::string featureGroup1[] = {"\"NodeDepth",
                                 "\"LeafNodeDepth",
                                 "\"ForgetNodeDepth",
                                 "\"IntroduceNodeDepth",
                                 "\"JoinNodeDepth",
                                 "\"NonLeafDepth",
                                 "\"NonEmptyDepth",
                                 "\"BagSize",
                                 "\"LeafNodeBagSize",
                                 "\"ForgetNodeBagSize",
                                 "\"IntroduceNodeBagSize",
                                 "\"JoinNodeBagSize",
                                 "\"NonLeafBagSize",
                                 "\"NonEmptyBagSize"};
  std::string featureGroup2[] = {"\"ContainerCount",
                                 "\"ItemLifeTime",
                                 "\"JoinNodeDistance",
                                 "\"BranchingFactor",
                                 "\"BagAdjacencyFactor",
                                 "\"BagConnectednessFactor",
                                 "\"BagNeighborhoodCoverageFactor",
                                 "\"ForgottenVertexNeighborCount",
                                 "\"IntroducedVertexNeighborCount",
                                 "\"ForgottenVertexConnectednessFactor",
                                 "\"IntroducedVertexConnectednessFactor"};
  std::string specialFeatures[] = {
      "\"CumulativeBagSize\"",           "\"LeafNodeCumulativeBagSize\"",
      "\"ForgetNodeCumulativeBagSize\"", "\"IntroduceNodeCumulativeBagSize\"",
      "\"JoinNodeCumulativeBagSize\"",   "\"LeafPercentage\"",
      "\"ForgetNodePercentage\"",        "\"IntroduceNodePercentage\"",
      "\"JoinNodePercentage\""};
  std::string header = "";

  // Write depth and bagsize features
  for (size_t i = 0; i < 14; i++) {
    for (size_t j = 0; j < 6; j++) {
      header += featureGroup1[i] + aggregates[j] + ",";
    }
  }

  // Write cumulative and percentage features
  for (size_t i = 0; i < 9; i++) {
    header += specialFeatures[i] + ",";
  }

  // Write features without counting different values.
  for (size_t i = 0; i < 11; i++) {
    for (size_t j = 1; j < 6; j++) {
      header += featureGroup2[i] + aggregates[j] + ",";
    }
  }

  header += "\"User time\"\n";
  out << header;
  out.close();
}

// Writes all features already known at runtime, user time must be added after
// run is completed
void writeFeatures(unsigned tdIndex, const htd::IMutableDirectedGraph &graph,
                   std::string path) {
  char delim = ',';
  std::ofstream out(path, std::ofstream::app);
  if (!out.is_open()) {
    std::cout << "Couldn't open output file!" << std::endl;
    return;
  }
  featureCalculation(tdIndex, graph);
  for (size_t i = 0; i < results[tdIndex].size(); i++)
    out << results[tdIndex][i] << delim;
  out.close();
}

void readDoubleVecFromCSVFile(std::vector<double> &vec, std::string file) {
  std::ifstream in(file);
  if (!in.is_open()) {
    std::cout << "Couldn't open CSV file!" << std::endl;
    return;
  }
  std::string line;
  int i = 0;

  while (std::getline(in, line)) {
    std::istringstream iss(line);
    std::istream_iterator<CSV> it(iss);
    while (it != std::istream_iterator<CSV>()) {
      vec[i] = std::stod(*it);
      it++;
      i++;
    }
  }
  in.close();
}

#endif  // !CountleIO_INCLUDED
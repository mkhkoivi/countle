# countle
Counting the linear extensions of a poset using dynamic programming over a tree-decomposition

Implements the algorithms in Kustaa Kangas, Mikko Koivisto, Sami Salonen: A faster tree-decomposition based algorithm for counting linear extensions. Algorithmica 82(8): 2156–2173 (2020).

Dependencies:
  countle:
      -htd
      -flint

Compiling files:
  Run setup.sh. It downloads and installs htd and flint, and then runs compile.sh to compile countle.

  Note: compile-scripts may NOT have correct paths for your installation, so you may have to set them up yourself!

Running countle:
  Before running, make sure that graphSizeLimit and setSize are set correctly for graph you are running the algorithm on. GraphSizeLimit must be at least the size of your graph and setSize must be set to graphSizeLimit/64 +1.

Countle can be run as

./countle -g pathToGraph [optional arguments]

Optional:
          -a write predictions to file specified by optarg
          -A use fixedNormalizationParameters. File with avgs in CSV-format passed as optarg. If no optarg given default to fixed avg and std.
    -B use fixedNormalizationParameters. File with stdevs in CSV-format passed as optarg. If no optarg given default to fixed avg and std.
          -d to debug.
    -e to use feature evaluation instead of memory based heuristic. Requires path to file with space separated coefficients for features in order.
          -E write seeds used to generate TDS that weren't selected to file specified by optarg. Requires -p and -R or does nothing.
          -f to write features from tree decomposition used by algorithm, requires path to file as argument
    -F to simply select feature file without writing features (Can be used to only write header).
          -g to select graph (necessary or seg fault will occur)
    -h to write header to feature file
          -n if you don't want algorithm to run but only run setup (used for generating features/tds)
          -p to use tdpool instead of only keeping track of single TD and using memory heuristic.
          -q quiet-mode (don't output anything, not even answer)
          -r to read a td from file and run algorithm on it (still requires original graph to be specified!), requires path to tdfile as next argument
          -R to reset process between each td generation and only take first tds suggested by htd over and over.
          -s to sample tds, requires seed to be selected
    -S to specify seed for random number generator (required if -s chosen)
          -t to write selected tree decomposition to file, requires path as next argument
          -T to evaluate expected running time of features in file optarg and then print selected line number.

Known bugs:
        -Writing tds to file and reading them from file consecutively fails.

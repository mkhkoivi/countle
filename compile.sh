#!/usr/bin/env bash
g++ -std=c++11 -Wall -O3 src/Countle.cpp -o countle -Idependencies/htd/usr/local/include -Idependencies/flint/include -L./dependencies/htd/usr/local/lib -L$PWD/dependencies/flint/lib -lhtd -lgmp -lflint -Wl,-rpath,$PWD/dependencies/htd/usr/local/lib:$PWD/dependencies/flint/lib

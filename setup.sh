#!/usr/bin/env bash

mkdir dependencies
cd dependencies

echo "Downloading HTD v1.1.0"
wget https://github.com/mabseher/htd/archive/1.1.zip
unzip 1.1.zip

echo "Compiling and installing HTD v1.1.0 to $PWD"
mkdir htd
cd htd
cmake ../htd-1.1
make
make DESTDIR=$PWD install

cd ..
echo "Downloading flint-2.5.2"
wget http://www.flintlib.org/flint-2.5.2.tar.gz
mkdir flint
tar -xzvf flint-2.5.2.tar.gz

echo "Compiling flint v2.5.2"
cd flint-2.5.2
./configure --prefix=$(dirname $PWD)/flint
make

echo "Checking that flint works properly"
make check
echo "Installing flint"
make install

#Fixing a bad include
sed 's/flint.h/..\/flint.h/' -i ../flint/include/flint/flintxx/flint_classes.h

cd ../..
echo "Compiling countle"
./compile.sh
#!/bin/bash

cd ..

if [ ! -d RibiClasses ]; then
 git clone https://github.com/richelbilderbeek/RibiClasses
fi

if [ ! -d RibiLibraries ]; then
 git clone https://github.com/richelbilderbeek/RibiLibraries
fi

if [ ! -d TestBinaryNewickVector ]; then
 git clone https://github.com/richelbilderbeek/TestBinaryNewickVector
fi
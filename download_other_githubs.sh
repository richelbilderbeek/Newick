#!/bin/bash

cd ..

if [ ! -d RibiClasses ]; then
 git clone https://github.com/richelbilderbeek/RibiClasses
fi

if [ ! -d RibiLibraries ]; then
 git clone https://github.com/richelbilderbeek/RibiLibraries
fi

if [ ! -d BinaryNewickVector ]; then
  git clone https://github.com/richelbilderbeek/BinaryNewickVector
fi

if [ ! -d NewickVector ]; then
  git clone https://github.com/richelbilderbeek/NewickVector
fi

if [ ! -d ManyDigitNewick ]; then
  git clone https://github.com/richelbilderbeek/ManyDigitNewick
fi

if [ ! -d SortedBinaryNewickVector ]; then
  git clone https://github.com/richelbilderbeek/SortedBinaryNewickVector
fi

if [ ! -d TwoDigitNewick ]; then
  git clone https://github.com/richelbilderbeek/TwoDigitNewick
fi

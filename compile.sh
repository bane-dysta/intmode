#!/bin/bash

mkdir -p build

if [ "$1" == "static" ]; then   
    g++ -O2 -static -static-libstdc++ -static-libgcc -o build/intmode intmode.cpp
else
    g++ -o build/intmode intmode.cpp
fi

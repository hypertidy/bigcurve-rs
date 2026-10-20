#!/bin/sh
# Regenerate tests/fixtures from the C++ engine and PROJ.
# Needs a C++ compiler and PROJ (libproj-dev); run from the repo root.
set -e
cd "$(dirname "$0")/.."
c++ -std=c++17 -O2 -o oracle/oracle oracle/oracle.cpp $(pkg-config --cflags --libs proj)
rm -f tests/fixtures/*.txt
./oracle/oracle tests/fixtures

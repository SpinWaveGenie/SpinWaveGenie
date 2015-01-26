#!/bin/bash

brew tap homebrew/science
# many formulas are updated, so we suppress the output.
brew update > /dev/null 2>&1

brew install --quiet eigen
brew install --quiet nlopt
brew install --quiet tbb --c++11
#cmake and boost are already installed
#brew install --quiet cmake
#brew install --quiet boost --c++11


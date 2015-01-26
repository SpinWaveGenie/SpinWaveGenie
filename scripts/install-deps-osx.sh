#!/bin/bash
brew update

brew tap homebrew/science

brew install eigen
brew install nlopt
brew install tbb --c++11
#cmake and boost are already installed
#brew install cmake
#brew install boost --c++11


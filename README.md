
[![BSD License](https://img.shields.io/badge/license-BSD-blue.svg)](http://opensource.org/licenses/BSD-3-Clause)
[![Build Status](https://travis-ci.org/SpinWaveGenie/SpinWaveGenie.svg?branch=master)](https://travis-ci.org/SpinWaveGenie/SpinWaveGenie) 
[![Build status](https://ci.appveyor.com/api/projects/status/2m7m8u685l1vqk4u?svg=true)](https://ci.appveyor.com/project/quantumsteve/spinwavegenie)
<a href="https://scan.coverity.com/projects/4034">
  <img alt="Coverity Scan Build Status"
       src="https://scan.coverity.com/projects/4034/badge.svg"/>
</a>
[![Coverage Status](https://coveralls.io/repos/SpinWaveGenie/SpinWaveGenie/badge.svg)](https://coveralls.io/r/SpinWaveGenie/SpinWaveGenie)
[![Github Releases](https://img.shields.io/github/downloads/SpinWaveGenie/SpinWaveGenie/latest/total.svg)]()

---

# SpinWaveGenie
Library for simplifying linear spin wave calculations.

* **Performant**. SpinWaveGenie is written in C++ using many C++11 features. Linear algebra operations utilize the Eigen library and the code is parallelized over Q-points using the Intel Threading Building Blocks.
* **Extensible**. New interactions can easily be added by inheriting an abstract base class. Additional post-processing effects can be added via composition.
* **Post-processing**. Convolute your model calculation with a resolution function and/or integrate each data point over a region in reciprocal space.
* **Cross-platform**. Our continuous integration platforms build on Linux, OS X and Windows. A Homebrew formula and RPM package simplify user installation.
* **Free**. SpinWaveGenie and all of its dependencies are freely available and open source.

## Documenation

[User Installation Instructions](https://github.com/SpinWaveGenie/SpinWaveGenie/wiki/User-Installation-Instructions)

[Developer Installation Instructions](https://github.com/SpinWaveGenie/SpinWaveGenie/wiki/Installing-SpinWaveGenie)

[Example Scripts](https://github.com/SpinWaveGenie/SpinWaveGenie/wiki/Examples)

## Community

Public User Conversations

[![Gitter chat](https://badges.gitter.im/SpinWaveGenie/Users.svg)](https://gitter.im/SpinWaveGenie/Users "Gitter chat")

[Follow @SpinWaveGenie on twitter](https://twitter.com/SpinWaveGenie)

Contributor Conversations

[![Gitter chat](https://badges.gitter.im/SpinWaveGenie.svg)](https://gitter.im/SpinWaveGenie "Gitter chat")

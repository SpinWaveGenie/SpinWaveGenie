#How to contribute

We are really glad you are reading this and look forward to receiving 
your contribution to SpinWaveGenie! This document is 
meant to outline the steps for contributing to SpinWaveGenie. These 
are guidelines, not rules, so use your best judgment and feel free 
to propose changes to this document in a pull request.

##Submitting changes

We aspire to have similar guidelines as GitHub. In general, we ask that you

* Fork the repository.
* Make changes as you see fit.
* Submit a [pull request](https://github.com/blog/1943-how-to-write-the-perfect-pull-request) to this branch. 
* A pull request is a start to the conversation. Once the automated tests pass, we may suggest some changes, improvements or alternatives.

##Testing
Pull requests are automatically built using travis-ci on os x and linux. 
* Pull requests must build without warnings and all tests must pass on gcc 4.4+ and clang 3.4+
* All new work must contain automated tests. Code coverage is automatically submitted 
to coveralls and a 5% decrease in coverage fails the build. 
* Dynamic analysis including the address sanitizer and undefined behavior sanitizer are run on all tests.

##Coding conventions
* Run clang-format on all commits. We use LLVM style with some modifications described in [.clang-format](https://github.com/SpinWaveGenie/SpinWaveGenie/blob/master/.clang-format).
* Use doxygen comments in header files.

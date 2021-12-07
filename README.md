# cleanpt

`cleanpt` is a cleaner and modernized, in the C++ sense, version of smallpt, a small path tracer originally published by Kevin Beason [here](https://www.kevinbeason.com/smallpt/).

The goal of the original code was to be fit under 100 lines of code.

This goal lead to code that is hard to read, understand and reason about.
It also doesn't current C++ best practices.

For `cleanpt` I started from the original code and then incrementally improved it.
I do this as an exercice for myself and as a way to show people what can be done with "legacy" code to make it better.

Each big step will be registered as one commit to make it easier to track the progress.

# Build & run

You can build the project with CMake as usual.
You probably want to build in release as the program is computationally intensive and thus will complete much faster.
```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build .
./cleanpt 500 # samples per pixel, the higher the better, but the longer it takes
```

# TODO
 - [x] Set up Git
 - [x] Use CMake to build the code
 - [x] Write a simple readme to describe the project and how to build and run it
 - [x] Use clang-format to make the code more readable
 - [x] Set-up clang-tidy & fix what needs to be fixed
 - [x] Turn on most relevant warnings on & fix what needs to be fixed
 - [x] Add sanitizers support & fix what needs to be fixed
 - [x] Set-up CI to make sure the code build without warnings on multiple platforms
 - [ ] Start reading the code and refactoring it. I guess I'll need to:
   - [ ] Rename variables and/or functions to better describe their meaning
   - [ ] Comment the unintuitive parts
   - [ ] Break up the code into smaller functions
   - [ ] Make things `const` when possible
   - [ ] Look up for possible modernizations (I guess original code is C++98): `constexpr`, `std::array`, (parallel) algorithms, `enum class`, etc.
 - [ ] Split the code into a library + an application
 - [ ] Provide more customizations on the library: `double` vs `float`, output format, etc.
 - [ ] Add features to the application: better CLI, GUI with live update, etc.

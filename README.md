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


# Notes on refactoring

## First refactoring step (a.k.a bad comments)
When starting the refactoring, the first thing I looked for are comments.

Comments are here to clarify the intent of the code, but most of the time it is better to make the code as clear as possible (i.e self-documenting) and then removing the comment.

The issue with such comments is that they tend to get out of date when the code is changing.
They also sometimes try to legitimate some bad design.
A good example of such a comment here is:
```cpp
struct Vec {
    double x, y, z; // position, also color (r,g,b)
    ...
};
```
Here the type `Vec` is used to represent both positions in 3D space and colors in RGB format.
These two things have nothing in common and so should each have a distinct type.
Making strong types like this reduces the chance of making the mistake of mixing unrelated things since doing so would lead to a compiler error rather than a runtime bug.

## Const, const, const
Whenever possible, make a variable `const`.
This make the code more correct and less prone to refactoring bugs (oops modified the wrong variable).
It also helps the compiler optimizing the code since it can be sure that the variable's value will never change.

## Auto, auto, auto?
`auto` is a great tool as it can make the code more readable by making it shorter and remove the possibilities of implicit conversions among other things.

`auto` is especially great for things like:
```cpp
// ptr is std::unique_ptr<int>, writing the type brings nothing
auto ptr = std::make_unique<int>();

std::vector<int> vec{1, 2, 3};
// who wants to write std::vector<int>::iterator instead of auto?
auto it = vec.begin();
```

But since the type is not spelled out anymore, it is probably best to make even greater efforts to name your variables, functions, etc so that the code remains easily understandable.
```cpp
// What is x?
auto x = get_x();

// Naming can make things much clearer
auto distance = get_distance_to_objective();
```

In the end, when and when not to use `auto` is a matter of taste but keep in mind the above when choosing whether to use it or not for a piece of code.
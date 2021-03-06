set(WARNING_LEVEL DEFAULT CACHE STRING "Set the warning level for targets of this project. DEFAULT/MORE/ALL")
option(WARNINGS_AS_ERRORS "Treat warnings as errors" OFF)

function(add_warnings target)
    set(MSVC_DEFAULT_WARNINGS)
    set(MSVC_MORE_WARNINGS
        /W4 # Baseline reasonable warnings
    )
    set(MSVC_ALL_WARNINGS
        ${MSVC_MORE_WARNINGS}
        /w14242 # 'identifier': conversion from 'type1' to 'type1', possible loss of data
        /w14254 # 'operator': conversion from 'type1:field_bits' to 'type2:field_bits', possible loss of data
        /w14263 # 'function': member function does not override any base class virtual member function
        /w14265 # 'classname': class has virtual functions, but destructor is not virtual instances of this class may not
                # be destructed correctly
        /w14287 # 'operator': unsigned/negative constant mismatch
        /we4289 # nonstandard extension used: 'variable': loop control variable declared in the for-loop is used outside
                # the for-loop scope
        /w14296 # 'operator': expression is always 'boolean_value'
        /w14311 # 'variable': pointer truncation from 'type1' to 'type2'
        /w14545 # expression before comma evaluates to a function which is missing an argument list
        /w14546 # function call before comma missing argument list
        /w14547 # 'operator': operator before comma has no effect; expected operator with side-effect
        /w14549 # 'operator': operator before comma has no effect; did you intend 'operator'?
        /w14555 # expression has no effect; expected expression with side- effect
        /w14619 # pragma warning: there is no warning number 'number'
        /w14640 # Enable warning on thread un-safe static member initialization
        /w14826 # Conversion from 'type1' to 'type_2' is sign-extended. This may cause unexpected runtime behavior.
        /w14905 # wide string literal cast to 'LPSTR'
        /w14906 # string literal cast to 'LPWSTR'
        /w14928 # illegal copy-initialization; more than one user-defined conversion has been implicitly applied
        /permissive- # standards conformance mode for MSVC compiler.
    )
    set(MSVC_WARNINGS_AS_ERRORS
        /WX
    )

    set(CLANG_DEFAULT_WARNINGS)
    set(CLANG_MORE_WARNINGS
        -Wall
        -Wextra # reasonable and standard
    )
    set(CLANG_ALL_WARNINGS
        ${CLANG_MORE_WARNINGS}
        -Wshadow # warn the user if a variable declaration shadows one from a parent context
        -Wnon-virtual-dtor # warn the user if a class with virtual functions has a non-virtual destructor. This helps
                            # catch hard to track down memory errors
        -Wold-style-cast # warn for c-style casts
        -Wcast-align # warn for potential performance problem casts
        -Wunused # warn on anything being unused
        -Woverloaded-virtual # warn if you overload (not override) a virtual function
        -Wpedantic # warn if non-standard C++ is used
        -Wconversion # warn on type conversions that may lose data
        -Wsign-conversion # warn on sign conversions
        -Wnull-dereference # warn if a null dereference is detected
        -Wdouble-promotion # warn if float is implicit promoted to double
        -Wformat=2 # warn on security issues around functions that format output (ie printf)
    )
    set(CLANG_WARNINGS_AS_ERRORS
        -Werror
    )

    set(GCC_DEFAULT_WARNINGS)
    set(GCC_MORE_WARNINGS
        ${CLANG_MORE_WARNINGS}
    )
    set(GCC_ALL_WARNINGS
        ${CLANG_ALL_WARNINGS}
        -Wmisleading-indentation # warn if indentation implies blocks where blocks do not exist
        -Wduplicated-cond # warn if if / else chain has duplicated conditions
        -Wduplicated-branches # warn if if / else branches have duplicated code
        -Wlogical-op # warn about logical operations being used where bitwise were probably wanted
        -Wuseless-cast # warn if you perform a cast to the same type
    )
    set(GCC_WARNINGS_AS_ERRORS
        ${CLANG_WARNINGS_AS_ERRORS}
    )

    target_compile_options(
        ${target} PRIVATE
        $<$<CXX_COMPILER_ID:MSVC>:${MSVC_${WARNING_LEVEL}_WARNINGS}>
        $<$<CXX_COMPILER_ID:Clang>:${MSVC_${WARNING_LEVEL}_WARNINGS}>
        $<$<CXX_COMPILER_ID:GNU>:${GCC_${WARNING_LEVEL}_WARNINGS}>
    )

    if(WARNINGS_AS_ERRORS)
        target_compile_options(
            ${target} PRIVATE
            $<$<CXX_COMPILER_ID:MSVC>:${MSVC_WARNINGS_AS_ERRORS}>
            $<$<CXX_COMPILER_ID:Clang>:${CLANG_WARNINGS_AS_ERRORS}>
            $<$<CXX_COMPILER_ID:GNU>:${GCC_WARNINGS_AS_ERRORS}>
        )
    endif()
endfunction(add_warnings)

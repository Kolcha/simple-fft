Simple FFT library
==================

Yet another small and simple FFT implementation, provides both FFT of complex and real input data. It was written with focus on embedded application (for audio processing), but certainly can be used in any other cases.

Key features (attractive for embedded):

- `float` type instead of `double`

  not so many chips have FPU with `double` support, so `float` is reasonable choice to get the maximum performance

- no dynamic memory allocation

  library doesn't allocate any buffers, only operates on provided by the user, this makes possible to use it even in "bare-metal" environments

- in-place calculations

  library uses the same buffer that was provided as input for output, this allows to reduce memory usage

- very minimal binary footprint

  there are only few functions and dependency only on standard library math functions, so no extra "bloat"

- no run-time initialization required

  run-time initialization is fully optional, required data structures can be fully initialized with static data at compile time, see header for details

In addition ("nice-to-have"):

- no hidden static variables

  there are no static variables inside the library, so it can be freely used for simultaneous processing of multiple data streams in different configurations (different number of FFTs)

- full ANSI C99 compliance

  no compiler-specific extensions are used, only standard C99 features

- very simple but flexible interface

  API consists only of 3 straightforward functions, but you may use only one of them, see header for details

Performance is comparable with FFTW implementation, at least for the cases of 512/1024/2048 FFTs.

What's behind
-------------

This library doesn't introduce anything new, it is even not an "just implementation" of some FFT algorithm, this is actually mix of some parts of 2 existing libraries, just with obvious optimizations:

- [Arduino FFT](https://github.com/lloydroc/arduino_fft)

  the whole FFT algorithm implementation (for complex data input) is used from this library almost as is, just added lookup table to avoid heavy trigonometric functions calls

- [KISS FFT](https://github.com/mborgerding/kissfft)

  the way how FFT is done on real data input is copied from this library C++ implementation, just slightly improved, just to avoid heavy trigonometric functions calls

Any requirements, implications, or limitations of algorithms in the libraries above are applicable to this library.

Usage
-----

There is no build/project file provided, just build it as part of the project. API usage is straightforward, see the comments in header for details.

The library is "pure C", there are no C++ guards in the header file, so to use it in C++ environment just wrap include statement with `extern "C"`.

```c++
extern "C" {
#include "simple_fft.h"
}
```

License
-------

Public domain, I don't claim any copyright.

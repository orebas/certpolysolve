// Compiled by the tests ieee_arithmetic_header, fast_math_refused and finite_math_refused (cmake/Tests.cmake): the
// header that refuses a build which may rewrite the library's floating-point arithmetic
// (src/platform/ieee_arithmetic.hpp) compiles under the ordinary options and stops with its message under
// -ffast-math or -ffinite-math-only.
#include "platform/ieee_arithmetic.hpp"
int main() { return 0; }

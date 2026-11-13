# Documentation

The pages for users of CertPolySolve. Each is kept current with the library; the test suite checks their links and
their quoted programs. The records of the project's own work (dated reports, plans, the history) are not part of
this repository.

## Start here

- [The README](../README.md): what CertPolySolve does, how to install it, a first solve.
- [quickstart-cpp.md](quickstart-cpp.md): the quick start for C++: five short programs with their output, and what
  each kind of stop means.
- [examples-cpp.md](examples-cpp.md): eight more programs by topic, and the C++ traps to know.
- [how-it-works.md](how-it-works.md): how the solver works, in two pages, with the measured speed by degree.

## Reference

- [api.md](api.md): the guide to every type, option, entry point and result, with the C name beside each C++ name.
- [api-spec-0.5.md](api-spec-0.5.md): the specification of the 0.5 API: inputs, options, goals, results, errors,
  lifetimes and thread guarantees.
- [mathematics.md](mathematics.md): what each certificate proves.
- [production-methods.md](production-methods.md): the three methods and how the solver chooses among them.
- [resource-policy.md](resource-policy.md): the optional limits on exact arithmetic and storage.
- [testing.md](testing.md): what the test suite checks and how to run it.
- [glossary.md](glossary.md): the words the documentation uses, with the plain phrase that means the same and the
  name in the headers; and the notation of the mathematics page.

## About these pages

- [WRITING.md](WRITING.md): how the documentation is written: three readers, plain words, no residue of the
  project's work.
- [figures/](figures/make_figures.py): the pictures, drawn by a script from the solver's own output.
- [The changelog](../CHANGELOG.md).

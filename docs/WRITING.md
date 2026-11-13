# How the documentation is written

A style sheet for the pages a user reads: the README, the quick start, the examples, the guide, the how-it-works
page, the mathematics, the spec and the changelog. It adapts the writing programme used for the project's papers to
documentation. Internal notes (the handoff, the dated reports, the plans, all under `dev/`) are written for the
people and agents working on CertPolySolve and follow their own conventions; none of their vocabulary or history belongs
here, and `docs/` holds nothing but user pages.

## Who is reading

Three readers, in layers, and a page serves one of them:

1. **Someone who wants roots.** They do not want to know what a certificate is. The README and the quick start are
   theirs: plain words, one task per program, a result they can paste and run.
2. **Someone who wants to know what is proved and how.** The how-it-works page, the guide, the methods page and the
   mathematics are theirs: precise, defined before used, complete, still in plain words where plain words exist.
3. **A program reading on someone's behalf.** The guide, the headers, the spec and the tested examples are its
   sources, and `llms.txt` at the root says which files are authoritative and which are records.

A page written for layer 1 never needs a word from layer 2 to be understood. A page for layer 2 may assume the
reader has seen layer 1.

## The first-screen contract

By the end of the first screen of any page, the reader knows what the page is for, has seen one concrete thing (a
program, a number, a picture) and knows where to go next. The README's first screen says what CertPolySolve does, how
to install it, shows one solve and its output, and links onward. Background that does not serve one of those goes
further down or on another page.

## Examples before machinery

Show the program, then its output, then say what the reader should notice. A definition appears where the reader
first needs it, never in a glossary at the top. Every program in a user page is a complete file under `examples/`,
compiled and run by the test suite, and quoted verbatim (`tests/test_examples_quoted.py` checks). The output shown
is what the program printed.

## Words

- Plain words first: "a box proved to hold the root", not "a certified enclosure"; "the roots as numbers", not
  "the centres of the regions"; "stops and says why", not "an honest stop". The project's own terms (regions,
  certificates, the methods' names) appear in layer 2 pages, defined once, and in the glossary of the mathematics.
- The reader is "you". The library is "CertPolySolve" or "the solver"; "we" is for the authors' choices, used rarely.
- Say what the library does. Never argue for a design, defend a choice or mention what it replaced, unless the
  reader needs the explanation to use the library ("the brace form reads the constant term first; print the
  polynomial to see which you built").
- No history: no dates, version numbers of the past, commit hashes, gate numbers, session or agent names, or
  quotations from the project's discussions. A lint (`tests/test_doc_residue.py`) rejects these in user pages.
- Numbers come from measurements that exist, with the setting stated ("on one core, with the defaults").
  Asymptotic claims are measured ones unless a theorem is cited.
- Credit early: the libraries and the ideas CertPolySolve is built on are named on the README, by the people who
  published them.

## Sentences and paragraphs

Short sentences, one idea each. A paragraph has one job, and its first sentence says what. Few dashes and
semicolons. No "A, not B" constructions where "A" alone says it. No hedging ("may", "should", "generally") where
the behaviour is definite. Emphasis is rare: no boxed slogans, no bold runs; headings and code carry the structure.

## Pictures

A picture earns its place by showing something the text cannot say as well: the boxes around the roots of a
polynomial, a box holding two roots and then two boxes after more precision, time against degree. Each is generated
from real output by a script kept beside it, so it can be regenerated when the solver changes, and is checked in as
SVG.

## Before a page is published

- A careful nonexpert reads it and marks the first place they stop understanding; that place is the defect.
- A fresh reader with only the page in front of them, human or model, is asked to write a program from it; where
  they fail, the page is wrong, not the reader.
- The checklist of the writing programme (its section 9) is run on pages of layer 2.
- The tests pass: `doc_links`, `examples_quoted`, and the residue lint.

Linas' Fractals
===============
Analytic Combinatorics and Dynamical Systems
--------------------------------------------

This git repo contains a large collection of software and dozens of
formal papers exploring a range of concepts in self-similarity.
This is the code that was used to create the pictures found in
[**Linas Art Gallery**](https://linas.org/art-gallery), as well as the
text source for the PDF's in
[**Linas Math page**](https://linas.org/math/sl2z.html) (subtitled
**The Modular Group and Fractals**).

Most of the content is closely related to one-another, exploring
the relationships between the following topics:

* Continued fractions
* Rational numbers
* The Modular Group SL(2,Z)
* The Cantor Set
* The Gauss-Kuzmin-Wirsing (GKW) operator
* The Minkowski Question Mark function
* The Minkowski measure (the "derivative" of the Question Mark)
* Symmetries of period-doubling maps
* de Rham curves
* The Bernoulli Map
* The Beta Transform bx mod 1
* Subshifts
* Analytic Number Theory
* Newton Series
* Riemann Hypothesis
* Hurwitz zeta function
* Ramanujan identities
* The Divisor Operator
* Mode locking in the Circle Map

All of the above are deeply inter-related with one-another, and, in a
certain sense, are all "isomorphic", exhibiting different manifestations
of "the same thing": the Cantor set. There's so much here, mostly
because the Cantor set is a very rich system.

There's also some unrelated content in here, having to do with either
linguistics, or with quantum mechanics.

Table of Contents
-----------------
Some of the directories contain a lot of content, and are deeply nested;
others contain a lot, and are shallow. Yet others are nearly empty.
There's a fair amount of source code in here; however, *it was never
intended to be used/usable by others*. It is sloppy, hacky and
minimalist, just enough to do what I wanted, but not presentable to
the public. This is my private workshop, and not my showroom. I am
publishing this in response to requests that I make certain parts
available. Well, here's all of it.

FYI, some of the code is very old, and might not actually compile.
Some of the libraries and dependencies have moved, breaking this old
code. I will accept bug reports and pull requests to fix what's broken.

In alphabetical order:

- [**generate**](generate) -- Code for filling in 2D pixmaps of a
  broad assortment of dynamical systems. This includes code for making
  short animated movies of the same. The "main" workhorse/driver is
  `brat.C`, which will fill in a pixmap with floating-point numbers,
  given a range of parameters and a function implementing the dynamical
  system.  The shell script `cvt` will then assign colors to the floats
  and create an RGB image.

- [**generate/samples**](generate/samples) -- A variety of images that
  were created with the above. Most of these are test runs, and are
  incomplete or damaged in some way. You will find final versions of
  these in the art gallery.

- [**image**](image) -- A single utility for converting `flo` files to
  `mtv` files. A `flo` file contains one floating point number per
  pixel. An `mtv` file contains one 24-bit RGB per pixel. Both formats
  include the pixmap width and height, and no other metadata.

- [**lang**](lang) -- Some very simple natural-language and AI
  experiments from 1997. This is perhaps my earliest work on AI.
  The [diary.txt](lang/diary.txt) file writes up my thoughts from
  that era.

- [**misc**](misc) -- This is a ***HUGE*** directory, containing various
  experiments in various subdirectories. See the README file there, for
  a summary.

- [**paper**](paper) -- Another ***HUGE*** directory, containing LyX
  source for all of the papers that I've published, and a bunch of
  others that I have not published. Please keep in mind that using this
  stuff is plagiarism. I've been plagiarized before. If you think that
  you can just grab this stuff, and publish it under your own name,
  well, that would be wrong. You will get caught, sooner or later.
  The stuff here is a historical record, and is here for your
  entertainment. It's not here for theft.

---
title: The `mintest` test support library: The little engine that could
published: false
description: Presentation of the testing support library `mintest`
tags: 
# cover_image: https://direct_url_to_image.jpg
# Use a ratio of 100:42 for best results.
# published_at: 2024-07-18 15:13 +0000
---

The testing is an important element to improve the quality of a software. It's an old concept, almost as old as the software development itself. Fortunately, ever since [JUnit](https://martinfowler.com/bliki/Xunit.html) many testing frameworks and libraries have been proposed that alleviate the chores of writing software tests.

There are many software testing frameworks. Some of them are very complex, and in exchange to their testing power, they can be intrusive to the software being adapted to their use. That's a problem particularly for legacy projects.

But there's [minunit](https://github.com/siu/minunit/tree/master) by [David Siñuela Pastor](https://github.com/siu), aka *siu*, that uses a minimalistic philosophy, with clever and elegant solutions, that makes it easy to insert the testing practice to any software project.

A basic documentation is already present at the repository page, but here follows a more detailed documentation to hopefully explain it a little more to the busy programmers that maybe would like to have more information before trying to use the interesting `minunit` code.

## Architecture

All of `minunit` is contained in a single header file, named [minunit.h](https://github.com/siu/minunit/blob/master/minunit.h). It's cross-platform and can be used in Microsoft Windows, Linux and many other Unix varieties -- for instance MacOS.

`minunit` lists the CPU and the clock wall timings in high resolution, and since this is OS-specific, preprocessor directives invoke the right OS primitives to time the tests.

By the way, the preprocessor is heavily used in `minunit`. That will disgust C++ purists, but there are no cryptic preprocessor tricks and minunit code is very easy to read and understand. Since it works so well, usually there's no need to look under the hood to learn how the tests are done.

The `minunit` tests are integrated with the conventional program code, and there's not an external runner program that will run the tests.

## Units

Two elements of minunit that are present in most testing frameworks and libraries are:

* **Fixture** a function that will create a sane environment for the tests;

* **Teardown** a function that will release the data allocated by the fixture.

Both elements are used as function pointers and of them elements can be the `NULL` pointer.

The basic units of `minunit` are

1. *test report* is a simple function, actually a C preprocessor macro, that will summarize the timings 

2. *test suite* is the largest unit in `minunit`. It allows the definition the fixture and teardown function pointers, and can generate a report of the passed and failure test sunder it. It will **not** release the resources allocated by the fixture, by calling the teardown function, if one was provided -- this is a *test* responsability.

    Any number of tests can be associated with a suite, that's defined using the macro `MU_TEST_SUITE()` that annotates a test suite function.

    Inside the test suite function, the macro `MU_SUITE_CONFIGURE()` will set the fixture and teardown function pointers that will be used by each individual test -- see below.

    In the calling program -- possibly `main()` -- the test suite will be run using -- you got it -- `MU_RUN_SUITE()`, that will set to NULL the fixture and 


    A test function has three functions (actually C preprocessor) associated with it: ; `MU_RUN_SUITE()`, that will run all tests .

3. *test* is the

## A sample program

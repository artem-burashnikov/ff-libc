# ff-libc

[![MIT License][license-shield]][license-url]

## Overview

This is a C library which prorvides arithmetic operations in finite fields, also known as Galois fields.

### Project structure

```ignorelang
ff-libc/
├── src/
│   └── (source files)
│
├── tests/
│    └── (test files)
│
├── CMakeLists.txt
├── LICENSE
└── README.md
```

## Features

A field can be initialized by defining the underlying structure or by passing an irreducible polynomial along with a field characteristic to the initializer. 

Validation (valid irreducible polynomial, valid prime number) of arguments is not performed. 

Polynomial can be initialized seperately. Coefficients are stored in little-endian order.

- **Basic arithmetic operations modulo Irreducible polynomial:** 
    - GF_elem_sum
    - GF_elem_diff
    - GF_elem_prod
    - GF_elem_div
    - GF_elem_get_complement
    - GF_elem_get_neutral
    - GF_elem_get_unity
    - GF_elem_get_inverse


## Getting Started

### Prerequisites

- CMake (version 3.10 or higher)
- GCC

### Building

Open the terminal and follow these steps:

1. Clone the repository:

    ```sh
    git clone git@github.com:artem-burashnikov/ff-libc.git
    ```

2. Navigate to the project root:

    ```sh
    cd ff-libc
    ```

3. Generate Makefile using CMake:

    ```sh
    mkdir build && cd build;
    cmake ..
    ```

4. Build a static library:

    ```sh
    cmake --build . --clean-first
    ```

## Licenses

The project is licensed under an [MIT License][license-url].

<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[license-shield]: https://img.shields.io/github/license/artem-burashnikov/ff-libc.svg?style=for-the-badge&color=blue
[license-url]: LICENSE

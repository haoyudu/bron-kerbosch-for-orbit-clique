# Modified Bron-Kerbosch Algorithm for the Orbit-Clique problem

The code is a galois-based program computing the size of maximal skew sets of lines on the Fermat surface. 

## Dependencies

The program requires python 3.

The following python packages are needed to run the program:

- [numpy](https://numpy.org/)
- [galois](https://pypi.org/project/galois/)
- [sympy](https://www.sympy.org/en/index.html)

Note: installing galois will simultaneously install numpy because of its dependence.

## Usage

Run the program on command line as:

```$ python3 LinesOnFermat```
(or `$ python LinesOnFermat` depending on your system)

Users will be prompted to enter the chosen characteristic $p$ and power $e$ of the field over which the Fermat surface of degree $q = p^e$ will be defined. The inputs $p$ and $e$ should be integers.

The efficiency of the code supports $q \leq 4$.

The example has hard-coded the choice of lines to be the indices $0,q+2,2q+4$.

Included in the files under `q=4 results`, `LineAuto4.py` contains automorphisms generated for $q = 4$ and `graph4.py` contains the corresponding graphs.

## Output

The program will output the following:

- There are at most $n$ orbits of maximal skew sets of size $m$ containing the set [list of indices] up to the action of the stablizer of that skew set.

- There are $k$ maximal skew sets of size $l$ containing the set [list of indices].

## Runtime

Currently, we expect the code to finish within a few minutes for $q = 3$. The runtime depends on number of cores used for $q \geq 4$.

## Example

By running main, the user should be able to achieve these outputs (note that 2 and 1 are user inputs, which the user must give to the program when prompted):

```
    What characteristic are you working in?
    2
    Which power of it do you wish to do use?
    1
    There are at most 1 orbits of maximal skew sets of size 5 containing the set [0, 4, 8] up to the action of the stablizer of that skew set.
    There are at most 1 orbits of maximal skew sets of size 6 containing the set [0, 4, 8] up to the action of the stablizer of that skew set.
    There are 3 maximal skew sets of size 5 containing the set [0, 4, 8]
    There are 2 maximal skew sets of size 6 containing the set [0, 4, 8]
```

## Contributors

Haoyu Du, Madhav Krishna, Janet Page, Timothy Ryan

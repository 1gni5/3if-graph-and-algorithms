# Simulated annealing applied to N-Queens

* Compile with `make nq`
* Execute with `./nq <input`

Generates a plot (`plot.ps`) of the cost against the temperature

Input must have the following format:

```
ALPHA = XXX
BETA = XXX
BETA0 = XXX
MAXTIME = XXX
Title of the first n-queens problem to solve
10
Title of the second n-queens problem to solve
20
END
```

* The first lines are the values of the parameters for the 
  simulated annealing strategy.
* Afterwards, each puzzle is described by two consecutive lines. 
  The first one is a title.
* The second one gives the number of queens.
  It is made of numbers separated by spaces.
  * The first number is the total number of queens.
  * Initially, queen i is on row i and column i.
  * A configuration is represented such that queen i stays on column i.
  * If a move swaps the rows of two queens, then queens can only be attacked on 
    diagonals. This property is used in the code.

# Paramless

## What is it?

We use a framework to evolve function valued traits. This is based off [the following code in python.](https://github.com/juliangarcia/paramless)

## Command-line useage

There is a simplified parameter input method which may be called as follows:

```
julia journal_paramless_simulation.jl version N c k quality iterations snapshot_iter
```

- $version \in \{3,4\}$
- $N \in \mathbb{N}$, it is the number of points $[0,1]$ gets discretised as 
- $c \in [0,1]$, cost to a scientist to submit
- $k \in [0,1]$, penalises harsh confusion curves for a journal
- quality is either true or false, false means cost is the rejection rate
- iterations is the number of iterations the simulation will run for 
- snapshot_iter will record the submission and confusion curve at the specified period

This will use a seed of $0$.

To specify the seed, add:

```
julia journal_paramless_simulation.jl version N c k quality iterations snapshot_iter seed
```

To specify mutation parameters, add the following four parameters:

```
julia journal_paramless_simulation.jl version N c k quality iterations snapshot_iter seed sci_epsilon journal_epsilon width journal_width
```

Here, the epsilon parameters indicate the size of the bump. By default, it is $0.0001$. The width indicates how far a noticeable gaussian mutation should be applied. By default, it is $0.001$

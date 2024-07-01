# Community Structure Testing by Counting Frequent Common Neighbor Sets

This repository contains the MATLAB implementation of the algorithm presented in the paper "Community Structure Testing by Counting Frequent Common Neighbor Sets". The code evaluates the community structure of a network by calculating a p-value using a Poisson distribution model based on common neighbor sets.

## Files Description

- `run.m`: This script is the main entry point to run the analysis.
- `pvalue_poisson.m`: Function to calculate the p-value given the test statistic and the parameter of the Poisson distribution.
- `pX.m`: Function to compute the probability as described in the paper by Kirsch et al. (2012).
- `RunPoisson.m`: Core function that computes the p-value by analyzing the network's adjacency matrix.
- `network_example.txt`: Sample dataset containing the network's adjacency matrix.

## How to Run

1. Clone this repository to your local machine.
2. Ensure you have MATLAB installed.
3. Open MATLAB and navigate to the directory containing these files.
4. Run the `run.m` script to execute the analysis.

## Example Usage

```matlab
run.m

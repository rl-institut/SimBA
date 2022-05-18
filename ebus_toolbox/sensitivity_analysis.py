""" Optimizer that evaluates inputs and outputs of every iteration, adapting scenario setup
    to optimize for specified metrics.
"""
from ebus_toolbox.util import read_arguments

def no_optimization():
    args = read_arguments()
    return "converged"

if __name__ == '__main__':
    no_optimization()

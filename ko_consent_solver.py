import random

import gurobipy as gp
import matplotlib.pyplot as plt
import numpy as np
from gurobipy import GRB
from copy import copy

m = gp.Model("simple_net")


class net_dynamics:
    def __init__(self, A, delta):
        assert A.shape[0] == A.shape[1]
        assert np.all(np.diag(A) == 0)
        self.nnodes = A.shape[0]
        self.dynA = copy(A)
        self.dynA[range(self.nnodes), range(self.nnodes)] = -np.sum(A, axis=1)
        self.dynD = self.dynA @ delta

    def __call__(self, x):
        x_dot = self.dynA @ x + self.dynD
        return x_dot


def run_knocked_out_network(system: net_dynamics, x0: np.ndarray, timestep, threshhold, knock_out_list):
    x0 = x0.reshape(-1, 1)
    x0[knock_out_list, :] = 0
    x = copy(x0)
    trace = copy(x0)
    while True:
        x_dot = system(x)
        x_dot[knock_out_list, :] = 0
        if np.all(np.abs(x_dot) < threshhold):
            break
        x += x_dot * timestep
        trace = np.concatenate((trace, x), axis=1)
    return trace


def get_sparse_network(nnodes, sparsity):
    A = np.random.random((nnodes, nnodes))
    A = (A + A.T) / 2  # making graph undirected
    A = A - np.diag(np.diag(A))  # setting diagonal to zero
    d = np.random.random((nnodes, 1))
    set_to_zero = np.random.random((nnodes, nnodes)) < sparsity
    set_to_zero = np.maximum(set_to_zero, set_to_zero.T)
    A[set_to_zero] = 0
    return A, d


nnodes = 5
A, d = get_sparse_network(nnodes=nnodes, sparsity=0.4)
network_dynamics = net_dynamics(A, d)
knock_outs = [[i, j] for i in range(nnodes) for j in range(nnodes) if i!=j]
# knock_outs.append([])
knock_outs = knock_outs[0:6]
experiment_numbers = np.zeros((0, ))
row_numbers = np.zeros((0, ))

X = np.zeros((nnodes, len(knock_outs)))

# Running the experiments with / without knockouts
for k in range(len(knock_outs)):
    knock_out_list = knock_outs[k]
    trace = run_knocked_out_network(network_dynamics, x0=np.zeros((nnodes, 1)), timestep=0.01, threshhold=0.001,
                                    knock_out_list=knock_out_list)
    X[:, k] = trace[:, -1]
    converged_nodes = set(range(nnodes)) - set(knock_out_list)
    converged_nodes = np.array(list(converged_nodes))
    row_numbers = np.concatenate((row_numbers, converged_nodes), axis=0)
    experiment_numbers = np.concatenate((experiment_numbers, k * np.ones(converged_nodes.shape)), axis=0)
    if len(knock_out_list) == 0:
        # d_tilde = trace[:, -1]
        plt.plot(trace.T)
        plt.show()
d_tilde = network_dynamics.dynD
row_numbers = row_numbers.astype(int)
experiment_numbers = experiment_numbers.astype(int)

plt.plot(range(nnodes), X, '*')
plt.show()
# Create variables
# d = m.addMVar(nnodes, vtype=GRB.CONTINUOUS)
a = m.addMVar(shape=(nnodes, nnodes), vtype=GRB.CONTINUOUS, lb=-np.inf, name='A')
err = m.addMVar(len(experiment_numbers), vtype=GRB.CONTINUOUS)


# Add constraint: sum of the rows is zero, off-diagonal elements are positive
for i in range(nnodes):
    vect = np.ones((1, nnodes))
    m.addConstr(vect @ a[i, :] == np.zeros((1, 1)))
    for j in range(nnodes):
        if i != j:
            m.addConstr(a[i, j] >= 0.0)
        if i < j:
            m.addConstr(a[i, j] == a[j, i])

# Add constraint: error is the difference
for i in range(len(experiment_numbers)):
    m.addConstr(a[row_numbers[i], :] @ X[:, experiment_numbers[i]] + d_tilde[row_numbers[i]] - err[i] <= 0)
    m.addConstr(a[row_numbers[i], :] @ X[:, experiment_numbers[i]] + d_tilde[row_numbers[i]] + err[i] >= 0)



# Set objective
# m.setObjective(err.sum() + , GRB.MINIMIZE)
rho = 0.0
# m.setObjective(-rho*gp.quicksum([a[i, i] for i in range(nnodes)]), GRB.MINIMIZE)
m.setObjective(gp.quicksum(err)-rho*gp.quicksum([a[i, i] for i in range(nnodes)]))
# Optimize model
m.optimize()

print('hi')

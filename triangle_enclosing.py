
import numpy as np
import gurobipy as gp
from gurobipy import GRB

def triangle_enclosing(A, L = 5, time_limit = 300):

    m = A.shape[0]

    lambda_indices = []

    for i in range(m):
        lambda_indices.append(i)

    alpha_minus_indices = lambda_indices
    alpha_plus_indices = lambda_indices

    beta_indices = lambda_indices

    lambda_prime_indices = lambda_indices

    alpha_prime_minus_indices = lambda_indices
    alpha_prime_plus_indices = lambda_indices

    beta_prime_indices = lambda_indices

    gamma_indices = lambda_indices

    dist_indices = []

    for i in range(3):
        dist_indices.append(i)

    MODEL = gp.Model('triangle_enclosing')

    landa = MODEL.addVars(lambda_indices, vtype = GRB.CONTINUOUS, lb = -10, ub = 10, name = 'lambda')
    w_landa = MODEL.addVars(lambda_indices, 2, vtype = GRB.CONTINUOUS, lb = -10, ub = 10, name = 'w_lambda')
    alpha_minus = MODEL.addVars(alpha_minus_indices, vtype = GRB.BINARY, name = 'alpha_minus')
    alpha_plus = MODEL.addVars(alpha_plus_indices, vtype = GRB.BINARY, name = "alpha_plus")
    beta = MODEL.addVars(beta_indices, vtype = GRB.BINARY, name = 'beta')

    landa_prime = MODEL.addVars(lambda_prime_indices, vtype = GRB.CONTINUOUS, lb = -10, ub = 10, name = 'lambda_prime')
    w_landa_prime = MODEL.addVars(lambda_prime_indices, 2, vtype = GRB.CONTINUOUS, lb = -10, ub = 10, name = 'w_lambda_prime')
    alpha_prime_minus = MODEL.addVars(alpha_prime_minus_indices, vtype = GRB.BINARY, name = 'alpha_prime_minus')
    alpha_prime_plus = MODEL.addVars(alpha_prime_plus_indices, vtype = GRB.BINARY, name = 'alpha_prime_plus')
    beta_prime = MODEL.addVars(beta_prime_indices, vtype = GRB.BINARY, name = 'beta_prime')

    gamma = MODEL.addVars(gamma_indices, vtype = GRB.BINARY, name = 'gamma')

    dif = MODEL.addVars(dist_indices, 2, vtype = GRB.CONTINUOUS, name = 'dif')
    
    x = MODEL.addVars(3, 2, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 10.0, name = 'x')

    bigM = 10

    for i in lambda_indices:
        MODEL.addConstrs(A[i, dim] == x[0, dim] + landa[i]*(x[1, dim] - x[0, dim]) + landa_prime[i]*(x[2, dim] - x[0, dim]) for dim in range(2))
        # MODEL.addConstrs(A[i, dim] == x[0, dim] + w_landa[i, dim] + w_landa_prime[i, dim] for dim in range(2))

        MODEL.addConstr(landa[i] <= bigM*(alpha_minus[i] + alpha_plus[i]))
        MODEL.addConstr(landa[i] >= -bigM*(1 - alpha_minus[i] + alpha_plus[i]))
        MODEL.addConstr(landa[i] - 1 <= bigM*(1 - alpha_minus[i] + alpha_plus[i]))
        MODEL.addConstr(landa[i] - 1 >= -bigM*(2 - alpha_minus[i] - alpha_plus[i]))

        MODEL.addConstr(alpha_minus[i] >= alpha_plus[i])

        MODEL.addConstr(beta[i] <= alpha_minus[i])
        MODEL.addConstr(beta[i] <= 1 - alpha_plus[i])
        MODEL.addConstr(beta[i] >= alpha_minus[i] - alpha_plus[i])



        MODEL.addConstr(landa_prime[i] <= bigM*(alpha_prime_minus[i] + alpha_prime_plus[i]))
        MODEL.addConstr(landa_prime[i] >= -bigM*(1 - alpha_prime_minus[i] + alpha_prime_plus[i]))
        MODEL.addConstr(landa_prime[i] - 1 + landa[i] <= bigM*(1 - alpha_prime_minus[i] + alpha_prime_plus[i]))
        MODEL.addConstr(landa_prime[i] - 1 + landa[i] >= -bigM*(2 - alpha_prime_minus[i] - alpha_prime_plus[i]))

        MODEL.addConstr(alpha_prime_minus[i] >= alpha_prime_plus[i])

        MODEL.addConstr(beta_prime[i] <= alpha_prime_minus[i])
        MODEL.addConstr(beta_prime[i] <= 1 - alpha_prime_plus[i])
        MODEL.addConstr(beta_prime[i] >= alpha_prime_minus[i] - alpha_prime_plus[i])

        MODEL.addConstr(gamma[i] >= beta[i] + beta_prime[i] - 1)
        MODEL.addConstr(2*gamma[i] <= beta[i] + beta_prime[i] )
        
    
    for i in dist_indices:
        if i < 2:
            MODEL.addConstrs(dif[i, dim] >=  x[i, dim] - x[i+1, dim] for dim in range(2))
            MODEL.addConstrs(dif[i, dim] >= -x[i, dim] + x[i+1, dim] for dim in range(2))
        else:
            MODEL.addConstrs(dif[i, dim] >=  x[i, dim] - x[0, dim] for dim in range(2))
            MODEL.addConstrs(dif[i, dim] >= -x[i, dim] + x[0, dim] for dim in range(2))
        
        
    
    MODEL.addConstr(gp.quicksum(dif[i, 0]*dif[i, 0] + dif[i, 1]*dif[i, 1] for i in dist_indices) <= L*L)

        

    objective = gp.quicksum(gamma[i] for i in gamma.keys())


    MODEL.setObjective(objective, GRB.MAXIMIZE)

    MODEL.Params.NonConvex = 2

    MODEL.optimize()

    print(landa)
    print(alpha_minus)
    print(alpha_plus)
    print(x)




    
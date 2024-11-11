import numpy as np
from gurobipy import *

""" The input parameter """
facility_num = 3 # m, indexed by i
customer_num = 3 # n, indexed by j
fixed_cost = [400, 414, 326] # f_i: fixed cost at location i if built
unit_capacity_cost = [18, 25, 20] # a_i
trans_cost = [[22, 33, 24],
              [33, 23, 30],
              [20, 25, 27]] # c_ij: transportation cost matrix, m by n matrix
max_capacity = 800 # K_i: maximal allowable capacity. In this case it is set identical for all locations

demand_nominal = [206, 274, 220] # d_j_low: basic demand for each j
demand_var = [40, 40, 40] # d_j_tilde: maximal demand deviation for each j

""" Create variables  """
model = Model('Location-transportation model')
x = {}
z = {}
y = {}

for i in range(facility_num):
    y[i] = model.addVar(vtype = GRB.BINARY, name = 'y_%s'%(i))
    z[i] = model.addVar(lb = 0, vtype = GRB.CONTINUOUS, name = 'y_%s'%(i))
    for j in range(customer_num):
        x[i, j] = model.addVar(lb = 0, vtype = GRB.CONTINUOUS, name = 'x_%s_%s'%(i, j))

""" Set objective """
obj = LinExpr()
for i in range(facility_num):
    obj.addTerms(fixed_cost[i], y[i])
    obj.addTerms(unit_capacity_cost[i], z[i])
    for j in range(customer_num):
        obj.addTerms(trans_cost[i][j], x[i, j])
model.setObjective(obj, GRB.MINIMIZE)

""" Add Constraints  """
# cons 1
for i in range(facility_num):
    model.addConstr(z[i] <= max_capacity * y[i])

# cons 2
for i in range(facility_num):
    expr = LinExpr()
    for j in range(customer_num):
        expr.addTerms(1, x[i, j])
    model.addConstr(expr <= z[i])

# cons 3
for j in range(facility_num):
    expr = LinExpr()
    for i in range(customer_num):
        expr.addTerms(1, x[i, j])
    model.addConstr(expr >= demand_nominal[j])

model.optimize()

print(f'Obj = {model.ObjVal}')
print('-----  location ----- ')
for key in z.keys():
    print(f'facility : {key}, location: {y[key].x}, capacity: {z[key].x}')
print(x)
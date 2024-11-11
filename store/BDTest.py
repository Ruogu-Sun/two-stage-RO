from gurobipy import *
import numpy as np

"""Parameter Setting"""
facility_num = 3
customer_num = 3

bigM = 10**3
LB = -GRB.INFINITY
UB = GRB.INFINITY
epsilon = 1e-5
iter_cnt = 1
max_iter=10
# benders_cuts = []  # List to store cuts


"""Cost Setting"""
fixed_cost = [400, 414, 326]
unit_capacity_cost = [18, 25, 20]
trans_cost = [[22, 33, 24],
              [33, 23, 30],
              [20, 25, 27]]
Gamma=0.6
Kapacity=[800,800,800]
# max_capacity = 800
demand_nominal = [206, 274, 220]
demand_var = [40, 40, 40]


"""Master Problem"""
MP = Model('Benders decomposition-MP')
# x_master = {}
z = {}
y = {}

# Set Variables
for i in range(facility_num):
    y[i] = MP.addVar(vtype=GRB.BINARY, name=f'y_{i}')
    z[i] = MP.addVar(lb=0, vtype=GRB.CONTINUOUS, name=f'z_{i}')
eta = MP.addVar(lb=0, vtype=GRB.CONTINUOUS, name='eta')

""" Set objective """
obj = LinExpr()
for i in range(facility_num):
    obj.addTerms(fixed_cost[i], y[i])
    obj.addTerms(unit_capacity_cost[i], z[i])
obj.addTerms(1, eta)

MP.setObjective(obj, GRB.MINIMIZE)


# Set Constraints
# cons1, each z_i below maximal capacity.
for i in range(facility_num):
    MP.addConstr(z[i] <= Kapacity[i] * y[i])

# Add Benders cuts from previous iterations
# for cut in benders_cuts:
#     MP.addConstr(cut)


MP.optimize()
# Get optimal value from MP
MP_obj = MP.ObjVal

# Get solution from MP
# y_sol = MP.getAttr('x', y)

z_sol = {}
for key in z.keys():
    z_sol[key] = z[key].X

LB = max(MP_obj, LB)


"""Demand Worst Problem"""
DP=Model('Benders decomposition-DP')
g={}
d={}

for j in range(customer_num):
    g[j] = DP.addVar(lb=0,ub=1,vtype=GRB.CONTINUOUS,name=f"g_{j}")
    d[j] = DP.addVar(lb=demand_nominal[j],ub=demand_nominal[j]+demand_var[j],vtype=GRB.CONTINUOUS,name=f'd_{j}')

# OBJ
obj=LinExpr()
for j in range(customer_num):
    obj.addTerms(1,d[j])
DP.setObjective(obj,GRB.MAXIMIZE)

# Constraint
expr=LinExpr()
for j in range(customer_num):
    expr.addTerms(1,g[j])
DP.addConstr(expr<=Gamma*customer_num)

for j in range(customer_num):
    expr=LinExpr()
    expr.addTerms(1,d[j])
    expr.addTerms(-demand_var[j],g[j])
    DP.addConstr(expr==demand_nominal[j],name=f"Demand_Uncertainty_{j}")

DP.optimize()

d_sol = {}
for key in d.keys():
    d_sol[key] = d[key].X
print(d_sol)



SP_dual = Model('SP_dual')
pi={}
lam={}
for i in range(facility_num):
    pi[i] = SP_dual.addVar(lb=0,vtype=GRB.CONTINUOUS, name=f'pi_{i}')

for j in range(customer_num):
    lam[j] = SP_dual.addVar(lb=0,vtype=GRB.CONTINUOUS,name=f"lam_{j}")
    # g[j] = SP_dual.addVar(lb=0,ub=1,vtype=GRB.CONTINUOUS,name=f"g_{j}")
    # d[j] = SP_dual.addVar(lb=demand_nominal[j],ub=demand_nominal[j]+demand_var[j],vtype=GRB.CONTINUOUS,name=f'd_{j}')
    # d[j] = SP_dual.addVar(lb=0, vtype=GRB.CONTINUOUS,
    #                       name=f'd_{j}')
    # d[j]=demand_nominal[j]+g[j]*demand_var[j]

""" Set objective """
obj = LinExpr()
for i in range(facility_num):obj.addTerms(-z_sol[i], pi[i])
for j in range(customer_num):obj.addTerms(d_sol[j],lam[j])

SP_dual.setObjective(obj, GRB.MAXIMIZE)

"""Set Constraints"""
# Duality Constraint
def getTotalTransCost(matrix):
    total_sum = 0
    for item in matrix:
        if isinstance(item, list):total_sum += sum(item)
        else:total_sum += item
    return total_sum
#
expr=LinExpr()
for i in range(facility_num):expr.addTerms(1,pi[i])
for j in range(customer_num):expr.addTerms(-1,lam[j])
# for i in range(facility_num):
#     sum(trans_cost[i])
#
SP_dual.addConstr(expr>=-getTotalTransCost(trans_cost),name="SP_Duality_Feasibility")

# Dual constraints
# for i in range(facility_num):
#     for j in range(customer_num):
#         SP_dual.addConstr(trans_cost[i][j] + pi[i] - lam[j] >= 0, name=f"DualConstraint_{i}_{j}")

SP_dual.setParam('InfUnbdInfo', 1)
SP_dual.optimize()
#

#

print(SP_dual.status)
SP_dual_obj = SP_dual.ObjVal
UB = min(UB, LB-eta.x+SP_dual_obj)

print(f"MP_obj:{MP_obj} \n LB:{LB} \n eta:{eta.x} \n SP_dual_obj: {SP_dual_obj} \n UB:{UB}")
#
#
while (UB - LB > epsilon and iter_cnt <= max_iter):
    iter_cnt+=1
    if(SP_dual.status == 2):
        print(f"At Iteration {iter_cnt}, an OPTIMALITY CUT is to be added")
        """Adding Optimality Cut & Reevaluate eta"""
        MP.reset()
        expr=LinExpr()
        cons=0
        for i in range(facility_num):expr.addTerms(-pi[i].x,z[i])
        for j in range(customer_num):cons+=lam[j].x*d_sol[j]
        MP.addConstr(eta>=cons+expr,name=f"optimality_cut_{iter_cnt}")
        MP.optimize()
        MP_obj=MP.ObjVal
        LB=max(LB,MP_obj)

        z_sol = {}
        for key in z.keys():
            z_sol[key] = z[key].x

        """SP_dual Reset"""
        SP_dual.reset()

        obj = LinExpr()
        for i in range(facility_num): obj.addTerms(-z_sol[i], pi[i])
        for j in range(customer_num):obj.addTerms(d_sol[j], lam[j])

        SP_dual.setObjective(obj, GRB.MAXIMIZE)

        # SP_dual.setParam('InfUnbdInfo', 1)
        SP_dual.optimize()
        SP_dual_obj = SP_dual.ObjVal
        UB = min(UB, LB - eta.x + SP_dual_obj)

    elif (SP_dual.status == 5):
        print(f"At Iteration {iter_cnt}, a FEASIBILITY CUT is to be added")
        """Adding Feasibility Cut"""
        MP.reset()
        expr = LinExpr()
        cons = 0
        for i in range(facility_num): expr.addTerms(-pi[i].unbdRay, z[i])
        for j in range(customer_num): cons += lam[j].unbdRay * d_sol[j]
        MP.addConstr(0 >= cons + expr, name=f"feasibility_cut_{iter_cnt}")
        MP.optimize()
        MP_obj = MP.ObjVal
        LB = max(LB, MP_obj)

        z_sol = {}
        for key in z.keys():
            z_sol[key] = z[key].x
        print(f"LB:{LB}" )

        """SP_dual Reset"""
        SP_dual.reset()

        obj = LinExpr()
        for i in range(facility_num): obj.addTerms(-z_sol[i], pi[i])
        for j in range(customer_num):obj.addTerms(d_sol[j], lam[j])

        SP_dual.setObjective(obj, GRB.MAXIMIZE)

        SP_dual.setParam('InfUnbdInfo', 1)
        SP_dual.optimize()
        SP_dual_obj = SP_dual.ObjVal
        UB = min(UB, LB - eta.x + SP_dual_obj)
        print(f"UB:{UB}")


print(f"LB{LB}")
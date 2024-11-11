from gurobipy import *
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

""" The input parameter """
facility_num = 25
customer_num = 25

def fixed_cost(m):return np.random.randint(100, 500, size=m)
# fixed_cost = [400, 414, 326]
fixed_cost=fixed_cost(m=facility_num)

def unit_capacity_cost(m):return np.random.randint(10, 60, size=m)
# unit_capacity_cost = [18, 25, 20]
unit_capacity_cost=unit_capacity_cost(m=facility_num)

def trans_cost(m,n):return np.random.randint(10, 40, size=m*n).reshape((m,n))
trans_cost = trans_cost(m=facility_num,n=customer_num)

# trans_cost = [[22, 33, 24],
#               [33, 23, 30],
#               [20, 25, 27]]

Gamma=0.6
# max_capacity = 800
Kapacity=np.random.randint(100000, 200000, size=facility_num)
demand_nominal=np.random.randint(50, 400, size=customer_num)
alpha=np.random.choice([0.1,0.2,0.3,0.4,0.5], size=customer_num)
demand_var = alpha*demand_nominal
feasibility=sum(demand_nominal+Gamma*demand_var)


""" build initial master problem """
""" Create variables """
master = Model('master problem')
x_master = {}
z = {}
y = {}

for i in range(facility_num):
    y[i] = master.addVar(vtype=GRB.BINARY, name=f'y_{i}')
    z[i] = master.addVar(lb=0, vtype=GRB.CONTINUOUS, name=f'z_{i}')
eta = master.addVar(lb=0, vtype=GRB.CONTINUOUS, name='eta')

""" Set objective """
obj = LinExpr()
for i in range(facility_num):
    obj.addTerms(fixed_cost[i], y[i])
    obj.addTerms(unit_capacity_cost[i], z[i])
obj.addTerms(1, eta)

master.setObjective(obj, GRB.MINIMIZE)

""" Add Constraints  """
# cons 1
for i in range(facility_num):
    master.addConstr(z[i] <= Kapacity[i] * y[i])

""" Add initial value Constraints  """
# create new variables x
iter_cnt = 0
for i in range(facility_num):
    for j in range(customer_num):
        x_master[iter_cnt, i, j] = master.addVar(lb=0
                                                 , ub=GRB.INFINITY
                                                 , vtype=GRB.CONTINUOUS
                                                 , name=f'x_{iter_cnt}_{i}_{j}')

# create new constraints: cutting plane
expr = LinExpr()
for i in range(facility_num):
    for j in range(customer_num):
        expr.addTerms(trans_cost[i][j], x_master[iter_cnt, i, j])
master.addConstr(eta >= expr)

"""Add feasibility constraint"""
expr = LinExpr()
for i in range(facility_num):
    expr.addTerms(1, z[i])
master.addConstr(expr >= feasibility)  # 206 + 274 + 220 + 40 * 1.8

""" solve the model and output """
master.optimize()
print(f'Obj = {master.ObjVal}')
print('-----  location ----- ')
for key in z.keys():
    print('facility : {}, location: {}, capacity: {}'.format(key, y[key].x, z[key].x))


""" Column-and-constraint generation """

LB = -np.inf
UB = np.inf
iter_cnt = 0
max_iter = 30
cut_pool = {}
eps = 0.001
Gap = np.inf

"""z^* for this step's MP"""
z_sol = {}
for key in z.keys():
    z_sol[key] = z[key].x
# print(z_sol)

""" solve the master problem and update bound """
# master.optimize()

""" 
 Update the Lower bound 
"""
LB = master.ObjVal
print('LB: {}'.format(LB))



''' create the subproblem '''
subProblem = Model('sub problem')
# transportation decision variables in subproblem
x = {}
# true demand
d = {}
# uncertainty part: var part
g = {}
# dual variable
pi = {}
# dual variable
theta = {}
# aux var
v = {}
# aux var
w = {}
# aux var
h = {}

big_M = 1000

for i in range(facility_num):
    pi[i] = subProblem.addVar(lb=-GRB.INFINITY, ub=0, vtype=GRB.CONTINUOUS, name=f'pi_{i}')
    v[i] = subProblem.addVar(vtype=GRB.BINARY, name=f'v_{i}')

for j in range(customer_num):
    w[j] = subProblem.addVar(vtype=GRB.BINARY, name=f'w_{j}')
    g[j] = subProblem.addVar(lb=0, ub=1, vtype=GRB.CONTINUOUS, name=f'g_{j}')
    theta[j] = subProblem.addVar(lb=-GRB.INFINITY, ub=0, vtype=GRB.CONTINUOUS, name=f'theta_{j}')
    d[j] = subProblem.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name=f'd_{j}')
for i in range(facility_num):
    for j in range(customer_num):
        h[i, j] = subProblem.addVar(vtype=GRB.BINARY, name=f'h_{i}_{j}')
        x[i, j] = subProblem.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name=f'x_{i}_{j}')




""" set objective """
sub_obj = LinExpr()
for i in range(facility_num):
    for j in range(customer_num):
        sub_obj.addTerms(trans_cost[i][j], x[i, j])
subProblem.setObjective(sub_obj, GRB.MAXIMIZE)

""" add constraints to subproblem """
# cons 1: not exceed inventory
for i in range(facility_num):
    expr = LinExpr()
    for j in range(customer_num):
        expr.addTerms(1, x[i, j])
    subProblem.addConstr(expr <= z_sol[i], name=f'sub_capacity_1_z_{i}')

# cons 2: demand be satisfied!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
for j in range(customer_num):
    expr = LinExpr()
    for i in range(facility_num):
        expr.addTerms(1, x[i, j])
    subProblem.addConstr(expr >= d[j])

# cons 3: duality feasibility constraint
for i in range(facility_num):
    for j in range(customer_num):
        subProblem.addConstr(pi[i] - theta[j] <= trans_cost[i][j])

""" demand constraints """
for j in range(customer_num):
    subProblem.addConstr(d[j] == demand_nominal[j] + g[j] * demand_var[j])

"""Gamma-related constraint"""
# subProblem.addConstr(g[0] + g[1] + g[2] <= 1.8)
subProblem.addConstr(sum(g.values()) <= Gamma*customer_num)
# subProblem.addConstr(g[0] + g[1] <= 1.2)


""" logic constraints """
# logic 1 complementary slackness condition for Lagragian Multipliers \pi
for i in range(facility_num):
    subProblem.addConstr(-pi[i] <= big_M * v[i])
    expr = LinExpr()
    for j in range(customer_num):
        expr.addTerms(1, x[i, j])
    subProblem.addConstr(z_sol[i] - expr <= big_M - big_M * v[i], name=f'sub_capacity_2_z_{i}')

# logic 2 complementary slackness condition for Lagragian Multipliers \theta
for j in range(customer_num):
    subProblem.addConstr(-theta[j] <= big_M * w[j])
    expr = LinExpr()
    for i in range(facility_num):
        expr.addTerms(1, x[i, j])
    subProblem.addConstr(expr - d[j] <= big_M - big_M * w[j])


# logic 3 First order derivative Stationarity condition \theta
for j in range(customer_num):
    for i in range(facility_num):
        subProblem.addConstr(x[i, j] <= big_M * h[i, j])
        subProblem.addConstr(trans_cost[i][j] - pi[i] + theta[j] <= big_M - big_M * h[i, j])


# subProblem.write('SP.lp')
subProblem.optimize()




d_sol = {}

print('\n\n\n *******C&CG starts ******* ')
print('\n **  Initial Solution  ** ')

def print_sub_sol(model, d, g, x):
    d_sol = {}
    if(model.status != 2):
        print('The problem is infeasible or unbounded!')
        print(f'Status: {model.status}')
        for j in range(customer_num):
            d_sol[j] = demand_nominal[j]
            # +g[j].x*demand_var[j]
    else:
        print(f'Obj(sub) : {model.ObjVal}', end='\t | ')
        for key in d.keys():
            # print('demand: {}, perturbation = {}'.format(d[key].x, g[key].x))
            d_sol[key] = d[key].x
    return d_sol


d_sol = print_sub_sol(subProblem, d, g, x)

UB = min(UB, subProblem.ObjVal + master.ObjVal - eta.x)

print(f'UB (iter {iter_cnt}): {UB}')

# close the outputflag

master.setParam('Outputflag', 0)

subProblem.setParam('Outputflag', 0)


lb=[]
ub=[]
kk=[]
lb.append(LB)
ub.append(UB)
kk.append(iter_cnt)

# master.reset() # discard solution information
"""Main loop of CCG algorithm"""
while (UB - LB )/UB> eps and iter_cnt <= max_iter:
    iter_cnt += 1
    # print('\n\n --- iter : {} --- \n'.format(iter_cnt))
    print(f'\n iter : {iter_cnt} ', end='\t | ')
    # create new variables x
    for i in range(facility_num):
        for j in range(customer_num):
            x_master[iter_cnt, i, j] = master.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name=f'x_{iter_cnt}_{i}_{j}')

    # if subproblem is feasible and bound, create variables xk+1 and add the new constraints
    if (subProblem.status == 2):
        master.reset()
        # create new constraints
        expr = LinExpr()
        for i in range(facility_num):
            for j in range(customer_num):
                expr.addTerms(trans_cost[i][j], x_master[iter_cnt, i, j])
        master.addConstr(eta >= expr)

        # create worst case related constraints
        # cons 2 not exceed inventory
        for i in range(facility_num):
            expr = LinExpr()
            for j in range(customer_num):
                expr.addTerms(1, x_master[iter_cnt, i, j])
            master.addConstr(expr <= z[i])

        # cons 3 satisfy demand d^*
        for j in range(customer_num):
            expr = LinExpr()
            for i in range(facility_num):
                expr.addTerms(1, x_master[iter_cnt, i, j])
            master.addConstr(expr >= d_sol[j])

        # solve the resulted master problem
        master.optimize()
        print(f'Obj(master): {master.ObjVal}', end='\t | ')
        """ Update the LB """
        LB = master.ObjVal
        print(f'LB (iter {iter_cnt}): {LB}', end='\t | ')

        """ Update the subproblem """
        subProblem.reset()
        # first, get z_sol from updated master problem
        for key in z.keys():
            z_sol[key] = z[key].x

        for i in range(facility_num):
            constr_name_1 = f'sub_capacity_1_z_{i}'
            constr_name_2 = f'sub_capacity_2_z_{i}'
            subProblem.remove(subProblem.getConstrByName(constr_name_1))
            subProblem.remove(subProblem.getConstrByName(constr_name_2))
        subProblem.update()

        # add new constraints
        # cons 1
        for i in range(facility_num):
            expr = LinExpr()
            for j in range(customer_num):
                expr.addTerms(1, x[i, j])
            subProblem.addConstr(expr <= z_sol[i], name=f'sub_capacity_1_z_{i}')

        # logic 1
        for i in range(facility_num):
            subProblem.addConstr(-pi[i] <= big_M * v[i])
            expr = LinExpr()
            for j in range(customer_num):
                expr.addTerms(1, x[i, j])
            subProblem.addConstr(z_sol[i] - expr <= big_M - big_M * v[i], name=f'sub_capacity_2_z_{i}')

        """ Update the lower bound """
        subProblem.optimize()
        # subProblem.reset() # new SP because new eta(cutting plane) will be added in
        d_sol = print_sub_sol(subProblem,d,g,x)

        """Update the Upper bound"""
        if (subProblem.status == 2):
            UB = min(UB,subProblem.ObjVal + master.ObjVal - eta.x)
        # print('eta = {}'.format(eta.x))
        print(f'UB (iter {iter_cnt}): {UB}', end='\t | ')
        Gap = round(100 * (UB - LB) / UB, 2)
        print(f'eta = {eta.x}', end='\t | ')
        print(f' Gap: {Gap} % ', end='\t')

    # If the subproblem is unbounded
    if (subProblem.status == 4):
        master.reset()
        # create worst case related constraints
        # cons 2
        for i in range(facility_num):
            expr = LinExpr()
            for j in range(customer_num):
                expr.addTerms(1, x_master[iter_cnt, i, j])
            master.addConstr(expr <= z[i])

        # cons 3
        for j in range(customer_num):
            expr = LinExpr()
            for i in range(facility_num):
                expr.addTerms(1, x_master[iter_cnt, i, j])
            master.addConstr(expr >= d_sol[j])

        # solve the resulted master problem
        master.optimize()
        print(f"Obj(master): {master.ObjVal}")

        """ Update the LB """
        LB = master.ObjVal
        print(f'LB (iter {iter_cnt}): {LB}')

        """ Update the subproblem """
        # first, get z_sol from updated master problem
        for key in z.keys():
            z_sol[key] = z[key].x

        subProblem.reset()
        # change the coefficient of subproblem
        for i in range(facility_num):
            constr_name_1 = f'sub_capacity_1_z_{i}'
            constr_name_2 = f'sub_capacity_2_z_{i}'
            subProblem.remove(subProblem.getConstrByName(constr_name_1))
            subProblem.remove(subProblem.getConstrByName(constr_name_2))
        subProblem.update()

        # add new constraints
        # cons 1
        for i in range(facility_num):
            expr = LinExpr()
            for j in range(customer_num):
                expr.addTerms(1, x[i, j])
            subProblem.addConstr(expr <= z_sol[i], name=f'sub_capacity_1_z_{i}')

        # logic 1
        for i in range(facility_num):
            # subProblem.addConstr(-pi[i] <= big_M * v[i])
            expr = LinExpr()
            for j in range(customer_num):
                expr.addTerms(1, x[i, j])
            subProblem.addConstr(z_sol[i] - expr <= big_M - big_M * v[i], name=f'sub_capacity_2_z_{i}')

        """ Update the lower bound """
        subProblem.optimize()
        d_sol = print_sub_sol(subProblem, d, g, x)

        """Update the Upper bound"""
        if (subProblem.status == 2):
            UB = min(UB,subProblem.ObjVal + master.ObjVal - eta.x)
            print(f'eta = {eta.x}')
        print(f'UB (iter {iter_cnt}): {UB}')
        Gap = round(100 * (UB - LB) / UB, 2)
        print(f'---- Gap: {Gap} % ---- ')
    for j in range(customer_num):
        print(f'\t actual demand: {d[j].varName}: {d[j].x}', end='')
        print(f'\t perturbation in worst case: {g[j].varName}: {g[j].x}')
    lb.append(LB)
    ub.append(UB)
    kk.append(iter_cnt)

# master.write('finalMP.lp')
print('\n\nOptimal solution found !')
print(f'Opt_Obj : {master.ObjVal}')

print(f' ** Final Gap: {Gap} % ** ')
print('\n** Solution ** ')
for i in range(facility_num):
    print(f' {y[i].varName}: {y[i].x},\t{z[i].varName}: {z[i].x} ', end='')
for j in range(customer_num):
    print(f'\t actual demand: {d[j].varName}: {d[j].x}', end='')
    print(f'\t perturbation in worst case: {g[j].varName}: {g[j].x}')
print('\n** Transportation solution ** ')
for i in range(facility_num):
    for j in range(customer_num):
        if (x[i, j].x > 0):
            print(f'trans: {x[i, j].varName}: {x[i, j].x}, cost : {trans_cost[i][j] * x[i, j].x} \t ',end='')
    print()

plt.figure(figsize=(10, 6))

sns.lineplot(x=kk, y=lb, label='Lower Bound (lb)', color='blue', marker='o')
sns.lineplot(x=kk, y=ub, label='Upper Bound (ub)', color='orange', marker='o')

plt.title('Line Plot of lb and ub')
plt.xlabel('k')
plt.ylabel('Values')
plt.legend()

plt.grid()
plt.show()
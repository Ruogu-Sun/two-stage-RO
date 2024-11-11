from gurobipy import *
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

""" The input parameter """
facility_num = 3
customer_num = 3
fixed_cost = [400, 414, 326]
unit_capacity_cost = [18, 25, 20]
trans_cost = [[22, 33, 24],
              [33, 23, 30],
              [20, 25, 27]]
max_capacity = 800

demand_nominal = [206, 274, 220]
demand_var = [40, 40, 40]

""" build initial master problem """
""" Create variables """
master = Model('master problem')
x_master = {}
z = {}
y = {}
eta = master.addVar(lb=0, vtype=GRB.CONTINUOUS, name='eta')

for i in range(facility_num):
    y[i] = master.addVar(vtype=GRB.BINARY, name='y_%s' % (i))
    z[i] = master.addVar(lb=0, vtype=GRB.CONTINUOUS, name='z_%s' % (i))

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
    master.addConstr(z[i] <= max_capacity * y[i])

""" Add initial value Constraints  """
# create new variables x
iter_cnt = 0
for i in range(facility_num):
    for j in range(customer_num):
        x_master[iter_cnt, i, j] = master.addVar(lb=0
                                                 , ub=GRB.INFINITY
                                                 , vtype=GRB.CONTINUOUS
                                                 , name='x_%s_%s_%s' % (iter_cnt, i, j))


# create new constraints
expr = LinExpr()
for i in range(facility_num):
    for j in range(customer_num):
        expr.addTerms(trans_cost[i][j], x_master[iter_cnt, i, j])
master.addConstr(eta >= expr)

expr = LinExpr()
for i in range(facility_num):
    expr.addTerms(1, z[i])
master.addConstr(expr >= 772)  # 206 + 274 + 220 + 40 * 1.8

""" solve the model and output """
master.optimize()
print('Obj = {}'.format(master.ObjVal))
print('-----  location ----- ')
for key in z.keys():
    print('facility : {}, location: {}, capacity: {}'.format(key, y[key].x, z[key].x))


def print_sub_sol(model, d, g, x):
    d_sol = {}
    if(model.status != 2):
        print('The problem is infeasible or unbounded!')
        print('Status: {}'.format(model.status))
        for j in range(customer_num):
            d_sol[j] = 0
    else:
        print(f'Obj(sub) : {model.ObjVal}', end='\t | ')
        for key in d.keys():
            # print('demand: {}, perturbation = {}'.format(d[key].x, g[key].x))
            d_sol[key] = d[key].x
    return d_sol


""" Column-and-constraint generation """

LB = -np.inf
UB = np.inf
iter_cnt = 0
max_iter = 30
cut_pool = {}
eps = 0.001
Gap = np.inf

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

big_M = 100000

for i in range(facility_num):
    pi[i] = subProblem.addVar(lb=-GRB.INFINITY, ub=0, vtype=GRB.CONTINUOUS, name='pi_%s' % i)

    v[i] = subProblem.addVar(vtype=GRB.BINARY, name='v_%s' % i)
for j in range(customer_num):
    w[j] = subProblem.addVar(vtype=GRB.BINARY, name='w_%s' % j)
    g[j] = subProblem.addVar(lb=0, ub=1, vtype=GRB.CONTINUOUS, name='g_%s' % j)
    theta[j] = subProblem.addVar(lb=-GRB.INFINITY, ub=0, vtype=GRB.CONTINUOUS, name='theta_%s' % j)
    d[j] = subProblem.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name='d_%s' % j)
for i in range(facility_num):
    for j in range(customer_num):
        h[i, j] = subProblem.addVar(vtype=GRB.BINARY, name='h_%s_%s' % (i, j))
        x[i, j] = subProblem.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name='x_%s_%s' %(i, j))




""" set objective """
sub_obj = LinExpr()
for i in range(facility_num):
    for j in range(customer_num):
        sub_obj.addTerms(trans_cost[i][j], x[i, j])
subProblem.setObjective(sub_obj, GRB.MAXIMIZE)

""" add constraints to subproblem """
# cons 1
for i in range(facility_num):
    expr = LinExpr()
    for j in range(customer_num):
        expr.addTerms(1, x[i, j])
    subProblem.addConstr(expr <= z_sol[i], name='sub_capacity_1_z_%s' % i)

# cons 2
for j in range(facility_num):
    expr = LinExpr()
    for i in range(customer_num):
        expr.addTerms(1, x[i, j])
    subProblem.addConstr(expr >= d[j])

# cons 3
for i in range(facility_num):
    for j in range(customer_num):
        subProblem.addConstr(pi[i] - theta[j] <= trans_cost[i][j])

""" demand constraints """
for j in range(customer_num):
    subProblem.addConstr(d[j] == demand_nominal[j] + g[j] * demand_var[j])

subProblem.addConstr(g[0] + g[1] + g[2] <= 1.8)
subProblem.addConstr(g[0] + g[1] <= 1.2)


""" logic constraints """
# logic 1
for i in range(facility_num):
    subProblem.addConstr(-pi[i] <= big_M * v[i])
    expr = LinExpr()
    for j in range(customer_num):
        expr.addTerms(1, x[i, j])
    subProblem.addConstr(z_sol[i] - expr <= big_M - big_M * v[i], name='sub_capacity_2_z_%s' % i)

# logic 2
for j in range(customer_num):
    subProblem.addConstr(-theta[j] <= big_M * w[j])
    expr = LinExpr()
    for i in range(facility_num):
        expr.addTerms(1, x[i, j])
    subProblem.addConstr(expr - d[j] <= big_M - big_M * w[j])


# logic 3
for j in range(customer_num):
    for i in range(facility_num):
        subProblem.addConstr(x[i, j] <= big_M * h[i, j])
        subProblem.addConstr(trans_cost[i][j] - pi[i] + theta[j] <= big_M - big_M * h[i, j])
# subProblem.write('SP.lp')
subProblem.optimize()
d_sol = {}

print('\n\n\n *******C&CG starts ******* ')
print('\n **  Initial Solution  ** ')

d_sol = print_sub_sol(subProblem, d, g, x)

UB = min(UB, subProblem.ObjVal + master.ObjVal - eta.x)

print('UB (iter {}): {}'.format(iter_cnt, UB))

# close the outputflag

master.setParam('Outputflag', 0)

subProblem.setParam('Outputflag', 0)


lb=[]
ub=[]
kk=[]
lb.append(LB)
ub.append(UB)
kk.append(iter_cnt)


"""Main loop of CCG algorithm"""
while (UB - LB > eps and iter_cnt <= max_iter):
    iter_cnt += 1
    # print('\n\n --- iter : {} --- \n'.format(iter_cnt))
    print('\n iter : {} '.format(iter_cnt), end='\t | ')
    # create new variables x
    for i in range(facility_num):
        for j in range(customer_num):
            x_master[iter_cnt, i, j] = master.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name='x_%s_%s_%s' % (iter_cnt, i, j))

    # if subproblem is feasible and bound, create variables xk+1 and add the new constraints
    if (subProblem.status == 2):
        # create new constraints
        expr = LinExpr()
        for i in range(facility_num):
            for j in range(customer_num):
                expr.addTerms(trans_cost[i][j], x_master[iter_cnt, i, j])
        master.addConstr(eta >= expr)

        # create worst case related constraints
        # cons 2
        for i in range(facility_num):
            expr = LinExpr()
            for j in range(customer_num):
                expr.addTerms(1, x_master[iter_cnt, i, j])
            master.addConstr(expr <= z[i])

        # cons 3
        for j in range(facility_num):
            expr = LinExpr()
            for i in range(customer_num):
                expr.addTerms(1, x_master[iter_cnt, i, j])
            master.addConstr(expr >= d_sol[j])

        # solve the resulted master problem
        master.optimize()
        print('Obj(master): {}'.format(master.ObjVal), end='\t | ')
        """ Update the LB """
        LB = master.ObjVal
        print('LB (iter {}): {}'.format(iter_cnt, LB), end='\t | ')

        """ Update the subproblem """
        # first, get z_sol from updated master problem
        for key in z.keys():
            z_sol[key] = z[key].x

        for i in range(facility_num):
            constr_name_1 = 'sub_capacity_1_z_' + str(i)
            constr_name_2 = 'sub_capacity_2_z_' + str(i)
            subProblem.remove(subProblem.getConstrByName(constr_name_1))
            subProblem.remove(subProblem.getConstrByName(constr_name_2))

        # add new constraints
        # cons 1
        for i in range(facility_num):
            expr = LinExpr()
            for j in range(customer_num):
                expr.addTerms(1, x[i, j])
            subProblem.addConstr(expr <= z_sol[i], name='sub_capacity_1_z_%s' % i)

        # logic 1
        for i in range(facility_num):
            # subProblem.addConstr(pi[i] <= big_M * v[i])
            expr = LinExpr()
            for j in range(customer_num):
                expr.addTerms(1, x[i, j])
            subProblem.addConstr(z_sol[i] - expr <= big_M - big_M * v[i], name='sub_capacity_2_z_%s' % i)

        """ Update the lower bound """
        subProblem.optimize()
        d_sol = print_sub_sol(subProblem,d,g,x)

        """Update the Upper bound"""
        if (subProblem.status == 2):
            UB = min(UB, subProblem.ObjVal + master.ObjVal - eta.x)
        # print('eta = {}'.format(eta.x))
        print('UB (iter {}): {}'.format(iter_cnt, UB), end='\t | ')
        Gap = round(100 * (UB - LB) / UB, 2)
        print('eta = {}'.format(eta.x), end='\t | ')
        print(' Gap: {} % '.format(Gap), end='\t')

    # If the subproblem is unbounded
    if (subProblem.status == 4):
        # create worst case related constraints
        # cons 2
        for i in range(facility_num):
            expr = LinExpr()
            for j in range(customer_num):
                expr.addTerms(1, x_master[iter_cnt, i, j])
            master.addConstr(expr <= z[i])

        # cons 3
        for j in range(facility_num):
            expr = LinExpr()
            for i in range(customer_num):
                expr.addTerms(1, x_master[iter_cnt, i, j])
            master.addConstr(expr >= d_sol[j])

        # solve the resulted master problem
        master.optimize()
        print(f"Obj(master): {master.ObjVal}")

        """ Update the LB """
        LB = master.ObjVal
        print('LB (iter {}): {}'.format(iter_cnt, LB))

        """ Update the subproblem """
        # first, get z_sol from updated master problem
        for key in z.keys():
            z_sol[key] = z[key].x

        # change the coefficient of subproblem
        for i in range(facility_num):
            constr_name_1 = 'sub_capacity_1_z_' + str(i)
            constr_name_2 = 'sub_capacity_2_z_' + str(i)
            subProblem.remove(subProblem.getConstrByName(constr_name_1))
            subProblem.remove(subProblem.getConstrByName(constr_name_2))

        # add new constraints
        # cons 1
        for i in range(facility_num):
            expr = LinExpr()
            for j in range(customer_num):
                expr.addTerms(1, x[i, j])
            subProblem.addConstr(expr <= z_sol[i], name='sub_capacity_1_z_%s' % i)

        # logic 1
        for i in range(facility_num):
            subProblem.addConstr(pi[i] <= big_M * v[i])
            expr = LinExpr()
            for j in range(customer_num):
                expr.addTerms(1, x[i, j])
            subProblem.addConstr(z_sol[i] - expr <= big_M - big_M * v[i], name='sub_capacity_2_z_%s' % i)

        """ Update the lower bound """
        subProblem.optimize()
        d_sol = print_sub_sol(subProblem, d, g, x)

        """Update the Upper bound"""
        if (subProblem.status == 2):
            UB = min(UB, subProblem.ObjVal + master.ObjVal - eta.x)
            print('eta = {}'.format(eta.x))
        print('UB (iter {}): {}'.format(iter_cnt, UB))
        Gap = round(100 * (UB - LB) / UB, 2)
        print('---- Gap: {} % ---- '.format(Gap))

    lb.append(LB)
    ub.append(UB)
    kk.append(iter_cnt)

# master.write('finalMP.lp')
print('\n\nOptimal solution found !')
print('Opt_Obj : {}'.format(master.ObjVal))

print(' ** Final Gap: {} % ** '.format(Gap))
print('\n** Solution ** ')
for i in range(facility_num):
    print(' {}: {},\t{}: {} '.format(y[i].varName, y[i].x, z[i].varName, z[i].x), end='')
    print('\t actual demand: {}: {}'.format(d[i].varName, d[i].x), end='')
    print('\t perturbation in worst case: {}: {}'.format(g[i].varName, g[i].x))
print('\n** Transportation solution ** ')
for i in range(facility_num):
    for j in range(customer_num):
        if (x[i, j].x > 0):
            print('trans: {}: {}, cost : {} \t '.format(x[i, j].varName, x[i, j].x,trans_cost[i][j] * x[i, j].x),end='')
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
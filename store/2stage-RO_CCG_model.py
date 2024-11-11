from unicodedata import name
import gurobipy as gp
from gurobipy import GRB
import numpy as np
import scipy.sparse as sp
import seaborn as sns
import matplotlib.pyplot as plt

# Constant Setting
f = [400, 414, 326]  # Coefficients of the objective function for variable y (binary variables)
a = [18, 25, 20]  # Coefficients of the objective function for variable z (continuous variables)
K = [800, 800, 800]
C = [[22, 33, 24],
     [33, 23, 30],
     [20, 25, 27]]  # Cost matrix for x variables
# D = [206+40, 274+40, 220+40]
dl = [206, 274, 220]  # Lower bounds for uncertainty parameters
du = [40, 40, 40]  # Upper bounds for uncertainty parameters
k = 0  # count iteration # Count of iterations for the cutting-plane method

############################# MASTER PROBLEM START ###########################################
MP = gp.Model()  # Creates a Gurobi model.
# Minimize 400 y[0] + 414 y[1] + 326 y[2] + 18 z[0] + 25 z[1] + 20 z[2] + \eta
# Decision variables for the master problem:
x = MP.addVars(3, 3, lb=0, vtype=GRB.CONTINUOUS, name='x_0')  # A 3x3 matrix of continuous variables.
y = MP.addVars(len(f), lb=0, ub=1, obj=f, vtype=GRB.BINARY,name='y')  # Binary variables (0 or 1) with objective coefficients f
z = MP.addVars(len(a), lb=0, obj=a, vtype=GRB.CONTINUOUS,name='z')  # Continuous variables with objective coefficients a.
g = MP.addVars(3, lb=0, ub=1, name='g')  # Continuous variables representing uncertainty in demand ranging from [0,1]
d = MP.addVars(3, lb=0, name='d')  # Continuous variables for demand.
eta = MP.addVar(obj=1, lb=0,
                name='η')  # A continuous variable in the master problem with an objective coefficient of 1.

# Constraints
MP_Cons_1 = MP.addConstrs((z[i] <= K[i] * y[i] for i in range(3)),
                          name='MP_Cons_1')  # For each y[i], it imposes a constraint that z[i] <= 800*y[i].
MP_Cons_2 = MP.addConstr((gp.quicksum(z[i] for i in range(3)) >= 772),
                         name='MP_Cons_2')  # Total z[i] values must be at least 772. (Ensures SP2 has feasible region even at the worst-case)

# iteration constraints
MP_Cons_3 = MP.addConstrs((gp.quicksum(x[i, j] for j in range(3)) <= z[i] for i in range(3)),
                          name='MP_Cons_3')  # For each x[i,j], the sum of x[i,j] values cannot exceed z[i]. Not over inventory
MP_Cons_4 = MP.addConstrs((gp.quicksum(x[i, j] for i in range(3)) >= d[j] for j in range(3)),
                          name='MP_Cons_4')  # For each demand variable d[j], the sum of x[i,j] must be greater than or equal to d[j]. Must satisfy demand
MP_Cons_eta = MP.addConstr(eta >= gp.quicksum(x[i, j] * C[i][j] for i in range(3) for j in range(3)),
                           name='MP_Cons_eta')  # Constraint on eta, requiring it to be greater than or equal to the weighted sum of x[i,j]*C[i,j]. # eta >= b^T x. See note. This is also the cutting plane

# Master-problem uncertainty constraints
MP_Cons_uncertainty_1 = MP.addConstrs((d[i] == dl[i] + du[i] * g[i] for i in range(3)),
                                      name='MP_Uncertainty_Cons1')  # Defines d[i] as dl[i] + du[i]*g[i], capturing uncertainty in demand.
MP_Cons_uncertainty_2 = MP.addConstr((gp.quicksum(g[i] for i in range(3)) <= 1.8), name='MP_Uncertainty_Cons2')
MP_Cons_uncertainty_3 = MP.addConstr((gp.quicksum(g[i] for i in range(2)) <= 1.2),
                                     name='MP_Uncertainty_Cons3')  # Constraints on the sum of g[i] variables. Predefined

MP.optimize()  # Solves the master problem.

LB = MP.objval  # Stores the lower bound from the master problem's objective value.
############################# MASTER PROBLEM END ###########################################


############################# SUB PROBLEM START ###########################################

SP = gp.Model()  # Creates the subproblem model.
# Decision variables for the subproblem, similar to the master problem but with different indices.
x_sub = SP.addVars(3, 3, lb=0, vtype=GRB.CONTINUOUS, name='x_sub')
d_sub = SP.addVars(3, lb=0, name='d_sub')
g_sub = SP.addVars(3, lb=0, ub=1, name='g_sub')
pi = SP.addVars(6, lb=0, vtype=GRB.CONTINUOUS,
                name='pi')  # The number of Lagragian Multiplier is 6, 3 for long constraints, 3 for x bounds constraints
v = SP.addVars(6, vtype=GRB.BINARY, name='v')
w = SP.addVars(3, 3, vtype=GRB.BINARY, name='w')
# A large constant used in the Big-M method for constraints.
M = 10000

# Constraints
# SP_Cons_1 and SP_Cons_2 handle transportation and demand satisfaction.
SP_Cons_1 = SP.addConstrs((gp.quicksum(x_sub[i, j] for j in range(3)) <= z[i].x for i in range(3)),
                          name='SP_Cons_1')  # [0:3] Not over inventory. z[i].x gives the optimized (current solution) value of ( z[i] ) after solving MP. This is a fixed number in the subproblem (SP), as it is based on the outcome of MP
SP_Cons_2 = SP.addConstrs((gp.quicksum(x_sub[i, j] for i in range(3)) >= d_sub[j] for j in range(3)),
                          name='SP_Cons_2')  # [3:6] Meet all demands
SP_Cons_3 = SP.addConstrs((-pi[i] + pi[j + 3] <= C[i][j] for i in range(3) for j in range(3)),
                          name='SP_Cons_3')  # [6:15] SP_Cons_3 deals with pricing constraints based on cost matrix C. \theta_j - \pi_i <= c_{ij}

# slack constraints part 1
# Complemantary slackness condition for /pi and /theta, with big-M and "v&w" for linearization.
SP_SLACK_CONS_1 = SP.addConstrs(
    (z[i].x - gp.quicksum(x_sub[i, j] for j in range(3)) <= M * (1 - v[i]) for i in range(3)),
    name='SP_SLACK_CONS_1')  # [15:18]
SP_SLACK_CONS_2 = SP.addConstrs(
    (gp.quicksum(x_sub[i, j] for i in range(3)) - d_sub[j] <= M * (1 - v[j + 3]) for j in range(3)),
    name='SP_SLACK_CONS_2')  # [18:21]
SP_SLACK_CONS_3 = SP.addConstrs((pi[i] <= M * v[i] for i in range(6)), name='SP_SLACK_CONS_3')  # [21:27]

# slack constraints part 2
# Stationary condition, with big-M and "h" for linearization
SP_SLACK_CONS_4 = SP.addConstrs((C[i][j] + pi[i] - pi[j + 3] <= M * (1 - w[i, j]) for i in range(3) for j in range(3)),
                                name='SP_SLACK_CONS_4')  # [27:36]
SP_SLACK_CONS_5 = SP.addConstrs((x_sub[i, j] <= M * w[i, j] for i in range(3) for j in range(3)),
                                name='SP_SLACK_CONS_5')  # [36:45]

# uncertainty
SP_Cons_uncertainty_1 = SP.addConstrs((d_sub[i] == dl[i] + du[i] * g_sub[i] for i in range(3)),
                                      name='SP_Uncertainty_Cons1')  # [45:48]
SP_Cons_uncertainty_2 = SP.addConstr((gp.quicksum(g_sub[i] for i in range(3)) <= 1.8),
                                     name='MP_Uncertainty_Cons2')  # [48]
SP_Cons_uncertainty_3 = SP.addConstr((gp.quicksum(g_sub[i] for i in range(2)) <= 1.2),
                                     name='MP_Uncertainty_Cons3')  # [49]

sub_obj = gp.quicksum(C[i][j] * x_sub[i, j] for i in range(3) for j in range(3))
SP.setObjective(sub_obj, GRB.MAXIMIZE)
SP.optimize()
SP_objval = SP.objval
UB = LB - eta.x + SP_objval
print(UB)
############################# SUB PROBLEM END ###########################################

lb=[]
ub=[]
kk=[]
lb.append(LB)
ub.append(UB)
kk.append(k)
############################# CCG START ###########################################
while np.abs(UB - LB) > 1e-5:
    print(f"The {k}th iteration")
    k = k + 1
    # Master-problem
    x_new = MP.addVars(3, 3, lb=0, vtype=GRB.CONTINUOUS, name='x_{0}'.format(k))
    MP_Cons_3_new = MP.addConstrs((gp.quicksum(x_new[i, j] for j in range(3)) <= z[i] for i in range(3)),
                                  name='MP_Cons_3_{}'.format(k))  # new "no more than inventory"
    MP_Cons_4_new = MP.addConstrs((gp.quicksum(x_new[i, j] for i in range(3)) >= d_sub[j].x for j in range(3)),
                                  name='MP_Cons_4_{}'.format(k))  # new "demand must be satisfied"
    MP_Cons_eta = MP.addConstr(eta >= gp.quicksum(x_new[i, j] * C[i][j] for i in range(3) for j in range(3)),
                               name='MP_Cons_eta_{}'.format(k))  # new "cutting plane added"
    MP.optimize()
    LB = max(LB, MP.objval)
    print(f"LB: {LB}")

    # Sub-problem update
    # delete old constraints which related to z
    SP.remove(SP.getConstrs()[0:3])
    SP.remove(SP.getConstrs()[15:18])
    # add new constraints which related to z
    SP_Cons_1 = SP.addConstrs((gp.quicksum(x_sub[i, j] for j in range(3)) <= z[i].x for i in range(3)),
                              name='SP_Cons_1')  # [0:3]
    SP_SLACK_CONS_1 = SP.addConstrs(
        (z[i].x - gp.quicksum(x_sub[i, j] for j in range(3)) <= M * (1 - v[i]) for i in range(3)),
        name='SP_SLACK_CONS_1')  # [15:18]

    SP.optimize()
    UB = LB - eta.x + SP.objval
    print(f"UB: {UB}")
    kk.append(k)
    lb.append(LB)
    ub.append(UB)
############################# CCG END ###########################################

# Some information
print("Iteration finished! We found the optimal solution!")
print("Final Objective:{0}".format(LB))
print(y[0], y[1], y[2])
print(z[0], z[1], z[2])
for i in range(3):
    for j in range(3):
        print(x_new[i, j])
print(k)

print(lb,ub,kk)

plt.figure(figsize=(10, 6))

# 绘制 lb 和 ub 的线
sns.lineplot(x=kk, y=lb, label='Lower Bound (lb)', color='blue', marker='o')
sns.lineplot(x=kk, y=ub, label='Upper Bound (ub)', color='orange', marker='o')

# 添加标题和标签
plt.title('Line Plot of lb and ub')
plt.xlabel('k')
plt.ylabel('Values')
plt.legend()

# 显示图形
plt.grid()
plt.show()
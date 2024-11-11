from gurobipy import *
import seaborn as sns
import matplotlib.pyplot as plt

# Problem data (replace with actual data)
f = [400, 414, 326]  # Fixed costs for opening each facility
a = [18, 25, 20]  # Cost coefficients for continuous decision variables z
K = [800,800,800]  # Capacity limit for each facility
c = [[22, 33, 24],
    [33, 23, 30],
    [20, 25, 27]]  # Cost matrix for transportation between facilities and demand points
d_hat = [206, 274, 220]  # Mean demand for each demand point
d_tilde = [40, 40, 40]  # Perturbation for demand points
Gamma = 0.60  # Ambiguity set budget
n = len(d_hat)  # Number of demand points
m = len(f)  # Number of facilities



# Initialize variables for convergence
LB = -float('inf')
UB = float('inf')
tolerance = 1e-4
benders_cuts = []  # List to store cuts

k=1


lb=[]
ub=[]
kk=[]
lb.append(LB)
ub.append(UB)
kk.append(k)

"""Demand Worst Problem"""
DP=Model('Benders decomposition-DP')
g={}
d={}

for j in range(n):
    g[j] = DP.addVar(lb=0,ub=1,vtype=GRB.CONTINUOUS,name=f"g_{j}")
    d[j] = DP.addVar(lb=d_hat[j],vtype=GRB.CONTINUOUS,name=f'd_{j}')

# OBJ
obj=LinExpr()
for j in range(n):
    obj.addTerms(1,d[j])
DP.setObjective(obj,GRB.MAXIMIZE)

# Constraint
expr=LinExpr()
for j in range(n):
    expr.addTerms(1,g[j])
DP.addConstr(expr<=1.8)

expr=LinExpr()
for j in range(2):
    expr.addTerms(1,g[j])
DP.addConstr(expr<=1.2)

for j in range(n):
    expr=LinExpr()
    expr.addTerms(1,d[j])
    expr.addTerms(-d_tilde[j],g[j])
    DP.addConstr(expr==d_hat[j],name=f"Demand_Uncertainty_{j}")

DP.optimize()

d_sol = {}
for key in d.keys():
    d_sol[key] = d[key].X
print(d_sol)

# d_sol={0: 206,1: 314,2:252}

# Step 1: Solve Master Problem (MP)
mp = Model("MasterProblem")

# Variables
y = mp.addVars(m, vtype=GRB.BINARY, name="y")  # Binary variables for facility opening
z = mp.addVars(m, vtype=GRB.CONTINUOUS, lb=0, name="z")  # Continuous variables for allocation
eta = mp.addVar(vtype=GRB.CONTINUOUS, name="eta")  # Eta for Benders cuts
print(f"z: {z}")
# Objective
mp.setObjective(quicksum(f[i] * y[i] + a[i] * z[i] for i in range(m)) + eta, GRB.MINIMIZE)

# Constraints
for i in range(m):
    mp.addConstr(z[i] <= K[i] * y[i], name=f"Capacity_{i}")

# Add Benders cuts from previous iterations
for cut in benders_cuts:
    mp.addConstr(cut)

mp.optimize()

# Get solution from MP
y_values = mp.getAttr('x', y)
z_values = mp.getAttr('x', z)
eta_value = eta.X
LB = max(mp.ObjVal,LB)

# Step 2: Solve Dual Subproblem (Dual SP)
sp = Model("SubProblem")

# Dual variables
lambda_ = sp.addVars(n, vtype=GRB.CONTINUOUS, lb=0, name="lambda")  # Dual for demand constraints
pi = sp.addVars(m, vtype=GRB.CONTINUOUS, lb=0, name="pi")  # Dual for facility capacity


# Objective: Maximize total profit minus penalty
sp.setObjective(quicksum(lambda_[j] * d_sol[j] for j in range(n)) - quicksum(pi[i] * z_values[i] for i in range(m)),GRB.MAXIMIZE)

# Dual constraints
for i in range(m):
    for j in range(n):
        sp.addConstr(c[i][j] + pi[i] - lambda_[j] >= 0, name=f"DualConstraint_{i}_{j}")
#
sp.setParam(GRB.Param.InfUnbdInfo, 1)
sp.optimize()
UB=min(UB,LB-eta.x+sp.ObjVal)

# Benders Decomposition Iterative Process
while abs(UB - LB) > tolerance and k<=20:
    print(benders_cuts)
    k+=1
#
    # Step 3: Check feasibility and add cuts
    if sp.Status == GRB.OPTIMAL:
        mp.reset()
        # Get dual values
        lambda_values = sp.getAttr('x', lambda_)
        pi_values = sp.getAttr('x', pi)
        # d_values = sp.getAttr('x', d)
# #
#         UB = sum(f[i] * y_values[i] + a[i] * z_values[i] for i in range(m)) + sum(lambda_values[j] * d_sol[j] for j in range(n)) - sum(pi_values[i] * z_values[i] for i in range(m))
# #
#         # Create optimality cut
        optimality_cut = (eta >= quicksum(lambda_values[j] * d_sol[j] for j in range(n)) - quicksum(pi_values[i] * z[i] for i in range(m)))
        # benders_cuts.append(optimality_cut)
        mp.addConstr(optimality_cut, name=f"optimality_cut_{k}")
        mp.optimize()
        MP_obj = mp.ObjVal
        LB = max(LB, MP_obj)

        z_sol = {}
        for key in z.keys():
            z_sol[key] = z[key].x

        """SP_dual Reset"""
        sp.reset()

        obj = LinExpr()
        for i in range(m): obj.addTerms(-z_sol[i], pi[i])
        for j in range(n): obj.addTerms(d_sol[j], lambda_[j])

        sp.setObjective(obj, GRB.MAXIMIZE)

        sp.setParam(GRB.Param.InfUnbdInfo, 1)
        # SP_dual.setParam('InfUnbdInfo', 1)
        sp.optimize()
        SP_dual_obj = sp.ObjVal
        UB = min(UB, LB - eta.x + SP_dual_obj)
#
    elif sp.Status == GRB.UNBOUNDED:
        mp.reset()
        # Get dual values for unbounded solution (extreme rays)

        expr = LinExpr()
        cons = 0
        for i in range(m): expr.addTerms(-pi[i].unbdRay, z[i])
        for j in range(n): cons += lambda_[j].unbdRay * d_sol[j]

        # Create feasibility cut
        feasibility_cut = (0 >= cons + expr)

        # benders_cuts.append(feasibility_cut)
        mp.addConstr(0 >= cons + expr, name=f"feasibility_cut_{k}")
        mp.optimize()
        MP_obj = mp.ObjVal
        LB = max(LB, MP_obj)

        z_sol = {}
        for key in z.keys():
            z_sol[key] = z[key].x
        print(f"LB:{LB}")

        """SP_dual Reset"""
        sp.reset()

        obj = LinExpr()
        for i in range(m): obj.addTerms(-z_sol[i], pi[i])
        for j in range(n): obj.addTerms(d_sol[j], lambda_[j])

        sp.setObjective(obj, GRB.MAXIMIZE)

        sp.setParam('InfUnbdInfo', 1)
        sp.optimize()
        SP_dual_obj = sp.ObjVal
        UB = min(UB, LB - eta.x + SP_dual_obj)
        print(f"UB:{UB}")

    lb.append(LB)
    ub.append(UB)
    kk.append(k)

# Solution output
print("Optimal value:", LB)
for i in range(m):
    print(f"Facility {i} - Open: {y[i].x}, Capacity Allocation: {z[i].x}")

print(d_sol,g)


plt.figure(figsize=(10, 6))

sns.lineplot(x=kk, y=lb, label='Lower Bound (lb)', color='blue', marker='o')
sns.lineplot(x=kk, y=ub, label='Upper Bound (ub)', color='orange', marker='o')

plt.title('Line Plot of lb and ub')
plt.xlabel('k')
plt.ylabel('Values')
plt.legend()

plt.grid()
plt.show()
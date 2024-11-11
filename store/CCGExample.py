# -*- coding: utf-8 -*-
"""
Created on Sun Apr  2 17:21:33 2023
@author: wyx
"""
from gurobipy import *
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Parameters
Ay = np.array([[350.,   0.,   0.,  -1.,   0.,   0.],
       [  0., 350.,   0.,   0.,  -1.,   0.],
       [  0.,   0., 350.,   0.,   0.,  -1.],
       [  0.,   0.,   0.,   1.,   1.,   1.]]) # 6x4 # 3x3 for capacity
by = np.array([  0.,   0.,   0., 772.]) # 1x4

G = np.array([[-1., -1., -1.,  0.,  0.,  0.,  0.,  0.,  0.],
       [ 0.,  0.,  0., -1., -1., -1.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  0.,  0.,  0., -1., -1., -1.],
       [ 1.,  0.,  0.,  1.,  0.,  0.,  1.,  0.,  0.],
       [ 0.,  1.,  0.,  0.,  1.,  0.,  0.,  1.,  0.],
       [ 0.,  0.,  1.,  0.,  0.,  1.,  0.,  0.,  1.]]) # 6x9
E = np.array([[0., 0., 0., 1., 0., 0.],
       [0., 0., 0., 0., 1., 0.],
       [0., 0., 0., 0., 0., 1.],
       [0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0.]]) # 6x6
M = np.array([[  0.,   0.,   0.],
       [  0.,   0.,   0.],
       [  0.,   0.,   0.],
       [-40.,   0.,   0.],
       [  0., -40.,   0.],
       [  0.,   0., -40.]]) #6x3
h = np.array([  0.,   0.,   0., 206., 274., 220.]) # 6x1

MP = Model('MP')

# construct main problem
f = np.array([400, 414, 326]) # Coefficients of the objective function for variable y (binary variables), fixed cost
a = np.array([18, 25, 20]) # Coefficients of the objective function for variable z (continuous variables), unit capacity cost
b = np.array([22, 33, 24, 33, 23, 30, 20, 25, 27]) # Cost matrix for x variables, transportation cost
bigM = 10**5
LB = -GRB.INFINITY
UB = GRB.INFINITY
epsilon = 1e-5
k = 1

y = MP.addMVar((3,), obj=f, vtype=GRB.BINARY) # facility location variable, binary
z = MP.addMVar((3,), obj=a, vtype=GRB.CONTINUOUS) # facility capacity hold variable, continuous
d = MP.addMVar((3,), lb=0, name='d') # demand requirements, uncertain
eta = MP.addMVar((1,), obj=1, vtype=GRB.CONTINUOUS)

# construct the original-MP, no eta(cutting plane), no x^l(added columns)
MP.addConstr(Ay[:, :3]@y+Ay[:, 3:]@z >= by) # 4x3@3x1 + 4x3@3x1 >= 4x1 # z_i not exceedings capacity K_i & total holding must have SP feasible.
MP.optimize()
MP_obj = MP.ObjVal
LB = max(MP_obj, LB)


SP = Model('SP')
x = SP.addMVar((9,), vtype=GRB.CONTINUOUS, name='x') # flattened transportation variables
pi = SP.addMVar(G.shape[0], vtype=GRB.CONTINUOUS, name='pi') # Lagrangian Multipliers, numbers equal to the row numbers : i+j
g = SP.addMVar((3,), ub=1, vtype=GRB.CONTINUOUS, name='g') # Uncertainty
v = SP.addMVar((G.shape[0],), vtype=GRB.BINARY, name='v') # big-M, row numbers
w = SP.addMVar((G.shape[1],), vtype=GRB.BINARY, name='w') # big-M, column numbers

G1 = SP.addConstr(G@x >= h-M@g-E@np.concatenate([y.x, z.x]), name="G1") # 6x9@9*1 >= 6x1 -6x3@3x1 - 6x6@(3x1~3x1); fisrt three rows: sumx_ij<=z_i^* inventory limit; last three rows: x_ij>=d_j demand limit
SP.addConstr(G.T@pi <= b, name='pi') # 9x6@6x1<=9x1 duality conditions on Lagragian Multipliers. \theta_j - \pi_i >=0 forall i,j

SP.addConstr(pi <= bigM*v, name='v') # big-M
G2 = SP.addConstr(
    G@x-h+E@np.concatenate([y.x, z.x])+M@g <= bigM*(1-v), name='G2') # i+j i: complementary slackness to \pi z_i^* - sumx_ij <= M(1-v_i). j:complemantary slackness to \theta sumx_ij-d_j<=M(v_j)

SP.addConstr(x <= bigM*w, name='w1') # big-M
SP.addConstr(b-G.T@pi <= bigM*(1-w), name='w2') # (c_ij + \pi_i -\theta_j)x_ij <= bigM(1-w) # second complementary slackness OR first order derivative stationarity

SP.addConstr(g[:2].sum() <= 1.2, name='g1')
SP.addConstr(g.sum() <= 1.8, name='g2')# Uncertainty issue

SP.setObjective(b@x, GRB.MAXIMIZE)
SP.optimize()
SP_obj = SP.ObjVal
UB = min(UB, f@y.x+a@z.x+SP_obj)
MP.reset() # discard solution information

lb=[]
ub=[]
kk=[]
lb.append(LB)
ub.append(UB)
kk.append(k)


while abs(UB-LB) >= epsilon:
    if SP_obj < GRB.INFINITY:
        MP.reset()
        # add x^{k+1}
        x_new = MP.addMVar((9,), vtype=GRB.CONTINUOUS)
        # eta>=bTx^{k+1}
        MP.addConstr(eta >= b.T@x_new) # b: cost matrix
        # Ey+Gx^{k+1}>=h-Mu_{k+1}
        MP.addConstr(E[:, :3]@y+E[:, 3:]@z+G@x_new >= h-M@g.x) # E[:, :3] all rows wanted, first three columns. first three rows: sumx_ij^l<=z_i; last three rows: sumx_ij^l>= d_j^l*
        SP.reset() # new SP because new eta(cutting plane) will be added in
        MP.optimize()
        MP_obj = MP.objval
        LB = max(LB, MP_obj)
    else:
        x_new = MP.addMVar((9,), vtype=GRB.CONTINUOUS)
        MP.addConstr(E[:, :3]@y+E[:, 3:]@z+G@x_new >= h-M@g.x) # add it, but not solve it
    # update the SP constrs according to the MP solution
    SP.remove(G1)
    SP.remove(G2)
    G1 = SP.addConstr(G@x >= h-M@g-E@np.concatenate([y.x, z.x]), name="G1")
    G2 = SP.addConstr(
        G@x-h+E@np.concatenate([y.x, z.x])+M@g <= bigM*(1-v), name='G2')
    SP.optimize()
    # obtain the optimal y^{k+1}
    SP_obj = SP.ObjVal
    UB = min(UB, f@y.x+a@z.x+SP_obj)
    k += 1
    # go back to the MP
    print("经过{}次迭代".format(k))
    print("上界为：{}".format(UB))
    print("下界为：{}".format(LB))
    kk.append(k)
    lb.append(LB)
    ub.append(UB)


print(y[0], y[1], y[2])
print(z[0], z[1], z[2])
print(x)

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

# -*- coding: utf-8 -*-
"""
Created on Nov 5 17:21:33 2024
@author: Ruogu Sun
"""
from gurobipy import *
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# np.random.seed(42)

m=70
n=70
Gamma=0.10
Kapacity=np.random.uniform(500, 700, size=m)
d_bot=np.random.randint(10, 500, size=n)
alpha=np.random.uniform(0.1, 0.5, size=n)
d_tilde=alpha*d_bot

feasibility=sum(d_bot+Gamma*d_tilde)

# Parameters
def generate_Ay(m,Kapacity):
    if len(Kapacity)==m:
        K=np.diag(Kapacity)
        block2=np.diag([-1 for _ in range(m)])
        block3=np.array(([0 for _ in range(m)]+[1 for _ in range(m)],))
        middle=np.concatenate((K,block2)).T
        print(middle.shape,block3.shape)
        return np.concatenate((middle,block3))

    else:raise ValueError
# Ay = np.array([[800.,   0.,   0.,  -1.,   0.,   0.],
#        [  0., 800.,   0.,   0.,  -1.,   0.],
#        [  0.,   0., 800.,   0.,   0.,  -1.],
#        [  0.,   0.,   0.,   1.,   1.,   1.]]) # 6x4 # 3x3 for capacity
Ay=generate_Ay(m=m,Kapacity=Kapacity)

def generate_by(m,feasibility):
    return np.array([0 for _ in range(m)]+[feasibility])
# by = np.array([  0.,   0.,   0., 772.]) # 4x1
by=generate_by(m=m,feasibility=feasibility)

def generate_h(m,n,d_bot):
    if len(d_bot)==n:
        block1=np.zeros(m)
        block2=np.array(d_bot)
        return np.concatenate((block1,block2))
    else:raise ValueError
# h = np.array([  0.,   0.,   0., 206., 274., 220.]) # 6x1
h=generate_h(m=m,n=n,d_bot=d_bot)


def generate_G(m, n):
    block1=[]
    for i in range(m):
        for _ in range(n):
            l1 = [0 for _ in range(m)]
            l1[i]=-1
            block1.append(l1)
    block1=np.array(block1).T

    identity_matrices = [np.eye(n) for _ in range(m)]

    block2 = np.vstack(identity_matrices).T

    return np.concatenate((block1,block2))
# G = np.array([[-1., -1., -1.,  0.,  0.,  0.,  0.,  0.,  0.],
#        [ 0.,  0.,  0., -1., -1., -1.,  0.,  0.,  0.],
#        [ 0.,  0.,  0.,  0.,  0.,  0., -1., -1., -1.],
#        [ 1.,  0.,  0.,  1.,  0.,  0.,  1.,  0.,  0.],
#        [ 0.,  1.,  0.,  0.,  1.,  0.,  0.,  1.,  0.],
#        [ 0.,  0.,  1.,  0.,  0.,  1.,  0.,  0.,  1.]]) # 6x9
G=generate_G(m=m,n=n)


def generate_E(m,n):

    zero_matrix1 = np.zeros((m,m))
    identity_matrices = np.eye(m)
    zero_matrix2 = np.zeros((n,2*m))

    return np.concatenate((np.concatenate((zero_matrix1,identity_matrices)).T,zero_matrix2))
# E = np.array([[0., 0., 0., 1., 0., 0.],
#        [0., 0., 0., 0., 1., 0.],
#        [0., 0., 0., 0., 0., 1.],
#        [0., 0., 0., 0., 0., 0.],
#        [0., 0., 0., 0., 0., 0.],
#        [0., 0., 0., 0., 0., 0.]]) # 6x6
E=generate_E(m=m,n=n)

def generate_M(m,n,d_tilde):
    if len(d_tilde)==n:
        block1=np.zeros((m,n))
        block2=np.diag(d_tilde)
        return np.concatenate((block1,block2))
    else:raise ValueError
# M = np.array([[  0.,   0.,   0.],
#        [  0.,   0.,   0.],
#        [  0.,   0.,   0.],
#        [-40.,   0.,   0.],
#        [  0., -40.,   0.],
#        [  0.,   0., -40.]]) #6x3
M=generate_M(m=m,n=n,d_tilde=d_tilde)



MP = Model('MP')

# construct main problem
def f(m):
    return np.random.uniform(100, 1000, size=m)
# f = np.array([400, 414, 326])# Coefficients of the objective function for variable y (binary variables), fixed cost
f=f(m=m)

def a(m):
    return np.random.uniform(10, 100, size=m)
# a = np.array([18, 25, 20]) # Coefficients of the objective function for variable z (continuous variables), unit capacity cost
a=a(m=m)

def b(mn):
    return np.random.uniform(1, 1000, size=mn)
# b = np.array([22, 33, 24, 33, 23, 30, 20, 25, 27]) # Cost matrix for x variables, transportation cost
b=b(mn=m*n)

bigM = 10**5
LB = -GRB.INFINITY
UB = GRB.INFINITY
epsilon = 1e-2
k = 1

y = MP.addMVar((len(f),), obj=f, vtype=GRB.BINARY) # facility location variable, binary
z = MP.addMVar((len(a),), ub=Kapacity,obj=a, vtype=GRB.CONTINUOUS) # facility capacity hold variable, continuous
# d = MP.addMVar((int(len(b)/len(f)),), lb=0, name='d') # demand requirements, uncertain
eta = MP.addMVar((1,), obj=1, vtype=GRB.CONTINUOUS)

# construct the original-MP, no eta(cutting plane), no x^l(added columns)
MP.addConstr(Ay[:, :len(f)]@y+Ay[:, len(f):]@z >= by) # 4x3@3x1 + 4x3@3x1 >= 4x1 # z_i not exceedings capacity K_i & total holding must have SP feasible.
MP.optimize()
MP_obj = MP.ObjVal
LB = max(MP_obj, LB)


SP = Model('SP')
x = SP.addMVar((len(b),), vtype=GRB.CONTINUOUS, name='x') # flattened transportation variables
pi = SP.addMVar(G.shape[0], vtype=GRB.CONTINUOUS, name='pi') # Lagrangian Multipliers, numbers equal to the row numbers : i+j
g = SP.addMVar((int(len(b)/len(f)),), ub=1, vtype=GRB.CONTINUOUS, name='g') # Uncertainty
v = SP.addMVar((G.shape[0],), vtype=GRB.BINARY, name='v') # big-M, row numbers
w = SP.addMVar((G.shape[1],), vtype=GRB.BINARY, name='w') # big-M, column numbers

G1 = SP.addConstr(G@x >= h-M@g-E@np.concatenate([y.x, z.x]), name="G1") # 6x9@9*1 >= 6x1 -6x3@3x1 - 6x6@(3x1~3x1); fisrt three rows: sumx_ij<=z_i^* inventory limit; last three rows: x_ij>=d_j demand limit
SP.addConstr(G.T@pi <= b, name='pi') # 9x6@6x1<=9x1 duality conditions on Lagragian Multipliers. \theta_j - \pi_i >=0 forall i,j

SP.addConstr(pi <= bigM*v, name='v') # big-M
G2 = SP.addConstr(
    G@x-h+E@np.concatenate([y.x, z.x])+M@g <= bigM*(1-v), name='G2') # i+j i: complementary slackness to \pi z_i^* - sumx_ij <= M(1-v_i). j:complemantary slackness to \theta sumx_ij-d_j<=M(v_j)

SP.addConstr(x <= bigM*w, name='w1') # big-M
SP.addConstr(b-G.T@pi <= bigM*(1-w), name='w2') # (c_ij + \pi_i -\theta_j)x_ij <= bigM(1-w) # second complementary slackness OR first order derivative stationarity

# SP.addConstr(g[:2].sum() <= 1.2, name='g1')
SP.addConstr(g.sum() <= Gamma*n, name='g2')# Uncertainty issue

SP.setObjective(b@x, GRB.MAXIMIZE)
SP.optimize()
SP_obj = SP.ObjVal
UB = min(UB, f@y.x+a@z.x+SP_obj)
# MP.reset() # discard solution information

lb=[]
ub=[]
kk=[]
lb.append(LB)
ub.append(UB)
kk.append(k)



while abs(UB-LB) >= epsilon:
    if SP_obj < GRB.INFINITY:
        # MP.reset()
        # add x^{k+1}
        x_new = MP.addMVar((len(b),), vtype=GRB.CONTINUOUS)
        # eta>=bTx^{k+1}
        MP.addConstr(eta >= b.T@x_new) # b: cost matrix
        # Ey+Gx^{k+1}>=h-Mu_{k+1}
        MP.addConstr(E[:, :len(f)]@y+E[:, len(f):]@z+G@x_new >= h-M@g.x) # E[:, :3] all rows wanted, first three columns. first three rows: sumx_ij^l<=z_i; last three rows: sumx_ij^l>= d_j^l*
        # SP.reset() # new SP because new eta(cutting plane) will be added in
        MP.optimize()
        MP_obj = MP.objval
        LB = max(LB, MP_obj)
    else:
        x_new = MP.addMVar((len(b),), vtype=GRB.CONTINUOUS)
        MP.addConstr(E[:, :len(f)]@y+E[:, len(f):]@z+G@x_new >= h-M@g.x) # add it, but not solve it
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
    print("With {} iterations".format(k))
    print("UB：{}".format(UB))
    print("LB：{}".format(LB))
    kk.append(k)
    lb.append(LB)
    ub.append(UB)


print(y)
print(z)
print(x)

plt.figure(figsize=(10, 6))

sns.lineplot(x=kk, y=lb, label='Lower Bound (lb)', color='blue', marker='o')
sns.lineplot(x=kk, y=ub, label='Upper Bound (ub)', color='orange', marker='o')

plt.title('Line Plot of lb and ub')
plt.xlabel('k')
plt.ylabel('Values')
plt.legend()

plt.grid()
plt.show()
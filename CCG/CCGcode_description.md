# Code Description for CCG in Location-Transportation

## Problem Setting
$m$: potential facility location, indexed by $i$
$n$: customers with different demands, indexed by $j$

$y_i$: binary variable. Whether the $i^{th}$ facility is open.
$z_i$: continuous variable. How many inventory the $i^{th}$ facility currently holds.
$d_j$: continuous variable. How many inventory does $j$ require
$x_{ij}$: continuous variable. How many inventory should be sent from facility $i$ to demand $j$

$g_j$: continuous variable. Uncertainty floating term, ranging from [0,1]
$\pi$: length is $m+n$. Lagrangian multipliers.
$v_i$: length is $m+n$. Binary varibales. Linear Relaxation term for SP's complementary slackness in KKT
$w_i$: length is $m\times n$. Binary varibales. Linear Relaxation term for SP's first-order derivative staionarity in KKT


## Master Problem Initialization
> Note: the example blocks show the numerical case in paper for simplicity and clarity.

```python
Ay = np.array([[800.,   0.,   0.,  -1.,   0.,   0.],
       [  0., 800.,   0.,   0.,  -1.,   0.],
       [  0.,   0., 800.,   0.,   0.,  -1.],
       [  0.,   0.,   0.,   1.,   1.,   1.]])

by = np.array([  0.,   0.,   0., 772.]) # 1x4
```
Ay and by serves for MP formulation. Suppose Ay is $P \times Q$ where $P=m+1$ and $Q=2m$; by is $(m+1)\times 1$
```python
def generate_Ay(m,Kapacity):
    if len(Kapacity)==m:
        K=np.diag(Kapacity)
        block2=-np.diag([-1 for _ in range(m)])
        block3=np.array(([0 for _ in range(m)]+[1 for _ in range(m)],))
        middle=np.concatenate((K,block2)).T
        print(middle.shape,block3.shape)
        return np.concatenate((middle,block3))

    else:raise ValueError
```


The $m\times m$ block is the (K)apacity matrix in diagonal form with $diag(i)=K_i$. $K_i$ stands for the maximal capacity a facility could hold.
The $m\times (m+3)$ block is $diag(i)=-1$
The last row is [0 for _ in range(m)]+[1 for _ in range(m)]
by[-1]: this element should be calculated to ensure SP has feasible region.

```python
f = np.array([400, 414, 326])
a = np.array([18, 25, 20]) 
b = np.array([22, 33, 24, 33, 23, 30, 20, 25, 27]) 
```
f: Fixed Cost. Length is $m$. Coefficients of the objective function for variable y (binary variables)
a: Unit Capacity Cost. Length is $m$. Coefficients of the objective function for variable z (continuous variables)
b: Transportation Cost. Length is $m\times n$. Cost matrix which is flattened for variable x

```python
y = MP.addMVar((m,), obj=f, vtype=GRB.BINARY) # facility location variable, binary
z = MP.addMVar((m,), obj=a, vtype=GRB.CONTINUOUS) # facility capacity hold variable, continuous
d = MP.addMVar((n,), lb=0, name='d') # demand requirements, uncertain
eta = MP.addMVar((1,), obj=1, vtype=GRB.CONTINUOUS)
MP.addConstr(Ay[:, :m]@y+Ay[:, m:]@z >= by) # 4x3@3x1 + 4x3@3x1 >= 4x1 # z_i not exceedings capacity K_i & total holding must have SP feasible.
MP.optimize()
MP_obj = MP.ObjVal
LB = max(MP_obj, LB)
```
Note that for the initial MP, d and eta are not incorporated into the constraints.
y: length should be len(f)=$m$
z: length should be len(a)=$m$
d: length should be $n$

## SubProblem Initialization
```python
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
```
G: shape is $(m+n) \times mn$. Formulation rule:
```python
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
```
E: shape is $(m+n)\times 2m$
```python
def generate_E(m,n):
    zero_matrix1 = np.zeros((m,m))
    identity_matrices = np.eye(m)
    zero_matrix2 = np.zeros((n,2*m))
    return np.concatenate((np.concatenate((zero_matrix1,identity_matrices)).T,zero_matrix2))
```
M: shape is $(m+n)\times n$
```python
def generate_M(m,n,d_tilta:list):
    if len(d_tilta)==n:
        block1=np.zeros((m,n))
        block2=np.diag(d_tilta)
        return np.concatenate((block1,block2))
    else:raise ValueError
```

h: length is $m+n$
```python
def generate_h(m,n,d_bot):
    if len(d_bot)==n:
        block1=np.zeros(m)
        block2=np.array(d_bot)
        return np.concatenate((block1,block2))
    else:raise ValueError
```

## Uncertainty
The value of Γ represents the maximum range of the uncertain demands that can simultaneously deviate from
their nominal values. As the uncertainty is on the right hand sides (demands), Γ belongs to [0, n].
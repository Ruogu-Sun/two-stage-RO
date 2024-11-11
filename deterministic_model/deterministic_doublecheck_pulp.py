import pulp

# Model setup
problem = pulp.LpProblem("Mixed_Integer_Programming", pulp.LpMinimize)

# Decision variables
y = [pulp.LpVariable(f"y_{i}", cat="Binary") for i in range(3)]
z = [pulp.LpVariable(f"z_{i}", lowBound=0) for i in range(3)]
x = [[pulp.LpVariable(f"x_{i}{j}", lowBound=0) for j in range(3)] for i in range(3)]

# objective
problem += (
    400 * y[0] + 414 * y[1] + 326 * y[2] +
    18 * z[0] + 25 * z[1] + 20 * z[2] +
    22 * x[0][0] + 33 * x[0][1] + 24 * x[0][2] +
    33 * x[1][0] + 23 * x[1][1] + 30 * x[1][2] +
    20 * x[2][0] + 25 * x[2][1] + 27 * x[2][2]
)

# supply side constraints
for i in range(3):
    problem += z[i] <= 800 * y[i], f"Supply_Limit_{i}"
    problem += pulp.lpSum(x[i][j] for j in range(3)) <= z[i], f"Production_Limit_{i}"

# demand side constraints
d = [206, 274, 220]
for j in range(3):
    problem += pulp.lpSum(x[i][j] for i in range(3)) >= d[j], f"Demand_Requirement_{j}"

# solve
problem.solve()

# model output
print("Status:", pulp.LpStatus[problem.status])
for var in y + z + [x[i][j] for i in range(3) for j in range(3)]:
    print(f"{var.name} = {var.varValue}")
print("Total Cost:", pulp.value(problem.objective))

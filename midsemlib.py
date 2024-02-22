# Regula Falsi method
def regula_falsi(f, a, b, tol=1e-6, max_iter=1000):
    iterations = 0
    while iterations < max_iter:
        c = (a*f(b) - b*f(a)) / (f(b) - f(a))
        if np.abs(f(c)) < tol:
            return c, iterations
        elif f(c) * f(a) < 0:
            b = c
        else:
            a = c
        iterations += 1
    return None, iterations

# Newton-Raphson method
def newton_raphson(f, df, x0, tol=1e-6, max_iter=1000):
    iterations = 0
    while iterations < max_iter:
        x_next = x0 - f(x0) / df(x0)
        if np.abs(x_next - x0) < tol:
            return x_next, iterations
        x0 = x_next
        iterations += 1
    return None, iterations

# RK4 method for solving the ODEs
def rk4(f, x0, y0, h):
    k1 = h * f(x0, y0)
    k2 = h * f(x0 + 0.5*h, y0 + 0.5*k1)
    k3 = h * f(x0 + 0.5*h, y0 + 0.5*k2)
    k4 = h * f(x0 + h, y0 + k3)
    return y0 + (k1 + 2*k2 + 2*k3 + k4) / 6

# Shooting method to find the initial slope
def shooting_method(slope_guess):
    T_initial = np.array([T_0, slope_guess])
    x_values = np.arange(0, L+dx, dx)

    for x in x_values:
        T_initial = rk4(ode_system, x, T_initial, dx)
        if T_initial[0] >= T_L:
            break

    return x, T_initial[0]

# Binary search to find the correct initial slope
def find_initial_slope():
    min_slope = 0
    max_slope = 100
    tolerance = 1e-6

    while True:
        mid_slope = (min_slope + max_slope) / 2
        x, temperature = shooting_method(mid_slope)
        
        if abs(temperature - T_L) < tolerance:
            return mid_slope, x
        
        if temperature < T_L:
            min_slope = mid_slope
        else:
            max_slope = mid_slope
# Perform LU decomposition without pivoting
def lu_decomposition(A):
    n = len(A)
    L = np.zeros((n, n))
    U = np.zeros((n, n))

    for j in range(n):
        L[j][j] = 1
        
        for i in range(j+1):
            s1 = sum(U[k][j] * L[i][k] for k in range(i))
            U[i][j] = A[i][j] - s1
            
        for i in range(j, n):
            s2 = sum(U[k][j] * L[i][k] for k in range(j))
            L[i][j] = (A[i][j] - s2) / U[j][j]

    return L, U


# Solve Ly = b using forward substitution
def forward_substitution(L, b):
    n = len(b)
    y = np.zeros(n)

    for i in range(n):
        y[i] = b[i]
        for j in range(i):
            y[i] -= L[i][j] * y[j]
    
    return y

# Solve Ux = y using backward substitution
def backward_substitution(U, y):
    n = len(y)
    x = np.zeros(n)

    for i in range(n-1, -1, -1):
        x[i] = y[i]
        for j in range(i+1, n):
            x[i] -= U[i][j] * x[j]
        x[i] /= U[i][i]

    return x
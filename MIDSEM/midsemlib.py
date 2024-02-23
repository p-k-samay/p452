import numpy as np
import matplotlib.pyplot as plt
def newton_raphson(f,d,x):
    ct=0
    while abs(f(x))>10**-6:
        ct+=1
        x= x - f(x)/d(x)
        print(x)
    print("No. of Iterations:",ct)
    print("The root is:",x)
    return
def regula_falsi(f,a,b):
    c=0.0
    ct=0
    c= b - ((b-a)*f(b))/(f(b)-f(a))
    e=10**-6
    while abs(f(c))>e and (b-a)>e and abs(f(a))>e and abs(f(b))>e:
        c= b - ((b-a)*f(b))/(f(b)-f(a))
        if f(a)*f(c)<0:
            b=c
        elif f(b)*f(c)<0:
            a=c
        else:
            a=a+0.1
        ct+=1
        print(c)
    print("No. of Iterations:",ct)
    print("The root is:",c)
    
    def shooter(y_0: float, y_L: float, y_x: float, N: int) -> float:
    z1 = -0.5 
    z1h, y1 = RK4_boundary_value(y_0, z1, 0, y_x, 10, 1) 
    z2 = -2 
    z2h, y2 = RK4_boundary_value(y_0, z2, 0, y_x, 10, 1)  
    iter = 0
    while abs(y1 - y_L) >= 0.001 and abs(y2 - y_L) >= 0.001 and iter <= 30:
        iter += 1
        znew = z2 + ((z1 - z2) * (y_L - y2)) / (y1 - y2)
        znew2, ynew = RK4_boundary_value(y_0, znew, 0, N, 10, 1)  # newton raphson
        if abs(ynew - y_L) < 0.001:
            z, y, x = RK4_boundary_value(y_0, znew, 0, N, 10, 0)  
            break
        else:
            if ynew < y_L:
                z2 = znew
                y2 = ynew
            else:
                z1 = znew
                y1 = ynew
    plt.plot(x, y)
    for i in range(0, len(y)):
        if abs(y[i] - 100) < 0.1:
            out = x[i]
            break
    return out

def RK4_boundary_value(y0, z0,x0, N, end, interpol):
    y_i = y0
    z_i = z0
    step = 0
    yl = [y0]
    zl = [z0]
    xl = [x0]
    h = (end-0)/N
    while step <= end:
        k1y = h*bRK4_dyx(y_i, z_i, step)
        k1z = h*bRK4_dzx(y_i, z_i, step)

        k2y = h*bRK4_dyx(y_i + k1y/2, z_i+k1z/2, step+h/2)
        k2z = h*bRK4_dzx(y_i + k1y/2, z_i+k1z/2, step+h/2)

        k3y = h*bRK4_dyx(y_i+k2y/2, z_i+k2z/2, step+h/2)
        k3z = h*bRK4_dzx(y_i+k2y/2, z_i+k2z/2, step+h/2)

        k4y = h*bRK4_dyx(y_i+k3y, z_i+k3z, step+h)
        k4z = h*bRK4_dzx(y_i+k3y, z_i+k3z, step+h)

        y_i += (k1y+2*k2y+2*k3y+k4y)/6
        z_i += (k1z+2*k2z+2*k3z+k4z)/6
        yl.append(y_i)
        zl.append(z_i)
        step += h
        xl.append(step)

    if interpol == 0:
        return zl, yl, xl
    else:
        return z_i, y_i
    
def bRK4_dyx(y,z,x):
    f = z
    return f

def bRK4_dzx(y,z,x):
    f = 0.01*(y-20)
    return f

def pde(Lx,Nx,Lt,Nt):
    hx=Lx/Nx
    ht=Lt/Nt
    alpha=ht/(hx**2)
    V0=[0]*(Nx+1)
    V1=[0]*(Nx+1)
    I=[]
    print("The value of alpha is", alpha)
    if alpha<0.5:
        print("alpha < 0.5 is satisfied")
    if alpha>0.5:
        print("unstable")
    for i in range(0,Nx+1):
        V0[i]=0.0
        I.append(hx*i)
    V0[int(Nx/2)]=573
    plt.plot(I,V0)
    for j in range(1,1001):
        for i in range(0,Nx+1):
            if i == 0:
                V1[i]=(1-2*alpha)*V0[i]+alpha*V0[i+1]
            elif i==Nx:
                V1[i]=alpha*V0[i-1]+(1-2*alpha)*V0[i]
            else:
                V1[i]=alpha*V0[i-1]+(1-2*alpha)*V0[i]+alpha*V0[i+1]
        for k in range(0,Nx+1):
            V0[k]=V1[k]
        if j==0: 
            plt.plot(I,V0,label='initial parameter')
        if j==10: 
            plt.plot(I,V0,label='10')
        if j==50: 
            plt.plot(I,V0,label='50')
        if j==100: 
            plt.plot(I,V0,label='100')
        if j==500: 
            plt.plot(I,V0,label='500')
        if j==1000: 
            plt.plot(I,V0,label='1000')
    plt.xlabel("Position")
    plt.ylabel("Temperature (K)")
    plt.legend()
    plt.show()


def simpsons(f, a, b, n):
    h = (b - a) / n
    sol = f(a) + f(b)

    for i in range(1, n):
        x_i = a + i * h
        factor = 4 if i % 2 == 1 else 2
        sol += factor * f(x_i)

    sol *= h / 3
    return sol

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

def forward_substitution(L, b):
    n = len(b)
    y = np.zeros(n)

    for i in range(n):
        y[i] = b[i]
        for j in range(i):
            y[i] -= L[i][j] * y[j]
    
    return y

def backward_substitution(U, y):
    n = len(y)
    x = np.zeros(n)

    for i in range(n-1, -1, -1):
        x[i] = y[i]
        for j in range(i+1, n):
            x[i] -= U[i][j] * x[j]
        x[i] /= U[i][i]

    return x
import numpy as np
import scipy as scipy
import math as math
from scipy.optimize import root
  
class Matrix_Operations:
   
    def __init__(self, A):
        self.A = A

    def print_matrix(self,B,flag):
        
        if flag==1:
            for row in B:
                for element in row:
                    print("\t",element, end="\t")
                print()
        elif flag==0:       
            for row in self.A:
                for element in row:
                    print("\t",element, end="\t")
                print()
        else:
            raise ValueError("The flag should be 0 or 1")        


    def matrix_copy(self,A1):
        
        for i in range(len(A1)):
            A1[i]=A1[i][:]
        return A1    



    def add(self,B1):
        
        A1=self.matrix_copy(self.A)
        if len(A1)!=len(B1) or len(A1[0])!=len(B1[0]):
            raise ValueError("The matrices are not of the same order")
        elif():
            pass
        else:
            for i in range(len(A1)):
                for j in range(len(A1[0])):
                    A1[i][j]=A1[i][j]+B1[i][j]
            return A1
        
    def multiply_by_scalar(self,scalar):
        
        A1=self.matrix_copy(self.A)
        for i in range(len(A1)):
            for j in range(len(A1[0])):
                A1[i][j]=A1[i][j]*scalar
        return A1    
        
    def multiply(self,B1):

        A1=self.matrix_copy(B1)
        del B1
        B1=self.matrix_copy(self.A)
        if len(A1[0])!=len(B1):
            raise ValueError("Length of the columns of A should be equal to the length of the rows of B")
        else:
            result=[]
            for i in range(len(A1)):
                row=[]
                for j in range(len(B1[0])):
                    sum=0
                    for k in range(len(A1[0])):
                        sum+=A1[i][k]*B1[k][j]
                    row.append(sum)
                result.append(row)
            return result

  
class Solve_Non_Linear_Equation_1D:
    def __init__(self,f,df,a0,b0,tol):
        self.f = f
        self.a0= a0
        self.b0 = b0
        self.tol = tol
        self.df = df
    def bracket(a,b,f):
        while f(a)*f(b)>=0:
            if abs(f(a))>abs(f(b)):
                b=b+1.5*(b-a)
            else:
                a=a-1.5*(b-a)       
        return(a,b)
    
    def bisection(self):
        a=self.a0
        b=self.b0
        fn=self.f
        epsilon=self.tol
        a,b=Solve_Non_Linear_Equation_1D.bracket(a,b,fn)
        count=0
        while (abs(b-a))>epsilon:
            c=(a+b)/2
            if fn(a)*fn(c)>0:
                a=c
            else:
                b = c 
            count+=1       
        return c,count 
    
    def secant(self,guess1,guess2):
        fn=self.f
        t=self.tol
        x0=guess1
        x1=guess2
        x2=x1-((fn(x1)*(x1-x0))/(fn(x1)-fn(x0)))
        step=1
        while abs(x2-x1)>t:
            if step>100:
                raise ValueError("The roots are not converging")
                break
            else:
                x0=x1
                x1=x2
                x2=x1-fn(x1)*(x1-x0)/(fn(x1)-fn(x0))
                step+=1
        return x2,step 
    
    def fixed_point(self,x0):
        g=self.f
        tol=self.tol
        x1=g(x0)
        step=1
        while abs(x1-x0)>tol:
            if step>100:
                print("The roots are not converging")
                break
            else:
                x0=x1
                x1=g(x0)
                step+=1
        return x1,step       

class Newton_Cotes:
    def __init__(self, f, a, b,N):
        self.a = a
        self.b = b
        self.N = N
        self.f = f

    def midpoint(self):
        h=(self.b-self.a)/self.N
        I=0
        x=self.a
        while x<=self.b:
            I+=h*self.f(x+h/2)
            x+=h
        return I
    
    def trapezoidal(self):
        h=(self.b-self.a)/self.N
        I=0
        x=self.a
        while x<=self.b:
            I+=h/2*(self.f(x)+self.f(x+h))
            x+=h
        return I


    def simpsons(self):
        h=(self.b-self.a)/self.N
        I=0
        x=self.a
        while x<=self.b:
            I+=h/3*(self.f(x)+4*self.f(x+h/2)+self.f(x+h))*0.5
            x+=h
        return I

class Gaussian_Quadrature: 
    def __init__(self,f,a,b,degree):
        self.a = a
        self.b = b
        self.N = degree
        self.f = f   

    def Pn(self,x,n):
        if n == 0:
            return 1
        elif n == 1:
            return x
        else:
            return ((2*n-1)*x*self.Pn(x,n-1)-(n-1)*self.Pn(x,n-2))/n

    def Pn_drvtv(self,x,n):
        if n == 0:
            return 0
        elif n == 1:
            return 1
        else:
            return (n*(x*self.Pn(x,n)-self.Pn(x,n-1)))/(1-x**2)
    
    def find_legendre_roots(self):
        n=self.N
        num_roots=self.N
        roots = []
        for i in range(1, num_roots + 1):
            guess = np.cos((2*i - 1) * np.pi / (2 * num_roots))
            result = root(self.Pn, guess, args=(n,), jac=self.Pn_drvtv, method='hybr')
            if result.success:
                roots.append(result.x[0])
        return roots
    

    def find_weights(self):
        n=self.N
        roots=self.find_legendre_roots()
        weights=[]
        for i in range(n):
            w=2/((1-roots[i]**2)*(self.Pn_drvtv(roots[i],n))**2)
            weights.append(w)
        return weights


    def integrate(self):
        a=self.a
        b=self.b
        n=self.N
        f=self.f
        sum=0
        weights=self.find_weights()
        roots=self.find_legendre_roots()
        for i in range(n):
            y=((b-a)*0.5*roots[i])+((b+a)*.5)
            weightf=weights[i]*f(y)
            sum+=weightf
        val=(b-a)*0.5*sum
        return val       

class ODE_Solve_XY:
    def __init__(self, dy, xi,yi,xf,N):
        self.xi = xi
        self.yi = yi
        self.xf = xf
        self.N = N
        self.dy = dy
        
    def forward_euler(self):
        x_ini=self.yi
        t_ini=self.xi
        t_final=self.xf
        n=self.N
        dx=self.dy
        dt=(t_final-t_ini)/n
        xlist=[]
        ylist=[]
        t=t_ini
        while t<=t_final:
            xlist.append(t)
            ylist.append(x_ini)
            x_ini+=dt*dx(t,x_ini)
            t+=dt
        return xlist,ylist        
    
    def backward_euler(self):
        x0=self.xi
        y0=self.yi
        xf=self.xf
        num_points=self.N
        fn=self.dy

        h = (xf - x0) / num_points
        x_values = np.linspace(x0, xf, num_points + 1)
        y_values = np.zeros(num_points + 1)
        y_values[0] = y0

        for i in range(1, num_points + 1):
            y_values[i] = y_values[i - 1] + h * fn(x_values[i], y_values[i - 1])

        return x_values, y_values

    def predictor_corrector(self):
        dybydx=self.dy
        x0=self.xi
        y0=self.yi
        x_f=self.xf
        n=self.N

        h=(x_f-x0)/n
        xlist=[]
        ylist=[]
        x=x0
        y=y0
        xlist.append(x)
        ylist.append(y)
        while x<x_f:
            k1=dybydx(x,y)*h
            k2=dybydx(x+h,y+k1)*h
            y=y+0.5*(k1+k2)
            x=x+h
            xlist.append(x)
            ylist.append(y)
        return xlist,ylist
    
    def RK2_solve(self):

        dybydx=self.dy
        x0=self.xi
        y0=self.yi
        xf=self.xf
        N=self.N

        h=(xf-x0)/N
        xlist=[]
        ylist=[]
        x=x0
        y=y0
        xlist.append(x)
        ylist.append(y)
        while x<xf:
            k1=h*dybydx(x,y)
            k2=dybydx(x+(h/2),y+(k1/2))*h
            y=y+k2
            x=x+h
            xlist.append(x)
            ylist.append(y)
        return xlist,ylist
    
    def RK4_solve(self):

        dybydx=self.dy
        x0=self.xi
        y0=self.yi
        x_f=self.xf
        N=self.N

        h=(x_f-x0)/N
        xlist=[]
        ylist=[]
        x=x0
        y=y0
        xlist.append(x)
        ylist.append(y)
        while x<x_f:
            k1=h*dybydx(x,y)
            k2=h*dybydx(x+(h/2),y+(k1/2))
            k3=h*dybydx(x+(h/2),y+(k2/2))
            k4=h*dybydx(x+h,y+k3)
            y=y+(k1+2*k2+2*k3+k4)/6
            x=x+h
            xlist.append(x)
            ylist.append(y)
        return xlist,ylist     

def semi_implicit_euler_solve(f,g,x0,v0,t0,t_max,step_size):
   
    h=step_size
    vlist=[]
    xlist=[]
    tlist=[]
    x=x0
    v=v0
    t=t0
    while t<=t_max:
        xlist.append(x)
        vlist.append(v)
        tlist.append(t)
        v=v + (h*g(x,t))
        x=x + (h*f(v,t))
        t=t+h
    return xlist,vlist,tlist    

class verlet_algorithm:

    def __init__(self,a, x0, v0, t0, t_max, step_size):
        self.a = a
        self.x0 = x0
        self.v0 = v0
        self.t0 = t0
        self.t_max = t_max
        self.step_size = step_size

    def verlet_solve(self):
        #Defining the variables
        h=self.step_size
        xlist=[]
        vlist=[]
        tlist=[]
        tm=self.t_max
        x=self.x0
        t=self.t0
        v=self.v0
        acc_fn=self.a

        #The first 
        xlist.append(x)
        vlist.append(v)
        tlist.append(t)
        x1=(x)+(h*v)+(0.5*h*h*self.a(x))
        v1=(x1-x)/h
        t=t+h
        xlist.append(x1)
        vlist.append(v1)
        tlist.append(t)

        #The rest of the steps
        while t<=tm:
            x2=(2*x1)-(x)+(h*h*acc_fn(x1))
            v=(x2-x)/(2*h)
            x=x1
            x1=x2
            t=t+h
            xlist.append(x)
            vlist.append(v)
            tlist.append(t)

        return xlist,vlist,tlist    

    def velocity_verlet(self):
        
        h=self.step_size
        x=self.x0
        v=self.v0
        t=self.t0

        xlist=[x]
        vlist=[v]
        tlist=[t]
        pass


class Explicit_solve_uxx_ut: #not complete
    def __init__(self,g,a,b,L,nx,T,nt,timestep):
   
        self.g = g
        self.a = a
        self.b = b
        self.L = L
        self.nx = nx
        self.nt = nt
        self.T = T
        self.timestep = timestep
    def check_stability(self):
 
        alpha = (self.T*(self.nx**2))/((self.L**2)*(self.nt))
        if alpha>0.5:
            raise ValueError("The method is not stable")
        else:
            pass

    def create_A_inv(self):
     
        self.check_stability()

        alpha = (self.T*(self.nx**2))/((self.L**2)*(self.nt))
        A=[]
        for i in range(self.nx):
            row=[]
            for j in range(self.nx):
                if i==j:
                    row.append(1+2*alpha)
                elif abs(i-j)==1:
                    row.append(-alpha)
                else:
                    row.append(0)
            A.append(row)
        A_inv = np.linalg.inv(A)
        return A_inv 

    def create_V0(self):
        f=self.g
        V0=[]
        for i in range(self.nx):
            if i==0:
                V0.append([self.a(0)])
            elif i==self.nx-1:
                V0.append([self.b(0)])
            else:
                V0.append([f(i*self.L/self.nx)])   
        return V0   

    def solve(self):     
        A=self.create_A()
        V0=self.create_V0()
        A_inv = np.linalg.inv(A)
        for i in range(self.timestep):
            V1=Matrix_Operations(V0).multiply(A_inv)
            # V1=np.matmul(A_inv,V0)
            del V0
            V0=V1
            del V1
            V0[0][0]=0
            V0[self.nx-1][0]=0
        return V0


        
def crank_nicholson(f, x_0, x_N, N_x, N_t, T, alpha):

    x = np.linspace(x_0, x_N, N_x+1)
    t = np.linspace(0, T, N_t+1)               
    u = []
    u.append([f(x[i]) for i in range(N_x+1)])  
    u[0][0] = 0
    u[0][N_x] = 0                          
    
    u = np.array(u)
    I2paB_inv_cols = []                        

    for j in range(N_x + 1):                  
        X = list([0] for i in range(N_x + 1))  
        for step in range(150):
            flag = 1
            for i in range(N_x + 1):
                sum = 0
                if i != 0:
                    sum += (- alpha * X[i - 1][0])
                if i != N_x:
                    sum += (- alpha * X[i + 1][0])
                if i == j:
                    temp = (1-sum) / (2+(2*alpha))
                else:
                    temp = (-sum) / (2+(2*alpha))
                if abs((temp) - (X[i][0])) > 1e-6: 
                    flag = 0
                X[i][0] = temp
            if flag == 1:
                break
        if flag == 0:
            print(f'Eqn not solved after 150 steps for {j}th col')
            return None
        I2paB_inv_cols.append(X)
    
    I2paB_inv = np.array(I2paB_inv_cols[0])
    for i in range(1,N_x + 1):
        I2paB_inv = np.append(I2paB_inv,I2paB_inv_cols[i],axis = 1)
    
    I2paB_inv_I2naB = []
    
    for row in range(N_x + 1):
        I2paB_inv_I2naB.append([])
        for col in range(N_x + 1):
            sum = 0
            sum += I2paB_inv[row][col] * 2 * (1 - alpha)
            if col != 0:
                sum += I2paB_inv[row][col - 1] * alpha
            if col != N_x:
                sum += I2paB_inv[row][col + 1] * alpha
            I2paB_inv_I2naB[row].append(sum)
    I2paB_inv_I2naB = np.array(I2paB_inv_I2naB)
    
    
    del I2paB_inv_cols, I2paB_inv
    
    for n in range(N_t):
        u = np.append(u, [I2paB_inv_I2naB @ u[-1,:]], axis = 0) 
        u[-1,0] = 0
        u[-1,N_x] = 0
  
    return x,t,u



def poisson_laplace(rho, x_i, x_f, y_i, y_f, u_iy, u_fy, u_xi, u_xf, N):
   
    x = np.linspace(x_i, x_f, N+2)
    y = np.linspace(y_i, y_f, N+2)
    hx = (x_f - x_i)/(N + 1)
    hy = (y_f - y_i)/(N + 1)
    if hx != hy:
        raise ValueError("The grid is not square")
    h = hx 
    A = np.zeros((N**2,N**2))
    
    for i in range(N**2):
        A[i, i] = 4
        if i == 0:
            A[i, i+1]=-1
            A[i, i+N]=-1
        elif i < N:
            A[i, i-1]=-1
            A[i, i+1]=-1
            A[i, i+N]=-1
        elif i < (N**2-N):
            A[i, i-1]=-1
            A[i, i+1]=-1
            A[i, i-N]=-1
            A[i, i+N]=-1
        elif i < (N**2-1):
            A[i, i-1]=-1
            A[i, i+1]=-1
            A[i, i-N]=-1
        else:
            A[i, i-1] = -1
            A[i, i-N] = -1
   
    B = []
    for i in range(1,N+1):
        for j in range(1,N+1):
            sum = rho(x[i], y[j]) * h**2
            if i == 0:
                sum += u_xi(y[j])
            if i == N:
                sum += u_xf(y[j])
            if j == 0:
                sum += u_iy(x[i])
            if j == N:
                sum += u_fy(x[i])
            B.append(sum)    
    B = np.array(B)[:,None]
    u = np.linalg.solve(A, B)
    u = np.array(u).reshape((N,N))
    u = np.append(u_iy(y[1:-1,None]), u, axis = 1)
    u = np.append(u, u_fy(y[1:-1, None]), axis = 1)
    u = np.append([u_xi(x)], u, axis = 0)
    u = np.append(u, [u_xf(x)], axis = 0)
    

   
    return x, y, u






                               

                    



        








        












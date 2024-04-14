def Pow(A,x0,y):
    Ax = [0,0,0,0]
    for i in range(len(x0)):
        for j in range(len(A)):
            Ax[i]+=(A[i][j])*x0[j]

    Axy=0
    for i in range(len(x0)):
        Axy+=(Ax[i])*(y[i])
    return Ax, Axy

def find(A,x0,y,e):
    mult=A
    k=0
    evl0 = 2
    v, evl1=Pow(mult,x0,y)
    vl0 = evl0
    vl1 = evl1 / evl0

    while abs(vl1-vl0)>e and k<50:
        k+=1
        evl0=evl1
        mult1=[[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
        for i in range(len(A)):
            for j in range(len(A)):
                for k in range(len(A)):
                    mult1[i][j]+=(mult[i][k])*(A[k][j])
        mult=mult1
        v, evl1=Pow(mult,x0,y)
        vl0 = vl1
        vl1=evl1/evl0
    sum=0

    for i in v:
        sum += i*i
    sum = math.sqrt(sum)
    for i in range(len(v)):
        v[i] =v[i]/sum
    return vl1, v, k
def QR_SCD(A):
    
    Q =  np.array([ [0  for i in range(len(A))] for j in range(len(A))])
    R = ([ [0  for i in range(len(A))] for j in range(len(A))])
    e = np.array( [0  for i in range(len(A))])
    
    for i in range(len(A)):
        substraction = 0
        for k in range(i):
            substraction +=  np.dot(np.dot( A[:,i],Q[:,k]),Q[:,k])
        ui =  copy.deepcopy(A[i,:] - substraction)
        norm_ui = np.linalg.norm(ui) 
        e = copy.deepcopy(ui / norm_ui)
        for j in range(len(A)):
            Q = Q.tolist() 
            e = e.tolist()                    
            Q[j][i] =e[j]
            Q = np.array(Q)
            e = np.array(e)             
        k = i
        while k < (len(A)):
            R[i][k] = np.dot(A[:,k],Q[:,i])
            k+=1      
    R = np.array(R)
    return Q,R

def find_eigenvalues(A):
    sol_ev_qr = copy.deepcopy(A)
    evalues_QR = []
    
    for i in range(20):
        Q, R = QR_SCD(sol_ev_qr)
        sol_ev_qr = np.dot(R, Q)
    for i in range(len(A)):
        temp_ev_QR = sol_ev_qr[i][i]
        evalues_QR.append(round(temp_ev_QR,6))
    return evalues_QR
def InverseGJ(l):
    for i in range(len(l)):
        for j in range(len(l),2*len(l)):
            if j==(i+len(l)):
                l[i].append(1)
            else:
                l[i].append(0)
    for i in range(len(l)-1,0,-1):
        if l[i-1][0]<l[i][0]:
            l[i-1],l[i]=l[i],l[i-1]
            
    for i in range(len(l)):
        for j in range(len(l)):
            if i!=j:
                p=l[j][i]/l[i][i]
                for k in range(2*len(l)):
                    l[j][k]=l[j][k]-p*l[i][k]
    for i in range(len(l)):
        p=l[i][i]
        for j in range(2*len(l)):
            l[i][j]=l[i][j]/p

    l_inv=[[0 for col in range(len(l))] for row in range(len(l))]
    for i in range(len(l)):
        for j in range(len(l),2*len(l)):
            l_inv[i][j-len(l)]=l[i][j]
    return l_inv


def LeastSquareFit(l):
    sum=0
    r_pt,h_pt=[],[]
    for i in range(len(l)):
        r_pt.append(l[i][0])
        h_pt.append(l[i][1])
    plt.scatter(r_pt,h_pt)
    x=[[0 for col in range(len(l[0])+2)] for row in range(len(l[0])+2)]
    y=[[0 for col in range(len(l[0])+2)] for row in range(len(l[0])+2)]
    
    for m in range(len(l[0])+2):
        for j in range(m,len(l[0])+m+2):
            sum=0
            for k in range(len(l)):
                sum=sum+(l[k][0])**j           
            x[m][j-m]=sum
    
    for i in range(len(l[0])+2):
        sum=0
        for k in range(len(l)):
            sum=sum+(l[k][1])*(l[k][0])**i
        y[i]=sum

    x_inv = InverseGJ(x)
    
    a=[0 for col in range(len(x_inv))]
    
    for i in range(len(x_inv)):
        for j in range(len(y)):
            a[i]=a[i]+(x_inv[i][j])*(y[j])
    
    print("Coefficients",a)
    x_plot, y_plot = [], []
    for i in range(len(l)):
        sum=0
        for j in range(len(a)):
            sum=sum+(a[j])*((l[i][0])**j)
        y_plot.append(sum)
        x_plot.append(l[i][0])
        
    plt.plot(x_plot,y_plot)
    plt.title("Least square fit with f(x)=a_0+a_1x+a_2x^2+a_3x^3")
    plt.ylabel("y axis")
    plt.xlabel("x axis")
    plt.show()
    return a
def UandL(A):
    u=[[0 for col in range(len(A))] for row in range(len(A))]
    l=[[0 for col in range(len(A))] for row in range(len(A))]
    for i in range(len(A)):
        u[0][i]=A[0][i]
        l[i][i]=1
    for j in range(len(A)):
        for i in range(1,j+1):
            sum=0
            for k in range(i):
                sum=sum+(l[i][k])*(u[k][j])
            u[i][j]=A[i][j]-sum

        for i in range(j,len(A)):
            sum=0
            for k in range(j):
                sum=sum+(l[i][k])*(u[k][j])
            if u[j][j]==0:
                continue
            else:
                l[i][j]=(A[i][j]-sum)/(u[j][j])
    return u,l

def LUForBack(l,U,L):
    y=[0 for i in range(len(L))]
    x=[0 for i in range(len(L))]
    y[0]=l[0]


    for i in range(1,len(L)):
        sum=0
        for j in range(i):
            sum=sum+L[i][j]*y[j]
        y[i]=(l[i]-sum)/L[i][i]
    x[len(L)-1]=y[len(L)-1]/U[len(L)-1][len(L)-1]
    
    for i in range(len(L)-2,-1,-1):
        sum=0
        for j in range(len(L)):
            sum=sum+U[i][j]*x[j]
        x[i]=(y[i]-sum)/U[i][i]
    return(x)
def chebyshev(x,order):
    if order==0:
        return np.ones_like(x)
    elif order==1:
        return (2*x)-1
    elif order==2:
        return (8*(x**2))-(8*x)+1
    elif order==3:
        return (32*(x**3)-(48*(x**2))+(18*x)-1)
    elif order==4:
        return (128*(x**4)-256*(x**3)+160*(x**2)-(32*x)+1)
    
def fit_chebyshev(l,order):
    x,y=[],[]
    for i in range(len(l)):
        x.append(l[i][0])
        y.append(l[i][1])
    parameters=order+1
    A=np.zeros((parameters,parameters))
    b = np.zeros(parameters)
    n=len(x)
    for i in range(parameters):
        for j in range(parameters):
            total=0
            for k in range(n):
                total+=chebyshev(x[k],j)*chebyshev(x[k],i)
            A[i,j]=total
    for i in range(parameters):
        total=0
        for k in range(n):
            total += chebyshev(x[k], i)*y[k]
        b[i]=total
    u,l=UandL(A)
    coeff=LUForBack(b,u,l)
    return coeff,A

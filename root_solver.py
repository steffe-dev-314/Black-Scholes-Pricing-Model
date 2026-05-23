import numpy as np 
from typing import Callable

def do_newton(f, fprime , x0 , args , maxiter , tol):

    run_num =0 
    while np.abs(f(x0 , *args)) > tol:
        x0 = x0 - f(x0 , *args)/fprime(x0 , *args)
        run_num += 1 
        if run_num == maxiter:
            break

def find_sign_change(f : Callable,args : tuple,a0:float , maxiter:int = 100, stepsize:float = 1):
    run_num = 0
    b0 = a0 + stepsize 

    # more stable and <= indicates sign chagne and near 0 too 
    while f(a0,*args)*f(b0,*args)>0:
        a0 = b0 
        b0 = b0 + stepsize
        run_num += 1 
        if run_num == maxiter:
            print('maxiter reached')
            break 
    return a0 , b0 

#stepsize as % of a0? 
def do_bisection(f:Callable , args : tuple , a0:float, maxiter:int = 100 , stepsize:float = 1 ,root_tol :float = 1e-6):
    a , b = find_sign_change(f,args,a0,maxiter,stepsize)
    print(f'initial a,b : ({a:.4f},{b:.4f})')
    c = (a+b)/2
    runnum = 0 
    while f(c,*args) >= root_tol:

        if f(a,*args)*f(c,*args)<=0:
            a = a 
            b = c 
        elif f(b,*args)*f(c,*args)<=0:
            a = c
            b = b 
        c = (a+b)/2 
        runnum +=1  
        if runnum == maxiter:
            break
    return c 
            
        
if __name__ == '__main__':


    def my_func(x):
        return x**3 - x # x*(x-1)*(x+1) 

    def der_my_func(x): 
        return 3*x**2 - 1

    def sign_func(x):
        return x*(x-2)*(x+2)*np.sqrt(x+1)

    a,b = find_sign_change(f = sign_func , args = () , a0 = 0.3 ,stepsize=-0.1)
    #print(np.round(a,2),np.round(b,2))
    #print(np.round(sign_func(a),2) , np.round(sign_func(b),2))
    print(do_bisection(f = sign_func , args = () , a0 = -0.5 , maxiter = 100 , stepsize= 1, root_tol=1e-6))

##### procedure ##### 
# f(xi) = x(i-1) - f(x(i-1))/f'(x(i-1))

## for IV:
# do this in class 
# def f(sigma) := C(sigma, args) - P_market
# def f'(sigma) := C'(sigma,args) 
# pass f , f' to newton OR 
# iter through sigma range, say 0 , 1 in 0.1 steps and find sign change --> pass interval , f , f' to bisect
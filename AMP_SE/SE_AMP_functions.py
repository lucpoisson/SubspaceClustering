import numpy as np
import random
from scipy import integrate
from scipy.optimize import fsolve
from scipy.signal import find_peaks


#### SE ###
def damp(x_new, x_old, damping=0.5):
    return damping * x_old + (1-damping) * x_new
def fv(A,B,rho):
    if rho == 1:
        ret = B/(1+A)
    else:
        den = (rho + (1-rho)*np.sqrt(1+A)*np.exp(-0.5*B**2/(A+1)))*(1+A)
        ret = rho*B/den
    return ret
def fu(B):
    return np.tanh(B)
def gaussian(z,mean=0,var=1):
    num = np.exp(-0.5*(z-mean)**2/var)
    den = np.sqrt(2*np.pi*var)
    ret = num/den
    return ret
def update_u(x):
    f = lambda w: 0.5*(fu(B=x + np.sqrt(x)*w)-fu(B=-x + np.sqrt(x)*w))*gaussian(z=w)
    I = integrate.quad(f, -np.inf,np.inf)
    return I[0]

def update_v(x,rho):
    f = lambda z: rho*gaussian(z=z,mean=0,var =1)*z*fv(A=x,B=np.sqrt(x*(1+x))*z,rho=rho)*np.sqrt(x**2/(x**2+x))
    I = integrate.quad(f, -np.inf,np.inf)
    return I[0]

def iterate_UV(init_u,init_v, tol, max_steps , *, snr: float, rho:float,alpha:float,damping:float):
    Mu = np.zeros(max_steps) ; Mv = np.zeros(max_steps)
    overlap_u = np.zeros(max_steps) ; mse_v = np.zeros(max_steps)
    Mu[0] = init_u ; Mv[0] = init_v ;   c=0
    print(f"Welcome to SE with snr={snr} and rho={rho}")
    for t in range(max_steps-1):
        xv = snr*Mv[t] ; Mutmp = update_u(x=xv)
        xu =   alpha*snr*Mutmp  ; Mvtmp = update_v(x=xu,rho=rho)
        diffu = np.sum(np.abs(Mutmp-Mu[t])) ; diffv = np.sum(np.abs(Mvtmp - Mv[t]))
        print(f"at time {t} we have deltau = {diffu} and deltav={diffv} and magnetization = {Mutmp} and Mv={Mvtmp} ")
        flag = (diffu + diffv < tol)
        # damping
        Mu[t+1] = damp(Mutmp,Mu[t],damping) ; Mv[t+1] = damp(Mvtmp,Mv[t],damping)
        if flag:
            c = 1
            break
    if  c:
        print(f"SE converged with well defined magnetization = {Mu[t]}")
    mse_cent = rho - Mv[t]
    mseu = 0.5 - Mu[t]
    return  mse_cent , Mv[t],  mseu , Mu[t]  


### AMP ###
def lowRAMP(truth,init, tol, max_steps, *, snr: float, rho:float,damping:float):
    u0,v0,X = truth   ; u, v = init
    d = len(v) ; n = len(u) ; alpha = n/d
    uold, vold = np.zeros(n), np.zeros(d) ; Bu, Bv = np.zeros(n), np.zeros(d)
    Au, Av  = 0,  0  ; sigmau, sigmav = np.zeros(n), np.zeros(d)
    conv = 0
    print(f"Welcome to LOWRAMP with snr={snr} and rho = {rho}")
    for iter_ in range(max_steps):

        Bvold = np.copy(Bv)  ;  Buold = np.copy(Bu)
        Avold = np.copy(Av)  ;  Auold = np.copy(Au)

        Bvnew = update_Bv(X,u,vold, sigmau, snr)
        Avnew =  update_Av(X,u, v, snr)
        Bv = damping * Bvold  + (1-damping) * Bvnew   ;  Av = damping * Avold  + (1-damping) * Avnew

        vold = np.copy(v)
        v  = evolve_v(Av, Bv, rho = rho)
        sigmav = evolve_sigmav(Av, Bv , rho = rho)


        Bunew = update_Bu(X, v, uold, sigmav, snr)
        Aunew =  update_Au(X, v, snr)
        Bu = damping * Buold  + (1-damping) * Bunew   ;  Au = damping * Auold  + (1-damping) * Aunew

        uold = np.copy(u)
        u    = evolve_u(Aunew, Bunew)
        sigmau = evolve_sigmau(Aunew, Bunew)

        diff = np.mean(((vold-v)**2)) + np.mean(((uold-u)**2))
        Mu = np.dot(u,u0)/n  ; Mv = np.dot(v,v0)/d
        if iter_ >=1:
            print(f"at iter = {iter_} we have diff = {diff} and magnetization={Mu}")
        if diff < tol:
            print("AMP converged")
            conv = 1
            break
        mseu = 0.5 - np.abs(Mu)
    return mseu , conv

def get_instance(rho, dimension, nsamples, snr):
    u0 = np.sign(np.random.randn(nsamples))
    v0 = np.random.normal(0,1,dimension)
    elim = np.random.binomial(n=1,p=rho,size =dimension)
    v0 *= elim
    noise =  np.random.randn(dimension,nsamples)
    Y = np.sqrt(snr/dimension) * np.outer(v0,u0) + noise
    return u0, v0, Y
def fprime_v(A,B,rho):
    if rho !=1:
        den = (rho + (1-rho)*np.sqrt(1+A)*np.exp(-0.5*B**2/(A+1)))*(1+A)
        der = -(1+A)*((1-rho)*B/np.sqrt(1+A)*np.exp(-0.5*B**2/(1+A)))
        ret = (rho*den - rho*B*der)/(den**2)
#         print(ret , 1/(1+A))
        return  ret
    if rho == 1:
        return 1/(1+A)
def fprime_u(B):
    ret = 1 - np.tanh(B)**2
    return ret
def update_Bv(X,u,v, sigmau, snr):
    d = len(v) ; n = len(u)
    term1 = np.sqrt(snr/d)*X@u ; term2 = snr/d*np.sum(sigmau)*v
    ret = term1 - term2
    return ret
def update_Av(X, u,v,snr):
    n = len(u) ; d = len(v)
    ret = snr/d*np.sum(u**2)
    return ret
def evolve_v(Av, Bv, rho):
    d = len(Bv) ; ret = np.zeros(d)
    for i in range(d):
        ret[i] = fv(Av,Bv[i],rho)
    return ret
def evolve_sigmav(Av, Bv , rho):
    d = len(Bv) ; ret = np.zeros(d)
    for i in range(d):
        ret[i] = fprime_v(Av,Bv[i],rho)
    return ret
def update_Bu(X, v, u, sigmav, snr):
    d = len(v) ; n = len(u)
    ret = np.sqrt(snr/d)*X.T@v - snr/d*np.sum(sigmav)*u
    return ret
def update_Au(X, v, snr):
    d = len(v)
    ret = snr/d*np.sum(v**2)
    return ret
def evolve_u(Au, Bu):
    n = len(Bu) ; ret = np.zeros(n)
    for i in range(n):
        ret[i] = fu(Bu[i])
    return ret
def evolve_sigmau(Au, Bu):
    n = len(Bu) ; ret = np.zeros(n)
    for i in range(n):
        ret[i] = fprime_u(Bu[i])
    return ret


### Thresholds ###
def aux(z,x,rho):
    E = np.exp(-x*z**2) ; A = 1+x ; R = 1-rho ; K = E*R*np.sqrt(A) ; D = rho+K
    term1 = -x*z/(A**2*D) ; term2 = z/(A*D)
    N3 = x*z*(E*R/(2*np.sqrt(A)) - E*R*np.sqrt(A)*z**2) ; D3 = A*D**2
    ret = -term1 + term2 - N3/D3
    ret*= rho**2
    return ret
def Vprime(x,rho):
    f = lambda z: aux(z,x=x,rho=rho)*z*gaussian(z=z)
    I = integrate.quad(f, -np.inf, np.inf)
    return I[0]
def integrand(t,alpha,rho,lamb):
    'Derivative of the free energy'
#     print("Computing the derivative of the free energy")
    a = Vprime(x=alpha*lamb*t,rho=rho)
    inp = update_v(x=alpha*lamb*t,rho=rho)
    b = t-update_u(lamb*inp)
    ret = a*b
    return ret
def obj(lamb,t,alpha,rho):
    inp = update_v(x=alpha*lamb*t,rho=rho)
    ret  = t-update_u(lamb*inp)
    return ret
def find_sol(t,alpha,rho,tol=10e-9,lamb_initial_guess = 0.1):
    eq = lambda lamb: obj(lamb,t=t,alpha=alpha,rho=rho)
    ret = fsolve(eq, lamb_initial_guess,xtol=tol)
    return ret

def deltaphi(val,lamb,alpha,rho):
    'Find the difference of free energy between trivial FP and non trivial FP at magnetization "val" '
    f = lambda t: integrand(t,alpha=alpha,rho=rho,lamb=lamb)
    I = integrate.quad(f, 0,val)
    ret = I[0]
    ret *= 0.5*(alpha*lamb)**2
#     print(f"deltaphi={ret}")
    return ret

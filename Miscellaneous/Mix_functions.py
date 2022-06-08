import numpy as np
import random
from sklearn.decomposition import SparsePCA
from scipy.linalg import eigh
### Principal Component Analysis ###
def get_instance(rho, dimension, nsamples, snr):
    u0 = np.sign(np.random.randn(nsamples))
    v0 = np.random.normal(0,1,dimension)
    elim = np.random.binomial(n=1,p=rho,size =dimension)
    v0 *= elim
    noise =  np.random.randn(dimension,nsamples)
    Y = np.sqrt(snr/dimension) * np.outer(v0,u0) + noise
    return u0, v0, Y
def oracle(X,rho,damp = 0.9,max_steps=200):
    'Finds the best lambda when fitting a data matrix X knowing that we should work with a sparsity rho'
    d,n = X.shape
    Y = X.T/(np.sqrt(n)) # Normalize such that both the spikes are of order one norm
    reg = 0.05 ; flag = True ; i=-1
    print("start the loop")
    while flag:
        i+=1 ; coef = 5
        transformer = SparsePCA(n_components=1,alpha=reg,random_state=0)
        transformer.fit(Y)
        shat = np.sum(transformer.components_!=0) ; strue = int(rho*d)
        obj =  shat - strue
        print(f"At iteration {i} we have regularization parameter = {reg}. \nThe estimated number of non-zero components is {shat} while the true number is {strue}")
        if i==max_steps:
            flag = False
        else:
            if obj < -1:
                reg = damp*reg + (1-damp)*reg/coef
            elif obj >1:
                reg = damp*reg + (1-damp)*reg*coef
            else:
                flag = False
    return reg
def SPCA(X,u0,v0,lam):
    d,n = X.shape
    Y = X.T/np.sqrt(n)  # Matrix on which we do SPCA
    transformer = SparsePCA(n_components=1,alpha=lam,random_state=0)
    transformer.fit(Y)
    Y_transformed = transformer.transform(Y)
    uhat = np.sign(Y_transformed[:,0]) ; vhat = transformer.components_
    retu = np.abs(np.dot(uhat,u0)/n)
    theta0 = v0!=0 ; theta_hat = vhat[0]!=0
    s = np.sum(v0!=0)
    # print(f"nonzerocomp estimated ={np.sum(transformer.components_!=0)} while true ={s}")
    retv = np.dot(theta0,theta_hat)/s
    mseu = 1-retu
    return mseu
### Diagonal Thresholding ###
def get_sample_cov(A):
    n,d = A.shape ; ret = np.zeros((d,d))
    for mu in range(n):
        ret += np.outer(A[mu,:],A[mu,:])
    return ret/n

def get_instance_normalized(rho, dimension, nsamples, snr):
    # Needed when the sparsity is sub-extensive to take into account normalization
    u0 = np.sign(np.random.randn(nsamples)) ; s = int(rho*dimension)
    aux1 = np.random.normal(0,1,s) ; v1 =aux1*np.sqrt(s)/(np.linalg.norm(aux1)) ; v2 = np.zeros(int(dimension - s))
    v0 = np.concatenate((v1,v2))
    noise =  np.random.randn(dimension,nsamples)
    Y = np.sqrt(snr/dimension) * np.outer(v0,u0) + noise
    return u0, v0, Y
def dtr(n,d,s0):
    alpha = n/d ;  rho = s0/d
    snr = 0.5/(rho*np.sqrt(alpha))  # Place yourself in the hard phase
    u0,v0,Ydata = get_instance_normalized(rho=rho, dimension=d, nsamples=n, snr=snr)
    ind_plant = np.sort(np.argpartition(np.abs(v0), -s0)[-s0:])
    X = Ydata.T
    muhat = 1/n*X.sum(axis=0)
    Xtild = X - np.outer(np.ones(n),muhat)
    Cov = get_sample_cov(Xtild) ; diag = np.diag(Cov) ;
    ind = np.sort(np.argpartition(np.abs(diag), -s0)[-s0:])
    Cov_tilda = np.zeros((d,d)) ; Cov_tilda1 = np.zeros((d,d)) ; Cov_tilda2 = np.zeros((d,d))
    Cov_tilda1[:,ind] = Cov[:,ind] ; Cov_tilda2[ind,:] = Cov[ind,:]
    Cov_tilda = Cov_tilda1*Cov_tilda2
    evals_large, evecs_large = eigh(Cov_tilda, eigvals=(d-1,d-1))
    vhat = evecs_large*np.sqrt(s0)/(np.linalg.norm(evecs_large))
    proj = vhat.T@Ydata  ;   uhat = np.sign(proj)
    ov = np.abs(np.dot(uhat,u0)/n)
    mseu = 1 - ov
    return ov

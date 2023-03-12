import numpy as np
from numpy import random




n = int(input("Enter the number of iterates: "))
print('n=',n)

q = float(input("Enter the value of q: "))
print('q=',q)

qstar = 1.0/(2-q)
print('qstar=',qstar)

beta = float(input("Enter the value of beta: "))
print('beta=',beta)


expr = int(input("Enter the number of experiments: "))
print('experiment=',expr)




for j in range(0, expr):

    fileq = open('q_estimates_for_q'+str(q)+'_beta'+str(beta)+'_with_'+str(expr)+'expr_n'+str(n)+'.txt','a')
    filebeta = open('beta_estimates_for_q'+str(q)+'_beta'+str(beta)+'_with_'+str(expr)+'expr_n'+str(n)+'.txt','a')

    u = random.rand(n)
    ln_qstar = (u**(1-qstar)-1)/(1-qstar)
    data = (1/beta)*(-qstar* ln_qstar)
    
    qval_n = np.zeros(expr)
    betaval_n = np.zeros(expr)
    
    x = np.zeros(n)
    nlist=np.arange(0,n,1)
    ii=0
    for i in nlist:
        x[i] = data[ii]
        ii=ii+1

    ss = 0.000231
    s  = 0.00023
    say = 0
    tol= 0.01
    theta = n*(sum(np.log(1+x/s)))**(-1)

    while tol > 0.0001 :
        say = say + 1
        theta = n*(sum(np.log(1+x/s)))**(-1)
        qval = 2 - 1/(1+ (1/theta))
        betaval = (1/ (s/theta)) / (2 - qval)
#        print(j,say,s,ss,tol,qval,betaval)
        if say>100000:
            s  = 200
            ss = 200
            tol= 0
            continue
        s   = s + 0.000001
        ss  = (1/n)*(1+n*(sum(np.log(1+x/s)))**(-1))*(sum(x/(1+x/s)))
        tol = np.abs(ss/s - 1)



#    print(j,say,s,ss,tol,qval,betaval)
    theta = n*(sum(np.log(1+x/s)))**(-1)

    qval_n[j] = 2 - 1/(1+ (1/theta))
    qval = 2 - 1/(1+ (1/theta))
#    print(j,qval)
    np.savetxt(fileq, np.c_[qval])

    betaval_n[j] = (1/ (s/theta)) / (2 - qval)
    betaval = (1/ (s/theta)) / (2 - qval)
#    print(j,betaval)
    np.savetxt(filebeta, np.c_[betaval])
    
    del(data)
    fileq.close()
    filebeta.close()
                


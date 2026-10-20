#Monte Carlo pricing
import numpy as np
def monte_carlo_put(S0,K,T,r,sigma,N):
    '''
    N: number of random sample path

    ------------
    return: list of floats and interval
        European put price via monte carlo simulation, std of sample path, se of estimator, confidence interval
    '''
    Z=np.random.randn(N)
    ST=S0*np.exp((r-1/2*sigma**2)*T)*np.exp((sigma*np.sqrt(T))*Z)
    payoff=np.exp(-r*T)*np.maximum(K-ST,0)

    # Monte Carlo price
    price=np.mean(payoff)

    #Sample standard deviation
    std=np.std(payoff,ddof=1)

    # Standard error of the estimator
    se=std/np.sqrt(N)

    # 95% confidence interval
    CI=(price-1.96*se,price+1.96*se) 

    return price,std,se,CI 


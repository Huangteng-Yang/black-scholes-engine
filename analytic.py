#Formulate the analytic BS formula of European put option 
import numpy as np
from scipy import stats 
def black_scholes_europut(S0,K,T,r,sigma):
    d1=(np.log(S0/K)+(r+1/2*sigma**2)*T)/(sigma*np.sqrt(T))
    d2=d1-sigma*np.sqrt(T)
    put_bs=K*np.exp(-r*T)*stats.norm.cdf(-d2)-S0*stats.norm.cdf(-d1)
    return put_bs

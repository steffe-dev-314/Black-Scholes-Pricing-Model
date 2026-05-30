
import numpy as np 
from scipy.stats import norm 
from scipy.optimize import newton , bisect
from dataclasses import dataclass
from typing import Callable
from root_solver import find_sign_change , do_bisection ,do_newton
# ------ Helpers ------ 

def _calc_PV(r : float ,x : float,tau : float)->float:
    return x*np.exp(-r*tau)

def _calc_intrinsic_value(S : float,K : float)->float:
    return max(0, S - K)

# no input validation since only used in classes where inputs already valid 
def _get_d1_d2(K :float , T : float , sigma :float , r : float, S : float , t:float)->tuple[float,float]:
    
    tau = T - t 

    tau_sigma = np.sqrt(tau)*sigma
    d1_fact_1 = tau_sigma**(-1) # tau = 0 --> div by 0 
    d1_fact_2 = np.log(S/K) + tau*(r + sigma**2/2) # S = 0 --> -inf 

    d1 = d1_fact_1 * d1_fact_2
    d2 = d1 - tau_sigma

    return d1 , d2   

def _get_d1_d2_prime(K :float , T : float , sigma :float , r : float, S : float , t:float)->tuple[float,float]:
    #d_1,2' = -log(S/K)*(1/sigma^2 sqrt(tau)) - r sqrt(tau)/sigma +- 1/2 sqrt(tau)
    tau = T - t
    s1 = -1*np.log(S/K)*(np.sqrt(tau)*sigma**2)**(-1) 
    s2 = -r * np.sqrt(tau)/sigma**2
    s3 = 1/2 * np.sqrt(tau)
    return s1 + s2 + s3 , s1 + s2 - s3 

# remove this ? 
def _get_normal_const(d1 : float,d2 : float) -> tuple[float,float]:
    d1 = float(d1)
    d2 = float(d2)
    n1 = float(norm.cdf(d1))
    n2 = float(norm.cdf(d2))
    return n1 , n2 


# ------ OptionData ------ 
@dataclass
class Option:
    K : float 
    T : float 
    is_call : bool 
    def __post_init__(self):
        self.K = float(self.K)
        self.T = float(self.T)
        if self.K <= 0:
            raise ValueError("Detected non-positive strike!")
        if self.T < 0: # handle T = 0 in calcs ! 
            raise ValueError("Detected negative maturity!")

# ------ MarketData ------ 
@dataclass
class Market:

    S : float 
    sigma : float 
    r : float 

    def __post_init__(self):
        self.S = float(self.S)
        self.sigma = float(self.sigma)
        self.r = float(self.r)

        if self.S < 0: # = 0 special case    
            raise ValueError("Detected negative underlying!")
        if self.sigma < 0: # = 0 special case 
            raise ValueError("Detected negative volatility!")
        if self.r < 0: # = 0 special case
            raise ValueError("Detected negative risk-free rate!")

# ----- Pricing Engine ------        

class Pricing:
            
    # ------ Pricing ------

    def _check_edge_cases(self,option: Option, market:Market , t :float)->str:

        tau = option.T - t

        if option.T == 0 or tau == 0: # implies t = tau = 0 due to above check 
            return 'AT_MATURITY' # self._get_intrinsic_value(option=option,market=market)
        if option.K == 0:
            return 'ZERO_STRIKE' #market.S 
        if market.S == 0:
            return 'ZERO_UNDERLYING' # 0
        if market.sigma == 0: # deterministic case --> price = S_0 - exp(-rtau)K
            return 'DETERMINISTIC' #max(0,market.S - np.exp(-market.r*tau)*option.K)
        else:
            return 'USE_FORMULA'

    def get_option_price(self, option :Option, market:Market , t:float)->float:

        tau = option.T - t

        #edge_cases = ['AT_MATURITY' , 'ZERO_STRIKE' , 'ZERO_UNDERLYING' , 'DETERMINISTIC']

        current_case = self._check_edge_cases(option,market,t)

        if current_case == 'AT_MATURITY':
            return _calc_intrinsic_value(K=option.K,S=market.S)
        if current_case == 'ZERO_STRIKE':
            return market.S 
        if current_case == 'ZERO_UNDERLYING':
            return 0
        if current_case == 'DETERMINISTIC':
            return max(0,market.S - np.exp(-market.r*tau)*option.K)

        d1 , d2 = _get_d1_d2(K = option.K , T = option.T , sigma = market.sigma , S= market.S  ,r = market.r, t = t)
        n1 , n2 = _get_normal_const(d1 = d1,d2 = d2)
        PV_K = _calc_PV(r = market.r, x = option.K,tau = tau)

        sum1 = n1*market.S
        sum2 = n2*PV_K
        price_call = sum1 - sum2
        if option.is_call: 
            return price_call 
        else: 
            return price_call - market.S + option.K * np.exp(-market.r * tau) # put price from put call partiry
# ----- analysis class ----- 

class Analysis:

    # ----- 1. order greeks ----- 
    def get_delta(self, option:Option , market:Market , model:Pricing ,t:float)->float:
        # del price / del underlying 
        d1 , d2 = _get_d1_d2(K = option.K , T = option.T , sigma = market.sigma , S= market.S  ,r = market.r, t = t)
        n1 , _ = _get_normal_const(d1 = d1,d2 = d2)  
        if option.is_call:
            return n1 
        else:
            return n1 - 1

    
    def get_vega(self , option : Option , market : Market   , model :Pricing , t:float) -> float:
        # del price / del vol 
        tau = option.T -t 
        d1 , _ = _get_d1_d2(K = option.K , T = option.T , sigma = market.sigma , S= market.S  ,r = market.r, t = t)
        n1_prime = norm.pdf(d1)
        S = market.S 
        return S * np.sqrt(tau) * n1_prime # no difference put or call 

    def get_theta(self, option : Option , market :Market ,model:Pricing, t: float)-> float:
        # - del price / del tau 
        S = market.S 
        sigma = market.sigma 
        tau = option.T - t 
        r = market.r 
        _ , d2 = _get_d1_d2(K = option.K , T = option.T , sigma = market.sigma , S= market.S  ,r = market.r, t = t)
        K = option.K 

        n2 =  norm.cdf(d2)
        n2_prime = norm.pdf(d2)

        fact_1 = S*sigma/(2*np.sqrt(tau)) 
        fact_2 = r 
        theta_call = -K*np.exp(-r*tau)*(fact_1 * n2_prime + fact_2 * n2)
        if option.is_call:
            return theta_call 
        else: 
            return theta_call + r*K*np.exp(-r*tau)
    
    def get_rho(self, option :Option , market :Market , model : Pricing , t:float) -> float:
        # del price / del r 
        tau = option.T - t 
        K = option.K 
        r = market.r 
        d1 , d2 = _get_d1_d2(K = option.K , T = option.T , sigma = market.sigma , S= market.S  ,r = market.r, t = t)
        _ , n2 = _get_normal_const(d1 = d1,d2 = d2)     
        rho_call = tau*K*np.exp(-r*tau)*n2 
        if option.is_call: 
            return rho_call 
        else: 
            return rho_call - tau*np.exp(-r*tau) 

    def get_lambda(self, option : Option , market : Market , model : Pricing , t : float)-> float:
        # underlying/price * Delta 
        delta= self.get_delta(option, market , model , t)
        return delta * (market.S / option.K)
    
    # ----- 2. order greeks ----- 



# ------ IV Solver ------ 

# this should later accept Pricing so we have access to its pricing function directly
class IVSolver:
    def __init__(self, option : Option , market:Market):
        self.r = market.r
        self.S = market.S 
        self.K = option.K 
        self.T = option.T

    # temporary : re-implement BS formula here again (redudant, maybe rework)
    def get_f(self,sigma , t , P_market): 
        K = self.K 
        T = self.T 
        tau = T - t 
        r = self. r 
        S= self.S 

        d1 , d2 = _get_d1_d2(sigma = sigma ,K = K , T  = T ,r = r , S = S , t = t)
        n1 , n2 = norm.cdf(float(d1)) , norm.cdf(float(d2))

        BS_price = S*n1 - K*np.exp(-r*tau)*n2 
        return BS_price - P_market

    def get_f_prime(self,sigma , t ,P_market): 
        K = self.K 
        T = self.T 
        tau = T - t 
        r = self. r 
        S= self.S 

        d1 , d2 = _get_d1_d2(sigma = sigma ,K = K , T  = T ,r = r , S = S , t = t)
        d1_prime , d2_prime = _get_d1_d2_prime(sigma = sigma ,K = K , T  = T ,r = r , S = S , t = t)
        n1 , n2 = norm.cdf(float(d1)) , norm.cdf(float(d2))
        n1_prime , n2_prime = norm.pdf(float(d1)) , norm.pdf(float(d2))

        BS_price_prime = S*n1_prime*d1_prime - K*np.exp(-r*tau)*n2_prime*d2_prime
        
        return BS_price_prime

    def get_IV(self,t ,method , P_market ,x0 = None):
        if x0 is None:
            x0 = 0.1
        if x0 == 0:
            raise ValueError(f'Got zero volatility which is not permitted!')

        if method == 'scipy_newton': 
            x_root = newton(func = self.get_f, fprime = self.get_f_prime , x0 = x0, args = (t,P_market) , maxiter = 100 ,tol = 1e-6)
        elif method == 'own_newton':
            x_root = do_newton(x0 = x0 , f = self.get_f , fprime = self.get_f_prime , maxiter = 100 , tol = 1e-6 , args = (t,P_market))

        return x_root 
        



if __name__ == '__main__':
    import time 

    start_scipy = time.perf_counter()
    # standard setting : S = K = 100 , T = 1 , r = 0.02 , sigma = 0.2 , t = 0, should give P = 8.9 
    market = Market(S= 100, sigma = 0.2, r = 0.02)
    #option = Option(K = 100 , T = 1 , option_type='call' , european=True)
    option = Option(K = 100 , T = 1 , is_call = True)

    model = Pricing()
    analysis = Analysis()
    solver = IVSolver(option=option, market = market)

    #print(model.get_option_price(option=option , market=market , t= 0))
    #print(analysis.get_delta(option , market , model , 0))
    iv_scipy = solver.get_IV(t = 0, method = 'scipy_newton' , P_market = 8.9, x0= 1)
    end_scipy = time.perf_counter()
    print(f'iv_scipy : {iv_scipy}')
    print(f'time :  {(end_scipy - start_scipy)*1e6:.2f} ms')
    #print(solver.get_f(sigma = 0.2, t = 0 , P_market= 8.9))
    #print(solver.get_f_prime(sigma = 0.2, t = 0 , P_market= 8.9))
    

    start_own = time.perf_counter()
    iv_own = solver.get_IV(t = 0 , method = 'own_newton' , P_market = 8.9 , x0 = 1)
    end_own = time.perf_counter()
    print(f'iv_own = {iv_own}')
    print(f'time :  {(end_own - start_own)*1e6:.2f} ms')

    a0 = 0.1
    maxiter = 1000
    s = 10/maxiter*max(a0,1) # a0 should somehow be bounded judging from a bit of testing 

    a,b = find_sign_change(f = solver.get_f , args = (0,8.9) , a0 = a0, maxiter = maxiter , stepsize= s) 
    #print(np.round(a,2),np.round(b,2))
    #print(solver.get_f(a,t=0 , P_market = 8.9) , solver.get_f(b,t=0 , P_market = 8.9))
    print(do_bisection(f = solver.get_f , args = (0,8.9) ,a0 = a0, maxiter = maxiter ,root_tol = 1e-6 , stepsize = s))


import numpy as np 
from scipy.stats import norm 
from scipy.optimize import minimize_scalar

# ------ Helpers ------ 

def _validate_float(x):
    try:
        float(x)
        return float(x)
    except ValueError:
        raise ValueError(f"Expected float or int not {type(x)}")

def _calc_PV(r : float ,x : float,tau : float)->float:
    return x*np.exp(-r*tau)

def _get_normal_const(d1 : float,d2 : float) -> tuple[float,float]:
    d1 = _validate_float(d1)
    d2 = _validate_float(d2)
    n1 = float(norm.cdf(d1))
    n2 = float(norm.cdf(d2))
    return n1 , n2 

# ------ OptionData ------ 

class Option:
    def __init__(self, K : float , T : float , option_type : str , european : bool):
        K = _validate_float(K)
        T = _validate_float(T)
        if K <= 0:
            raise ValueError("Detected non-positive strike!")
        if T < 0: # handle T = 0 in calcs ! 
            raise ValueError("Detected negative maturity!")
        self.K = K 
        self.T = T 

# ------ MarketData ------ 

class Market:
    def __init__(self , S : float , sigma : float , r : float):
        S = _validate_float(S)
        sigma = _validate_float(sigma)
        r = _validate_float(r)

        if S < 0: # = 0 special case    
            raise ValueError("Detected negative underlying!")
        if sigma < 0: # = 0 special case 
            raise ValueError("Detected negative volatility!")
        if r < 0: # = 0 special case
            raise ValueError("Detected negative risk-free rate!")
        
        self.S = S 
        self.sigma = sigma
        self.r = r 

# ----- Pricing Engine ------        

class Pricing:

    # ------ Helper Functions for const ------

    def _get_d1_d2(self,option :Option, market:Market , t:float)->tuple[float,float]:
        
        tau = option.T - t # tau here implicit from option.T and t 

        if tau < 0:
            raise ValueError(f"Negative time to maturity detected, please make sure t < T, not {t} > {option.T}")


        sigma = market.sigma
        S = market.S
        r = market.r
        K =option.K
        tau_sigma = np.sqrt(tau)*sigma
        d1_fact_1 = tau_sigma**(-1) # tau = 0 --> div by 0 
        d1_fact_2 = np.log(S/K) + tau*(r + sigma**2/2) # S = 0 --> -inf 
        d1 = d1_fact_1 * d1_fact_2
        d2 = d1 - tau_sigma

        return d1 , d2        
    
    def _get_intrinsic_value(self,option:Option,market:Market):
        return max(0, market.S - option.K)

    # ------ Pricing ------

    def get_option_price(self, option :Option, market:Market , t:float)->float:

        tau = option.T - t

        if option.T == 0 or tau == 0: # implies t = tau = 0 due to above check 
            return self._get_intrinsic_value(option=option,market=market)
        if option.K == 0:
            return market.S 
        if market.S == 0:
            return 0 
        if market.sigma == 0: # deterministic case --> price = S_0 - exp(-rtau)K
            return max(0,market.S - np.exp(-market.r*tau)*option.K)

        d1 , d2 = self._get_d1_d2(option , market , t)
        n1 , n2 = _get_normal_const(d1 = d1,d2 = d2)
        PV_K = _calc_PV(r = market.r, x = option.K,tau = tau)

        sum1 = n1*market.S
        sum2 = n2*PV_K

        return sum1 - sum2

# ----- analysis class ----- 

class Analysis:

    # ----- 1. order greeks ----- 
    def get_delta(self, option:Option , market:Market , model:Pricing ,t:float)->float:
        # del price / del underlying 
        d1 , d2 = model._get_d1_d2(option , market , t)
        n1 , _ = _get_normal_const(d1 = d1,d2 = d2)  
        return n1 
    
    def get_vega(self , option : Option , market : Market   , model :Pricing , t:float) -> float:
        # del price / del vol 
        pass 

    def get_theta(self, option : Option , market :Market ,model:Pricing, t: float)-> float:
        # - del price / del tau 
        pass 

    def get_rho(self, option :Option , market :Market , model : Pricing , t:float) -> float:
        # del price / del interest 
        pass 

    def get_lambda(self, option : Option , market : Market , model : Pricing , t : float)-> float:
        # underlying/price * Delta 
        delta= self.get_delta(option, market , model , t)
        return delta * (market.S / option.K)
    
    # ----- 2. order greeks ----- 



    # ------ IV Solver ------ 


if __name__ == '__main__':
    # standard setting : S = K = 100 , T = 1 , r = 0.02 , sigma = 0.2 , t = 0, should give P = 8.9 
    market = Market(S= 100, sigma = 0.2, r = 0.02)
    option = Option(K = 100 , T = 1 , option_type='call' , european=True)
    model = Pricing()
    analysis = Analysis()


    print(model.get_option_price(option=option , market=market , t= 0))
    print(analysis.get_delta(option , market , model , 0))



import numpy as np 
from scipy.stats import norm 
from scipy.optimize import minimize_scalar
from dataclasses import dataclass

# ------ Helpers ------ 

def _calc_PV(r : float ,x : float,tau : float)->float:
    return x*np.exp(-r*tau)

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

    # ------ Helper Functions for const ------

    def _get_d1_d2(self,option :Option, market:Market , t:float)->tuple[float,float]:
        
        tau = option.T - t # tau here implicit from option.T and t 
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
            return self._get_intrinsic_value(option=option,market=market)
        if current_case == 'ZERO_STRIKE':
            return market.S 
        if current_case == 'ZERO_UNDERLYING':
            return 0
        if current_case == 'DETERMINISTIC':
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
        tau = option.T -t 
        d1 , _ = model._get_d1_d2(option , market , t)
        n1_prime = norm.pdf(d1)
        S = market.S 
        return S * np.sqrt(tau) * n1_prime

    def get_theta(self, option : Option , market :Market ,model:Pricing, t: float)-> float:
        # - del price / del tau 
        S = market.S 
        sigma = market.sigma 
        tau = option.T - t 
        r = market.r 
        _ , d2 = model._get_d1_d2(option , market , t)
        K = option.K 

        n2 =  norm.cdf(d2)
        n2_prime = norm.pdf(d2)

        fact_1 = S*sigma/(2*np.sqrt(tau)) 
        fact_2 = r 

        return -K*np.exp(-r*tau)*(fact_1 * n2_prime + fact_2 * n2)
    
    def get_rho(self, option :Option , market :Market , model : Pricing , t:float) -> float:
        # del price / del interest 
        return 0

    def get_lambda(self, option : Option , market : Market , model : Pricing , t : float)-> float:
        # underlying/price * Delta 
        delta= self.get_delta(option, market , model , t)
        return delta * (market.S / option.K)
    
    # ----- 2. order greeks ----- 



    # ------ IV Solver ------ 


if __name__ == '__main__':
    # standard setting : S = K = 100 , T = 1 , r = 0.02 , sigma = 0.2 , t = 0, should give P = 8.9 
    market = Market(S= 100, sigma = 0.2, r = 0.02)
    #option = Option(K = 100 , T = 1 , option_type='call' , european=True)
    option = Option(K = 100 , T = 1)

    model = Pricing()
    analysis = Analysis()


    print(model.get_option_price(option=option , market=market , t= 0))
    print(analysis.get_delta(option , market , model , 0))


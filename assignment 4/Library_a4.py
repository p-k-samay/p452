import numpy as np
import scipy as scipy
import math as math
from scipy.optimize import root

class rng():
    def __init__(self,seed, a = 1103515245, c = 12345 ,m = 32768):
        self.term = seed
        self.a = a
        self.c = c
        self.m = m
        
    def gen(self):
        self.term = (((self.a * self.term) + self.c) % self.m)
        return self.term/self.m
    
    def genlist(self,length):
        RNs = []
        for i in range(length):
            self.term = (((self.a * self.term) + self.c) % self.m)
            RNs.append(self.term / self.m)
        return RNs 
    
def monte_carlo_integrate(f: float,a: float,b: float,N: int,seed: int,multiplier=1103515245,m=32768,c=12345):
    p=rng(seed,m=m,c=c,a=multiplier)
    F=0
    for i in range(N):
        k=p.gen()
        k=((b-a)*(k))+a
        F+=((b-a)*f(k))/N   
    return F   


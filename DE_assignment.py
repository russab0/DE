#Author: Ruslan Sabirov, BS17-02
#Variant: 23

"""
Task:
y' = f(x, y)
y(x0) = y0
x ∈ [x0, X]
f(x, y) = y^2 * e^x + 2 * y
"""

import math
import numpy as np
import plotly
import plotly.graph_objs as go



class Numeric_methods():
    h = n = None # h - step, n - number of grid steps
    EPS = 3 * 10 ** (-3) # epsilon
    x_discont = 1.38847 # when denominator of exact solution is equal to 0
    e = math.e # e constant
    
    
    # Returns True if point x lies around point of discontinuity and False otherwise
    def lies_around_discont(self, x):
        return self.x_discont - self.EPS < x < self.x_discont + self.EPS
    
    
    # Given function
    def f(self, x, y):
        e = self.e
        return (y ** 2) * (e ** x) + 2 * y    
    
    
    # Exact solution of given function
    def exact(self, x):
        e = self.e
        solution = lambda x: -3 * e ** (2 * x) / (e ** (3 * x) - 64.4199)
        if type(x) != type(np.arange(0, 1)):
            return solution(x)
        
        y = [0] * self.n
        for i in range(self.n):
            if self.lies_around_discont(x[i]):
                y[i] = None
            else:
                y[i] = solution(x[i])
        return y
    
    
    def __init__(self, x0, y0, X, n):
        self.n = n
        self.h = (X - x0) / n
        EPS = self.EPS
        x_discont = self.x_discont        
        
        # Calculating x-es and y-s        
        x = np.arange(x0, X, self.h)
        es_y = self.euler_standart(x, x0, y0, X)
        ei_y = self.euler_improved(x, x0, y0, X)
        rk_y = self.runge_kutta(x, x0, y0, X)
        ex_y = self.exact(x)
        
        # Creating traces
        lm = "lines+markers"
        es_trace = go.Scatter(x = x, y = es_y, name = "Euler", mode = lm)
        ei_trace = go.Scatter(x = x, y = ei_y, name = "Improved Euler", mode = lm)        
        rk_trace = go.Scatter(x = x, y = rk_y, name = "Runge Kutta", mode = lm)
        ex_trace = go.Scatter(x = x, y = ex_y, name = "Exact Solution", mode = lm)        
        
        # Drawing a graph
        data = [es_trace, ei_trace, rk_trace, ex_trace]
        plotly.offline.plot(data, filename="solutions.html")        
        
    
    # Euler method
    def euler_standart(self, x, x0, y0, xf): 
        h = self.h
        f = self.f
        y = [0] * len(x)
        y[0] = y0
        
        for i in range(1, len(x)):
            if self.lies_around_discont(x[i]):
                y[i] = None
                continue
            if len(y) > 1 and y[i - 1] is None:
                y[i] = self.exact(x[i])
                continue                
            
            y[i] = y[i - 1] + h * f(x[i - 1], y[i - 1])
        return y    
    
    
    # Improved Euler method
    def euler_improved(self, x, x0, y0, xf):
        h = self.h
        f = self.f
        y = [0] * len(x)
        y[0] = y0
        
        for i in range(1, len(x)):
            if self.lies_around_discont(x[i]):
                y[i] = None
                continue
            if len(y) > 1 and y[i - 1] is None:
                y[i] = self.exact(x[i])
                continue   
            
            delta_y = h * f(x[i - 1] + h / 2, y[i - 1] + h / 2 * f(x[i - 1], y[i - 1]))
            y[i] = y[i - 1] + delta_y
        return y    
    
    
    # Runge-Kutta method
    def runge_kutta(self, x, x0, y0, xf):
        h = self.h
        f = self.f
        y = [0] * len(x)
        y[0] = y0
        
        for i in range(1, len(x)):
            if self.lies_around_discont(x[i]):
                y[i] = None
                continue
            if len(y) > 1 and y[i - 1] is None:
                y[i] = self.exact(x[i])
                continue               
            
            x_prev = x[i - 1]
            y_prev = y[i - 1]
            k1 = f(x_prev, y_prev)
            k2 = f(x_prev + h / 2, y_prev + h * k1 / 2)
            k3 = f(x_prev + h / 2, y_prev + h * k2 / 2)
            k4 = f(x_prev + h, y_prev + h * k3)
            
            delta_y = h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            
            y[i] = y[i - 1] + delta_y
        return y    
    
    
Numeric_methods(x0 = 1, y0 = 0.5, X = 1.7769400000000002, n = 600)
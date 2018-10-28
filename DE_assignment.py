
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
import matplotlib.pyplot as plt
import plotly
import plotly.graph_objs as go


def func(x, y):
    e = math.e
    return (y ** 2) * (e ** x) - 2 * y #+++++++++ change to plus

class Numeric_methods():
    x0 = y0 = xf = f = None
    
    
    def __init__(self, x0, y0, X, f, h = 1/100):
        self.x0 = x0
        self.y0 = y0
        self.xf = f
        self.f = f
        self.h = h
        
        x = np.arange(x0, X, h)
        e_y = self.euler(x)
        ei_y = self.euler_improved(x)
        rk_y = self.runge_kutta(x)
        print(x, len(x))
        # Create a trace
        e_trace = go.Scatter(
            x = x,
            y = e_y,
            name = "Euler",
            mode = "lines+markers"
        )
        ei_trace = go.Scatter(
            x = x,
            y = ei_y,
            name = "Improved Euler",
            mode = "lines+markers"
        )        
        rk_trace = go.Scatter(
            x = x,
            y = rk_y,
            name = "Runge Kutta",
            mode = "lines+markers"
        )                
        
        data = [e_trace, ei_trace]#, rk_trace]
        
        plotly.offline.plot(data, filename="basic-line.html")        
        
    
    def euler( x0, y, h, x ): 
        temp = -0
      
        # Iterating till the point at which we 
        # need approximation 
        while x0 < x: 
            temp = y 
            y = y + h * func(x0, y) 
            x0 = x0 + h         
    
    def euler(self, x): 
        y0, h, f = self.y0, self.h, self.f
        y = [0] * len(x)
        y[0] = y0
        
        for i in range(1, len(x)):
            y[i] = y[i - 1] + h * f(x[i - 1], y[i - 1])
        return y    
    
    
    def euler_improved(self, x):
        y0, h, f = self.y0, self.h, self.f
        y = [0] * len(x)
        y[0] = y0
        
        for i in range(1, len(x)):
            delta_y = h * f(x[i - 1] + h / 2, y[i - 1] + h / 2 * f(x[i - 1], y[i - 1]))
            y[i] = y[i - 1] + delta_y
        return y    
    
    
    def runge_kutta(self, x):
        y0, h, f = self.y0, self.h, self.f
        y = [0] * len(x)
        y[0] = y0
        
        for i in range(1, len(x)):
            x_prev = x[i - 1]
            y_prev = y[i - 1]
            k1 = f(x_prev, y_prev)
            k2 = f(x_prev + h / 2, y_prev + h * k1 / 2)
            k3 = f(x_prev + h / 2, y_prev + h * k2 / 2)
            k4 = f(x_prev + h, y_prev + h * k3)
            
            delta_y = h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            
            y[i] = y[i - 1] + delta_y
        return y    
    
    
#TODO change X to 7 
X = 7
Numeric_methods(x0 = 1, y0 = 0.5, X = X, f = func, h = 1/50)
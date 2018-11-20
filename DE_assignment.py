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
        
        y = [0] * len(x)
        for i in range(len(x)):
            if self.lies_around_discont(x[i]):
                y[i] = None
            else:
                y[i] = solution(x[i])
        return y
    
    
    def __init__(self, x0, y0, X, n):
        h = (X - x0) / n
        self.EPS = h / 1.1
        EPS = self.EPS
        x_discont = self.x_discont 
        
        # Calculating x-es and y-s        
        x = np.arange(x0, X, h)
        es_y = self.euler_standart(x, x0, y0, X, h)
        ei_y = self.euler_improved(x, x0, y0, X, h)
        rk_y = self.runge_kutta(x, x0, y0, X, h)
        ex_y = self.exact(x)
        
        # Creating traces
        lm = "lines+markers"
        es_trace = go.Scatter(x = x, y = es_y, name = "Euler", mode = lm)
        ei_trace = go.Scatter(x = x, y = ei_y, name = "Improved Euler", mode = lm)        
        rk_trace = go.Scatter(x = x, y = rk_y, name = "Runge Kutta", mode = lm)
        ex_trace = go.Scatter(x = x, y = ex_y, name = "Exact Solution", mode = lm)        
        
        # Specifying layout
        layout = dict(title = "Graphs of solutions | X = {}, n = {}".format(X, n),
                      xaxis = dict(title = "x"),
                      yaxis = dict(title = "y"))        
        
        # Drawing a graph
        data = [es_trace, ei_trace, rk_trace, ex_trace]
        plotly.offline.plot(dict(data=data, layout=layout), filename="solutions.html")        
        
        self.truncation_error(x0, y0, X)
        
    
    # Plots graph of truncation errors over steps number
    def truncation_error(self, x0, y0, xf):
        #change!!!
        x0 = 2
        y0 = self.exact(x0)
        
        
        start_step = 20 # could be changed
        end_step = 200 # could be changed
        steps = np.arange(start_step, end_step, 1)
        es_error = [0] * len(steps)
        ei_error = [0] * len(steps)
        rk_error = [0] * len(steps)

        for i in range(len(steps)):
            n = steps[i]
            h = (xf - x0) / n
            self.EPS = h / 2
            
            x = np.arange(x0, xf, h)
            es_y = self.euler_standart(x, x0, y0, xf, h)
            ei_y = self.euler_improved(x, x0, y0, xf, h)
            rk_y = self.runge_kutta(x, x0, y0, xf, h)
            ex_y = self.exact(x)
            
            es_error[i] = max([abs(ex_y[i] - es_y[i]) for i in range(n) if es_y[i] and ex_y[i]])
            ei_error[i] = max([abs(ex_y[i] - ei_y[i]) for i in range(n) if ei_y[i] and ex_y[i]])
            rk_error[i] = max([abs(ex_y[i] - rk_y[i]) for i in range(n) if rk_y[i] and ex_y[i]])
            
        
        lm = "lines+markers"
        es_err_trace = go.Scatter(x = steps, y = es_error, name = "Euler", mode = lm)
        ei_err_trace = go.Scatter(x = steps, y = ei_error, name = "Improved Euler", mode = lm)  
        rk_err_trace = go.Scatter(x = steps, y = rk_error, name = "Runge Kutta", mode = lm)   
        
        layout = dict(title = "Graph of global truncation error depending on number of steps",
                      xaxis = dict(title = "n — number of steps"),
                      yaxis = dict(title = "Global error — Maximum of local truncation errors"))
        
        data = [es_err_trace, ei_err_trace, rk_err_trace]
        plotly.offline.plot(dict(data=data, layout=layout), filename="trunc_errors.html")          
            
    
    # Euler method
    def euler_standart(self, x, x0, y0, xf, h): 
        f = self.f
        y = [0] * len(x)
        y[0] = y0
        
            
        for i in range(1, len(x)):
            if self.lies_around_discont(x[i]): # current point lies around discontinuity
                y[i] = None
                continue
            if len(y) > 1 and y[i - 1] is None: # previous point lies around discontinuity
                y[i] = self.exact(x[i])
                continue                
            
            y[i] = y[i - 1] + h * f(x[i - 1], y[i - 1])
        return y    
    
    
    # Improved Euler method
    def euler_improved(self, x, x0, y0, xf, h):
        f = self.f
        y = [0] * len(x)
        y[0] = y0
        
        for i in range(1, len(x)):
            if self.lies_around_discont(x[i]): # current point lies around discontinuity
                y[i] = None
                continue
            if len(y) > 1 and y[i - 1] is None: # previous point lies around discontinuity
                y[i] = self.exact(x[i])
                continue       
            
            delta_y = h * f(x[i - 1] + h / 2, y[i - 1] + h / 2 * f(x[i - 1], y[i - 1])) # augmentation
            y[i] = y[i - 1] + delta_y
        return y    
    
    
    # Runge-Kutta method
    def runge_kutta(self, x, x0, y0, xf, h):
        f = self.f
        y = [0] * len(x)
        y[0] = y0
        
        for i in range(1, len(x)):
            if self.lies_around_discont(x[i]): # current point lies around discontinuity
                y[i] = None
                continue
            if len(y) > 1 and y[i - 1] is None: # previous point lies around discontinuity
                y[i] = self.exact(x[i])
                continue                  
            
            # Assigning coordinates of previous point to variables
            x_prev = x[i - 1]
            y_prev = y[i - 1]
            k1 = f(x_prev, y_prev)
            k2 = f(x_prev + h / 2, y_prev + h * k1 / 2)
            k3 = f(x_prev + h / 2, y_prev + h * k2 / 2)
            k4 = f(x_prev + h, y_prev + h * k3)
            
            delta_y = h / 6 * (k1 + 2 * k2 + 2 * k3 + k4) # augmentation
            
            y[i] = y[i - 1] + delta_y
        return y    
    
# x0, y0, X, n - could be changed    
Numeric_methods(x0 = 1, y0 = 0.5, X = 7, n = 500)
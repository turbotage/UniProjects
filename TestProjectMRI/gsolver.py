# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 16:34:50 2021

@author: TurboTage
"""

import torch

# SETUP MODEL    

def get_model(model_expr):
    # Replace math expr in order of pytorch docs
    model_expr = model_expr.replace('abs', 'torch.abs')
    model_expr = model_expr.replace('acos', 'torch.acos')
    model_expr = model_expr.replace('acosh', 'torch.acosh')
    model_expr = model_expr.replace('asin', 'torch.sin')
    model_expr = model_expr.replace('asinh', 'torch.asinh')
    model_expr = model_expr.replace('atan', 'torch.atan')
    model_expr = model_expr.replace('atanh', 'torch.atanh')
    model_expr = model_expr.replace('atan2', 'torch.atan2')
    #model_expr = model_expr.replace('conj', 'torch.conj')
    model_expr = model_expr.replace('cos', 'torch.cos')
    model_expr = model_expr.replace('cosh', 'torch.cosh')
    model_expr = model_expr.replace('exp', 'torch.exp')
    model_expr = model_expr.replace('log', 'torch.log')
    model_expr = model_expr.replace('pow', 'torch.pow')
    model_expr = model_expr.replace('rsqrt', 'torch.rsqrt')
    model_expr = model_expr.replace('sin', 'torch.sin')
    model_expr = model_expr.replace('sinh', 'torch.sinh')
    model_expr = model_expr.replace('sqrt', 'torch.sqrt')
    model_expr = model_expr.replace('tan', 'torch.tan')
    model_expr = model_expr.replace('tanh', 'torch.tanh')
    
    # Paramaters
    for i in range(50):
        model_expr = model_expr.replace('P' + str(i), 'params[:][' + str(i) + ']')
    
    # Variables
    for i in range(50):
        model_expr = model_expr.replace('X' + str(i), 'dependent[' + str(i) + ']')
    
    #print(model_expr)        
    return lambda dependent, params: eval(model_expr)


def _get_func(model,dependent):
    return lambda params: model(dependent, params)

# END OF MODEL


# START OF SOLVER

def _get_jacobians(f, params):
    JT = torch.autograd.functional.jacobian(f, params).transpose(0,2)
    J = JT.transpose(1,2)
    return J, JT

def _get_diff(func, params, data):
    diff = (data - func(params)).transpose(0,1).reshape(data.size()[1], 1, data.size()[0]).transpose(1,2)
    return diff

def _get_step(J, JT, JT_diff, guess_order):
    # AX=B
    JTJ = torch.bmm(JT,J) # (JTJ)
    
    ill_behaved = (torch.logical_or(torch.det(JTJ) < 1, torch.det(JTJ) > 1e30)).reshape(JTJ.size()[0],1,1)
    
    eye = torch.eye(JTJ.size()[1]).repeat(JTJ.size()[0],1).reshape(JTJ.size()[0],JTJ.size()[1],JTJ.size()[1]).cuda()
    
    perturb = eye * guess_order * ill_behaved * 1e10
    
    JTJ = JTJ *((ill_behaved * 1e-21) + torch.logical_not(ill_behaved)) + perturb
    
    JT_diff = JT_diff * ((ill_behaved * 1e-21) + torch.logical_not(ill_behaved))
    
    
    # Solve via Cholskey decomposition, in theory better than LU for this problem
    u = torch.cholesky(JTJ)
    param_step = torch.cholesky_solve(JT_diff,u).transpose(0,2)[0]
    
    return param_step


def solve(model, dependent, data, guess, guess_order, threshold = 0.0001, max_iter=10000):
    
    func = _get_func(model,dependent)
    f = lambda p: func(p).sum(dim=1)
    
    n_params = guess.size()[0]
    n_pixels = guess.size()[1]
    
    non_convergence_list = torch.arange(0,n_pixels).cuda()
    
    
    # Initial Step
    J, JT = _get_jacobians(f, guess)
    diff = _get_diff(func, guess, data)
    JT_diff = torch.bmm(JT, diff)
    
    old_S = torch.zeros(1,n_pixels).cuda()
    S = (diff * diff).sum(dim=1).transpose(0,1)
    
    #print(S.shape)
    #print(old_S.shape)
    
    old_step = torch.zeros(n_params, n_pixels).cuda()
    step = _get_step(J, JT, JT_diff, guess_order)
    
    #print((torch.abs((old_S - S) / S) < threshold).shape)
    #print((torch.abs((old_step - step) / guess) < threshold * 10).shape)
    
    converges = ((torch.abs((old_step - step) / guess) < threshold * 10).sum(dim=0) + (torch.abs((old_S - S) / S) < threshold))[0,:]

    params = guess
    #print("initial step done")
    
    iteration = 0
    while (converges.sum() != guess.size()[1] * (guess.size()[0] + 1)) and iteration < max_iter:
        iteration += 1
        
        guess = guess + step
        
        J, JT = _get_jacobians(f, guess)
        diff = _get_diff(func,guess,data)
        JT_diff = torch.bmm(JT,diff)
        
        
        old_step = step
        step = _get_step(J, JT, JT_diff, guess_order)
        
        old_S = S
        S = (diff * diff).sum(dim=1).transpose(0,1)
        
        converges = ((torch.abs((old_step - step) / guess) < threshold * 10).sum(dim=0) + (torch.abs((old_S - S) / S) < threshold))[0,:]
        
        
        if iteration % 10 == 0:
            not_converging = converges < 3
            converging = torch.logical_not(not_converging)
            
            # copy converging pixels
            params[:,non_convergence_list[converging]] = guess[:,converging]
            
            
            # rebuild solver to solve with non-converging pixels
            
            dependent = dependent[:,:,not_converging]
            data = data[:,not_converging]
            guess = guess[:,not_converging]
            step = step[:,not_converging]
            S = S[:,not_converging]
            converges = converges[not_converging]
            
            func = _get_func(model,dependent)
            f = lambda p: func(p).sum(dim=1)
            
            # rebuild non convergence list
            non_convergence_list = non_convergence_list[not_converging]
            
            
            convergence_percentage = float(n_pixels - (converges < 3).sum().cpu().numpy()) / float(n_pixels)
            print((1-convergence_percentage) * n_pixels)
            #print(step)
        
        
    params[:,non_convergence_list[converges > 2]] = guess[:,converges > 2]
    
    
    convergence_percentage = float(n_pixels - (converges < 3).sum().cpu().numpy()) / float(n_pixels)
    
    return params, convergence_percentage, iteration


# END OF SOLVER



# PARAMETER GUESSING

def get_uniform_guess(data,uniform_param_guess):
    guess = torch.zeros(len(uniform_param_guess),data.size()[1]).cuda()
    for i in range(len(uniform_param_guess)):
        guess[i,:] = uniform_param_guess[i]
        
    return guess



# END PARAMETER GUESSING



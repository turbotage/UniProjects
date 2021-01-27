# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 17:05:26 2021

@author: TurboTage
"""

def get_mean_guess(diffusion_model,data,b_vals):
    n_pixels = data.shape[0]
    arr = [random.randint(0,n_pixels-1) for i in range(10)]
    
    param = [0,0]
    initial_params = np.array([1000, 0.002])
    
    for i in arr:
        s = data[i, :]
        result = curve_fit(diffusion_model, b_vals, s, initial_params)
        param = param + result[0]
        #print(result[0])
        
    param = param / 5
    
    print(param)
    
    guess = torch.zeros(2,data.shape[0]).cuda()
    guess[0,:] = param[0].item()
    guess[1,:] = param[1].item()
    return guess


def get_guess(data, imgs, b_vals):
    log_S = torch.log(data)
    max_index = np.argmax(b_vals)
    min_index = np.argmin(b_vals)
    
    ADC = (log_S[max_index,:] - log_S[min_index,:])/(b_vals[min_index]-b_vals[max_index])
    
    S0 = torch.exp(log_S[max_index,:] + b_vals[max_index]*ADC)
    
    guess = torch.zeros(2,data.size()[1]).cuda()
    guess[0,:] = S0
    guess[1,:] = ADC
    return guess
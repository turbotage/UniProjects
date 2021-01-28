# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 23:34:19 2021

@author: TurboTage
"""

import numpy as np
import time

import torch
import matplotlib.pyplot as plt #for plotting
from scipy.optimize import curve_fit

from torchvision import datasets, transforms
import random


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


def get_func(model,dependent):
    return lambda params: model(dependent, params)

# END OF MODEL


# START OF SOLVER

def get_jacobians(f, params):
    JT = torch.autograd.functional.jacobian(f, params).transpose(0,2)
    J = JT.transpose(1,2)
    return J, JT

def get_diff(func, params, data):
    diff = (data - func(params)).transpose(0,1).reshape(data.size()[1], 1, data.size()[0]).transpose(1,2)
    return diff

def get_step(J, JT, JT_diff, guess_order):
    # AX=B
    JTJ = torch.bmm(JT,J) # (JTJ)
    
    ill_behaved = (torch.logical_or(torch.det(JTJ) < 1, torch.det(JTJ) > 1e30)).reshape(JTJ.size()[0],1,1)
    
    eye = torch.eye(JTJ.size()[1]).repeat(JTJ.size()[0],1).reshape(JTJ.size()[0],JTJ.size()[1],JTJ.size()[1]).cuda()
    
    perturb = eye * guess_order * ill_behaved * 1e10
    
    JTJ = JTJ *((ill_behaved * 1e-30) + torch.logical_not(ill_behaved)) + perturb
    
    JT_diff = JT_diff * ((ill_behaved * 1e-30) + torch.logical_not(ill_behaved))
    
    
    # Solve via Cholskey decomposition, in theory better than LU for this problem
    u = torch.cholesky(JTJ)
    param_step = torch.cholesky_solve(JT_diff,u).transpose(0,2)[0]
    
    return param_step


def solve(model, dependent, data, guess, guess_order, threshold, max_iter):
    
    func = get_func(model,dependent)
    f = lambda p: func(p).sum(dim=1)
    
    n_params = guess.size()[0]
    n_pixels = guess.size()[1]
    
    non_convergence_list = torch.arange(0,n_pixels).cuda()
    
    
    # Initial Step
    J, JT = get_jacobians(f, guess)
    diff = get_diff(func, guess, data)
    JT_diff = torch.bmm(JT, diff)
    step = get_step(J, JT, JT_diff, guess_order)
    
    old_step = torch.zeros(n_params, n_pixels).cuda()
    converges = torch.abs((old_step - step) / guess) < threshold
    not_converging = converges.sum(dim=0) < 2

    params = guess
    #print("initial step done")
    
    iteration = 0
    while (converges.sum() != guess.size()[1]*guess.size()[0]) and iteration < max_iter:
        iteration += 1
        
        guess = guess + step
        
        J, JT = get_jacobians(f, guess)
        diff = get_diff(func,guess,data)
        
        JT_diff = torch.bmm(JT,diff)
        
        old_step = step
        step = get_step(J, JT, JT_diff, guess_order)
        
        converges = torch.abs((old_step - step) / guess) < threshold
        
        if iteration % 10 == 0:
            
            # copy converging pixels
            params[:,non_convergence_list[converges.sum(dim=0) > 1]] = guess[:,converges.sum(dim=0) > 1]
            
            
            # rebuild solver to solve with non-converging pixels
            not_converging = converges.sum(dim=0) < 2
            
            dependent = dependent[:,:,not_converging]
            data = data[:,not_converging]
            guess = guess[:,not_converging]
            step = step[:,not_converging]
            converges = converges[:,not_converging]
            
            func = get_func(model,dependent)
            f = lambda p: func(p).sum(dim=1)
            
            # rebuild non convergence list
            non_convergence_list = non_convergence_list[not_converging]
            
            #print(converges)
            #print(step)
        
        
    params[:,non_convergence_list[converges.sum(dim=0) > 1]] = guess[:,converges.sum(dim=0) > 1]
    
    
    convergence_percentage = float(n_pixels - (converges.sum(dim=0) < 2).sum().cpu().numpy()) / float(n_pixels)
    
    return params, convergence_percentage, iteration

       
# END OF SOLVER
        
       



# GET PARAMETER GUESS

def get_uniform_guess(data,uniform_param_guess):
    guess = torch.zeros(len(uniform_param_guess),data.size()[1]).cuda()
    for i in range(len(uniform_param_guess)):
        guess[i,:] = uniform_param_guess[i]
        
    return guess



# END PARAMETER GUESS





# LOAD AND TIDY DATA

def load_data(file_name: str):
    data = np.load(file_name)
    return data['img'], data['b_vals'].astype(np.float32)


def project_to_vector(img: np.ndarray, mask_image: np.ndarray) -> np.ndarray:
    """
    Project an image down to a vector. Inverse of 'project_to_image'.
    :param img: An image of the format [image_shape, n_images] where image shape can be any
    :param mask_image: A binary mask of shape [image_shape]
    :return: A set of vectors as a matrix of shape [sum(mask_image), n_images]
    """
    return img[mask_image, :]


def project_to_image(vector: np.ndarray, mask_image: np.ndarray) -> np.ndarray:
    """
    Project a vector onto an image. Inverse of 'project_to_vector'.
    :param vector: A set of vectors as a matrix of shape [sum(mask_image), n_images]
    :param mask_image: A binary mask of shape [image_shape]
    :return: An image of the format [image_shape, n_images] where the shape is derived from the mask_image
    """
    shape = [*mask_image.shape, *vector.shape[1:]]
    img = np.zeros(shape, dtype=vector.dtype)
    img[mask_image, ...] = vector
    return img




def plot(params, imgs):
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
    im1 = ax1.imshow(imgs[:,:,13,0], vmin=0, vmax=1320)
    im1 = ax2.imshow(imgs[:,:,13,1], vmin=0, vmax=1320)
    im1 = ax3.imshow(imgs[:,:,13,2], vmin=0, vmax=1320)
    im1 = ax4.imshow(imgs[:,:,13,3], vmin=0, vmax=1320)
    plt.show()
    
    # View the result
    fig, (ax1, ax2) = plt.subplots(1, 2)
    im1 = ax1.imshow(parameters[:, :, 13, 0], vmin=0, vmax=1000)
    fig.colorbar(im1, ax=ax1, orientation='horizontal')
    im2 = ax2.imshow(1000 * parameters[:, :, 13, 1], vmin=0, vmax=4)
    fig.colorbar(im2, ax=ax2, orientation='horizontal')
    ax1.set_title("S0 [a.u.]")
    ax2.set_title("ADC [x 10^-3 mm^2/s]")
    plt.show()


if __name__ == "__main__":
    # Setup CUDA
    torch.cuda.set_device(0)
    print('Used GPU Name:', torch.cuda.get_device_name(torch.cuda.current_device()))
    print()
    
    
    # Load the data
    imgs, b_vals = load_data("img_1.npz")

    # If the signal is too low it is not the subject. It is only the background
    bg_threshold = 150
    mask = imgs[:, :, :, 0] > bg_threshold

    # The data is now shaped as (92, 92, 25, 5). It is more convenient to have it on the shape (n_pixels x 5)
    # where n_pixels are the number of pixels in the mask.
    
    initial_guess = [500,0.001]
    
    data = torch.from_numpy(project_to_vector(imgs, mask)).transpose(0,1).cuda()
    guess = get_uniform_guess(data, initial_guess)
    dependent = torch.from_numpy(b_vals).repeat(data.size()[1],1).transpose(0,1).reshape(1, data.size()[0], data.size()[1]).cuda()
    
    
    model_expr = "P0 * exp(-X0 * P1)"
    model = get_model(model_expr)
    
    guess_order = torch.diag(torch.tensor(initial_guess)).cuda()
    
    start = time.time() 

    found_params, conv_perc, iterations = solve(model, dependent, data, guess, guess_order, 0.001, 50)
    
    end = time.time()
    
    print(conv_perc * 100, "% of pixels converges")
    print()
    
    print("Number of iterations:", iterations)
    print()
    # Place the result back as a 4D image with the 3 first dims as xyz and the last dim as the two
    # parameters s0 and adc.
    parameters = found_params.transpose(0,1).cpu().numpy()
    parameters = project_to_image(parameters, mask)
    
    plot(parameters,imgs)
    

    time_elapsed = end - start

    print(time_elapsed, "s")
    



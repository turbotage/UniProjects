# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 23:34:19 2021

@author: TurboTage
"""

import numpy as np
import time

import torch
import matplotlib.pyplot as plt #for plotting

from torchvision import datasets, transforms


def setup_device():
    torch.cuda.set_device(0)
    print('Used GPU Name:', torch.cuda.get_device_name(torch.cuda.current_device()))
    

    
def model(dependent):
    return lambda params: params[:][0] * torch.exp(-dependent * params[:][1])

def get_guess(data, imgs, b_vals):
    log_S = torch.log(data)
    max_index = np.argmax(b_vals)
    min_index = np.argmin(b_vals)
    
    ADC = (log_S[max_index,:] - log_S[min_index,:])/(b_vals[min_index]-b_vals[max_index])
    print(ADC.size())
    
    S0 = torch.exp(log_S[max_index,:] + b_vals[max_index]*ADC)
    print(S0.size())
    
    guess = torch.zeros(2, data.size()[1]).cuda()
    guess[0,:] = S0
    guess[1,:] = ADC
    return guess


def get_jacobians(func, params):
    t1 = time.time()
    #JT = torch.autograd.functional.jacobian(func,params).diagonal(dim1=-1).transpose(0,2)
    JT = torch.autograd.functional.jacobian(func, params).sum(dim=1).transpose(0,2)
    J = JT.transpose(1,2)
    t2 = time.time()
    print(t2 - t1)
    return J, JT

def get_diff(func, params, data):
    diff = (data - func(params)).transpose(0,1).reshape(data.size()[1], 1, data.size()[0]).transpose(1,2)
    return diff

def get_step(J, JT, error):
    # AX=B
    A = torch.bmm(JT,J) # (JTJ)
    
    # Solve via Cholskey decomposition, in theory better than LU for this problem
    u = torch.cholesky(A)
    param_step = torch.cholesky_solve(error,u).transpose(0,2)[0]
    
    return param_step


def solve(func,guess,data):
    iteration = 0
    threshold = 10^(-4)

    J, JT = get_jacobians(func, guess)
    diff = get_diff(func,guess,data)
    
    error = torch.bmm(JT,diff)
    
    error_sum = (error*error).sum(dim=1)
    
    while (error_sum < threshold).sum() != data.size()[1] and iteration < 5:
        iteration += 1
        
        step = get_step(J, JT, error)
        
        guess = guess + step
        
        J, JT = get_jacobians(func, guess)
        diff = get_diff(func,guess,data)
        
        error = torch.bmm(JT,diff)
        error_sum = (error*error).sum(dim=1)
        
        print(iteration)
        
    return guess
        
        
    
def generate_data(func, params):
    data = func(params)
    
    max_value = data.abs().max() * 0.08
    noise = (torch.rand(data.size()[0],data.size()[1]) - 0.5) * max_value
    
    return data + noise



"""
P = 2
M = 5
D = 10
    
# (data points, pixels)
params = torch.rand(P,M)*10

dependent = (torch.rand(D,M)-0.5)*0.1

func = model(dependent)

data = generate_data(func,params)

guess = torch.ones(P,M) * params.mean()
"""



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




if __name__ == "__main__":
    # Setup CUDA
    setup_device()
    
    
    # Load the data
    imgs, b_vals = load_data("img_1.npz")

    # If the signal is too low it is not the subject. It is only the background
    bg_threshold = 150
    mask = imgs[:, :, :, 0] > bg_threshold

    # The data is now shaped as (92, 92, 25, 5). It is more convenient to have it on the shape (n_pixels x 5)
    # where n_pixels are the number of pixels in the mask.
    data = torch.from_numpy(project_to_vector(imgs, mask)).transpose(0,1).cuda()
    
    guess = get_guess(data, imgs, b_vals)
    guess2 = guess[:,:500]
    
    data2 = data[:,:500]

    b_vals_torch = torch.from_numpy(b_vals).cuda()
    b_vals_torch = b_vals_torch.repeat(data.size()[1],1).transpose(0,1)
    
    b_vals_torch2 = b_vals_torch[:,:500]
    
    func = model(b_vals_torch2)
    
    start = time.time() 
    # Do the fitting
    found_params = solve(func, guess2, data2)
    
    end = time.time()

    # Place the result back as a 4D image with the 3 first dims as xyz and the last dim as the two
    # parameters s0 and adc.
    #parameters = project_to_image(parameters, mask)

    # View the result
    #fig, (ax1, ax2) = plt.subplots(1, 2)
    #im1 = ax1.imshow(parameters[:, :, 13, 0], vmin=0, vmax=1000)
    #fig.colorbar(im1, ax=ax1, orientation='horizontal')
    #im2 = ax2.imshow(1000 * parameters[:, :, 13, 1], vmin=0, vmax=4)
    #fig.colorbar(im2, ax=ax2, orientation='horizontal')
    #ax1.set_title("S0 [a.u.]")
    #ax2.set_title("ADC [x 10^-3 mm^2/s]")
    #plt.show()

    time_elapsed = end - start

    print(time_elapsed)





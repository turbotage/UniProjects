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


# START OF SOLVER

def get_jacobians(f, params):
    JT = torch.autograd.functional.jacobian(f, params).transpose(0,2)
    J = JT.transpose(1,2)
    return J, JT

def get_diff(func, params, data):
    diff = (data - func(params)).transpose(0,1).reshape(data.size()[1], 1, data.size()[0]).transpose(1,2)
    return diff

def get_step(J, JT, JT_diff):
    # AX=B
    JTJ = torch.bmm(JT,J) # (JTJ)
    #print(JTJ)
    
    #print((torch.det(JTJ) < 0.0000001).sum())
    
    ill_behaved = (torch.logical_or(torch.det(JTJ) < 1, torch.det(JTJ) > 1e30)).reshape(JTJ.size()[0],1,1)
    
    eye = torch.eye(JTJ.size()[1]).repeat(JTJ.size()[0],1).reshape(JTJ.size()[0],JTJ.size()[1],JTJ.size()[1]).cuda()
    
    perturb = eye * ill_behaved * 2
    #print(perturb[14056])
    
    print(JTJ[14056])
    
    
    JTJ = JTJ *((ill_behaved * 1e-20) + torch.logical_not(ill_behaved)) - perturb
    
    
    print(JTJ[14056])
    
    #print(JTJ.size())
    
    #print(JTJ)
    
    #print((torch.det(JTJ) < 0.0000001).sum())
    
    # Solve via Cholskey decomposition, in theory better than LU for this problem
    #u = torch.cholesky(JTJ)
    #param_step = torch.cholesky_solve(JT_diff,u).transpose(0,2)[0]
    
    # Use LU for now
    param_step, LU = torch.solve(JT_diff, JTJ)
    
    return param_step


def solve(func,data,guess):
    
    f = lambda params: func(params).sum(dim=1)
    
    iteration = 0
    threshold = 10^(-6)
    
    # Initial Step
    J, JT = get_jacobians(f, guess)
    diff = get_diff(func,guess,data)
    
    JT_diff = torch.bmm(JT,diff)
    
    error_sum = (JT_diff*JT_diff).sum(dim=1)
    
    #print("initial step done")
    
    while (error_sum < threshold).sum() != data.size()[1] and iteration < 10:
        iteration += 1
        
        step = get_step(J, JT, JT_diff)
        
        guess = guess + step
        
        J, JT = get_jacobians(f, guess)
        diff = get_diff(func,guess,data)
        
        JT_diff = torch.bmm(JT,diff)
        error_sum = (JT_diff*JT_diff).sum(dim=1)
        
        print(iteration)
        
    return guess
        
       
# END OF SOLVER
        
       

    
def model(dependent):
    return lambda params: params[:][0] * torch.exp(-dependent * params[:][1])



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


def get_bad_guess(data,imgs,b_vals):
    guess = torch.zeros(2,data.size()[1]).cuda()
    guess[0,:] = 100
    guess[1,:] = 0.001
    return guess


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




def setup_device():
    torch.cuda.set_device(0)
    print('Used GPU Name:', torch.cuda.get_device_name(torch.cuda.current_device()))
    


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
    guess = get_bad_guess(data, imgs, b_vals)
    dependent = torch.from_numpy(b_vals).repeat(data.size()[1],1).transpose(0,1).cuda()
    
    
    
    func = model(dependent)
    
    start = time.time() 

    found_params = solve(func, data, guess)
    
    end = time.time()
    
    
    # Place the result back as a 4D image with the 3 first dims as xyz and the last dim as the two
    # parameters s0 and adc.
    parameters = found_params.transpose(0,1).cpu().numpy()
    parameters = project_to_image(parameters, mask)
    
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

    time_elapsed = end - start

    print(time_elapsed)
    



class Ob:
    def __init__(self,num):
        self.num = num


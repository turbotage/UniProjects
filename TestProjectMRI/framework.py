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

import random

import gsolver_lm

        
       



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


# END OF LOAD DATA


def exp_diffusion_model():
    # Load the data
    imgs, b_vals = load_data("img_1.npz")

    # If the signal is too low it is not the subject. It is only the background
    bg_threshold = 150
    mask = imgs[:, :, :, 0] > bg_threshold

    # The data is now shaped as (92, 92, 25, 5). It is more convenient to have it on the shape (n_pixels x 5)
    # where n_pixels are the number of pixels in the mask.
    
    initial_guess = [800,0.004]
    
    data = torch.from_numpy(project_to_vector(imgs, mask)).transpose(0,1).cuda()
    
    guess = gsolver_lm.get_uniform_guess(data, initial_guess)
    
    dependent = torch.from_numpy(b_vals).repeat(data.size()[1],1).transpose(0,1).reshape(1, data.size()[0], data.size()[1]).cuda()
    
    model_expr = "P0 * exp(-abs(X0 * P1))"
    model = gsolver_lm.Model(model_expr, dependent, data, guess)
    
    start = time.time() 
    
    
    
    found_params, conv_perc, iterations = model.solve(0.0001, 300)
    
    end = time.time()
    
    time_elapsed = end - start

    print(time_elapsed, "s")
    print()
    
    
    print(conv_perc * 100, "% of pixels converges")
    print()
    
    print("Number of iterations:", iterations)
    print()
    
    
    # Place the result back as a 4D image with the 3 first dims as xyz and the last dim as the two
    # parameters s0 and adc.
    parameters = found_params.transpose(0,1).cpu().numpy()
    parameters = project_to_image(parameters, mask)

    
    # Plot
    #fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
    #im1 = ax1.imshow(imgs[:,:,13,0], vmin=0, vmax=1320)
    #im1 = ax2.imshow(imgs[:,:,13,1], vmin=0, vmax=1320)
    #im1 = ax3.imshow(imgs[:,:,13,2], vmin=0, vmax=1320)
    #im1 = ax4.imshow(imgs[:,:,13,3], vmin=0, vmax=1320)
    #plt.show()
    
    # View the result
    fig, (ax1, ax2) = plt.subplots(1, 2)
    im1 = ax1.imshow(np.abs(parameters[:, :, 13, 0]), vmin=0, vmax=1000)
    fig.colorbar(im1, ax=ax1, orientation='horizontal')
    im2 = ax2.imshow(1000 * np.abs(parameters[:, :, 13, 1]), vmin=0, vmax=4)
    fig.colorbar(im2, ax=ax2, orientation='horizontal')
    ax1.set_title("S0 [a.u.]")
    ax2.set_title("ADC [x 10^-3 mm^2/s]")
    plt.show()
    


def vfa_model():
    # Load the data
    imgs, fa = load_data("VFA_dataset.npz")

    # If the signal is too low it is not the subject. It is only the background
    bg_threshold = 150
    mask = imgs[:, :, :, 0] > bg_threshold

    # The data is now shaped as (92, 92, 25, 5). It is more convenient to have it on the shape (n_pixels x 5)
    # where n_pixels are the number of pixels in the mask.
    
    initial_guess = [300,0.005]
    
    data = torch.from_numpy(project_to_vector(imgs, mask)).transpose(0,1).cuda()
    
    guess = gsolver_lm.get_uniform_guess(data, initial_guess)
    
    dependent = torch.zeros(2, data.size()[0], data.size()[1]).cuda()
    
    
    
    model_expr = "P0 * exp(-X0 * P1)"
    model = gsolver_lm.get_model(model_expr)
    
    guess_order = torch.diag(torch.tensor(initial_guess)).cuda()
    
    start = time.time() 

    found_params, conv_perc, iterations = gsolver_lm.solve(model, dependent, data, guess, guess_order, 0.0001, 300)
    
    end = time.time()
    
    time_elapsed = end - start

    print(time_elapsed, "s")
    print()
    
    
    print(conv_perc * 100, "% of pixels converges")
    print()
    
    print("Number of iterations:", iterations)
    print()
    
    
    # Place the result back as a 4D image with the 3 first dims as xyz and the last dim as the two
    # parameters s0 and adc.
    parameters = found_params.transpose(0,1).cpu().numpy()
    parameters = project_to_image(parameters, mask)

    
    # Plot
    #fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
    #im1 = ax1.imshow(imgs[:,:,13,0], vmin=0, vmax=1320)
    #im1 = ax2.imshow(imgs[:,:,13,1], vmin=0, vmax=1320)
    #im1 = ax3.imshow(imgs[:,:,13,2], vmin=0, vmax=1320)
    #im1 = ax4.imshow(imgs[:,:,13,3], vmin=0, vmax=1320)
    #plt.show()
    
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
    
    
    exp_diffusion_model()
    
    



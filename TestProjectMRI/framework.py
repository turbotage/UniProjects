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
    
    
def load_data(file_name: str):
    data = np.load(file_name)
    return data['img'], data['b_vals']



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
    


def model(dependent):
    return lambda parameters: parameters[0] * torch.exp(-dependent * parameters[1]) + 2*parameters[2]
        

def get_jacobian(model, current_param, n_params, m_points):
    JT = torch.autograd.functional.vjp(func,current_param, torch.ones(m_points).cuda())[1]
    J = torch.transpose(JT, 0, 1)
    return (J, JT)
        
    
#number of datapoints
m = 4


func = model(-1)

param1 = torch.ones(10).cuda()
param2 = torch.ones(10).cuda() * 2
param3 = torch.ones(10).cuda() * 3

params = torch.rand(3,10).cuda()
#params[0] = param1
#params[1] = param2
#params[2] = param3

v1 = torch.ones(4).cuda()

#l = torch.autograd.functional.vjp(func, d, torch.ones(3).cuda())


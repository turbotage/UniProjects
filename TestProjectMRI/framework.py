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
    return lambda parameters: parameters[:][0] * torch.exp(-dependent * parameters[:][1]) + 2*parameters[:][2]
        


def get_jacobian(model, current_param, m_cuda_ones):
    """
    Generates the jacobian and it's transpose to be used for calculating next step,
    in an iterative solver, example Newton-Rapson...

    Parameters
    ----------
    model : function
        The function to calculate the jacobian for
    current_param : TYPE
        the values of the parameters that the jacobian shall be differentiated with respect to
        a K x M tensor where M is the number of data points and K is the number of parameters 
    m_cuda_ones : TYPE
        a M (number of data points) long tensor of ones, it is passed so that no 
        new tensors are created unessesarily.

    Returns
    -------
    J : TYPE
        The jacobian evaluated at current_param
    JT : TYPE
        The transpose of the jacobian at current_param
    """
    JT = torch.autograd.functional.vjp(func,current_param, m_cuda_ones)[1]
    J = torch.transpose(JT, 0, 1)
    return (J, JT)
   

def solve(A, b):
    """
    Uses x = inv(A)*b to solve Ax=b

    Parameters
    ----------
    A : TYPE
        A in Ax=b
    b : TYPE
        DESCRIPTION.
        b in Ax=b
    Returns
        x in Ax=b
    """
    return torch.mm(A.inverse(), b)


def solve_LU(A, b):
    """
    Uses torch.solve for now which utilizes LU decomposition

    Parameters
    ----------
    A : TYPE
        A in Ax=b
    b : TYPE
        DESCRIPTION.
        b in Ax=b
    Returns
        x in Ax=b
    """
    return torch.solve(b, A)
    

def solve_Cholesky(A, b):
    """
    Uses torch. for now which utilizes LU decomposition, in future
    consider using Cholesky instead since it can be about 2 times faster and
    the problem to be solved is real hermitian

    Parameters
    ----------
    A : TYPE
        A in Ax=b
    b : TYPE
        b in Ax=b

    Returns
        x in Ax=b
    """
    u = torch.cholesky(A)
    return torch.cholesky_solve(b, u)
    

def get_next_step(func, current_params,data, m_cuda_ones):
    J,JT = get_jacobian(func,current_params, m_cuda_ones)
    LM = (J*J).sum(dim=1)
    DIFF = data - func(current_params)
    return JT*DIFF/LM
   

def solve_iter(func, guess, data):
    return 2
    
    

#number of datapoints




data = torch.rand(5,2*2*2).cuda()
dependent = torch.rand(5,2*2*2).cuda()

func = model(dependent)

params = torch.ones(3,2*2*2).cuda()
#params[0] = param1
#params[1] = param2
#params[2] = param3

#l = torch.autograd.functional.vjp(func, d, torch.ones(3).cuda())


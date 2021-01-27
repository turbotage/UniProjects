# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 13:15:49 2021

@author: TurboTage
"""

import numpy as np
import time

import torch
import matplotlib.pyplot as plt #for plotting

from torchvision import datasets, transforms


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
    
    data = torch.from_numpy(project_to_vector(imgs, mask)).transpose(0,1).cuda()
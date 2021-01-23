import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import time


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


# The model use in the curve-fitting
def diffusion_model(b, s0, adc):
    return s0 * np.exp(-b * adc)


# Fit a single pixel
def fit_diffusion_pixel(b, s, initial_params):
    curve_fit(diffusion_model, b, s, initial_params)


# Print iterations progress - from Stack-overflow
def printProgressBar(iteration, total, prefix='', suffix='', decimals=1, length=100, fill='█', printEnd="\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end='', flush=True)
    # Print New Line on Complete
    if iteration == total:
        print()


# Fit all pixels one-by-one
def fit_diffusion_data(b, data):
    n_pixels = data.shape[0]
    params = np.zeros((n_pixels, 2))

    # Uniform initial guess for all pixels
    initial_params = np.array([1000, 0.002])

    for i in range(n_pixels):
        s = data[i, :]
        result = curve_fit(diffusion_model, b, s, initial_params)
        params[i, :] = result[0]
        if i % 1000 == 0:
            printProgressBar(i, n_pixels, 'Progress', 'Complete', length=40)
    return params


def load_data(file_name: str):
    data = np.load(file_name)
    return data['img'], data['b_vals']


if __name__ == "__main__":
    # Load the data
    imgs, b_vals = load_data("img_1.npz")

    # If the signal is too low it is not the subject. It is only the background
    bg_threshold = 150
    mask = imgs[:, :, :, 0] > bg_threshold

    # The data is now shaped as (92, 92, 25, 5). It is more convenient to have it on the shape (n_pixels x 5)
    # where n_pixels are the number of pixels in the mask.
    data = project_to_vector(imgs, mask)
    
    print(b_vals)
    print(data)

    start = time.time()
    # Do the fitting
    parameters = fit_diffusion_data(b_vals, data)
    end = time.time()

    # Place the result back as a 4D image with the 3 first dims as xyz and the last dim as the two
    # parameters s0 and adc.
    parameters = project_to_image(parameters, mask)

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
    print(f"Speed: {data.shape[0] / time_elapsed} fit voxels/s")

import matplotlib.pyplot as plt
import numpy as np


def plot_correction_factors(coeff: np.ndarray, uncorrected: np.ndarray, corrected: np.ndarray):
    ''' plot correction factors from `front_tracking_toolkit.preprocess.preprocess.calculate_correction_factors()`

    Parameters:
    coeff (np.ndarray): correction coefficients
    uncorrected (np.ndarray): uncorrected data
    corrected (np.ndarray): corrected data

    Returns:
    fig, axs: matplotlib figure and array of axes
    '''
    fig, axs = plt.subplots(2, 1)

    axs[0].plot(uncorrected, label='Uncorrected')
    axs[0].plot(corrected, label='Corrected')
    axs[0].set_ylabel('Intensity Sum')
    axs[0].legend()

    axs[1].axhline(y=1.0, color='gray', linestyle=':')
    axs[1].plot(coeff)
    axs[1].set_ylabel('Correction Factor')

    return fig, axs

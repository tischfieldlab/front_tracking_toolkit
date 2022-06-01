from typing import Tuple, Union

import numpy as np
import statsmodels.api as sm
from pystackreg import StackReg
from tqdm import tqdm

from front_tracking_toolkit.io.image import load_image
from front_tracking_toolkit.stats import max_exclude_outliers



def crop_images(images: np.ndarray, x1: int = None, y1: int = None, x2: int = None, y2: int = None) -> np.ndarray:
    ''' Simple rectangular crop of images

    Parameters:
    images (np.ndarray): Array of image data, of shape (nframes, y, x)
    x1 (int): x position of top left corner of box to keep
    y1 (int): y position of top left corner of box to keep
    x2 (int): x position of bottom right corner of box to keep
    y2 (int): y position of bottom right corner of box to keep

    Returns:
    np.ndarray: cropped image series
    '''
    x1 = 0 if x1 is None else x1
    y1 = 0 if y1 is None else y1
    x2 = images.shape[2] if x2 is None else x2
    y2 = images.shape[1] if y2 is None else y2
    return images[:, y1:y2, x1:x2]


def mask_images(images: np.ndarray, mask: Union[str, np.ndarray]):
    ''' Apply a mask to images

    Pixels where mask is zero are set to zero in `images`, otherwise values are retained

    Parameters:
    images (np.ndarray): Array of image data, of shape (nframes, y, x)
    mask_path (str|np.ndarray): Path to mask image or directly a array of mask data.
    '''
    if isinstance(mask, str):
        mask = load_image(mask)
    elif isinstance(mask, np.ndarray):
        pass
    else:
        raise ValueError(f'Could not interpret parameter mask: {mask}')

    mask[mask > 0] = 1
    images = images[:, ...] * mask
    return images


def calculate_correction_factors(images) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    ''' Calculate correction factors for intensity fluctuations

    Calculates deviation from a lowess fit, and normalization factors
    that would correct the deviation

    Parameters:
    images (np.ndarray): array of image data, of shape (nframes, y, x)

    Returns
    coefficients (np.ndarray), image means (np.ndarray), corrected image means (np.ndarray)
    '''
    sums = np.mean(np.array(images), axis=(1, 2))
    lowess = sm.nonparametric.lowess(sums, np.arange(images.shape[0]), frac=0.2)
    coeff = lowess[:, 1] / sums
    coeff[coeff < 0] = 0

    return (coeff, sums, sums * coeff)


def apply_correction_factors(coeff: np.ndarray, images: np.ndarray) -> np.ndarray:
    ''' Apply factors calculated by `calculate_correction_factors()`
    '''
    corrected_images = images * coeff[:, None, None]
    corrected_images[corrected_images > 255] = 255
    corrected_images[corrected_images < 0] = 0
    return corrected_images.astype('uint8')


def correct_intensity_fluctuations(images: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    '''
    In the case there are artifactual intensity fluctuations

    '''
    factors = calculate_correction_factors(images)
    #plot_norm_factors(*factors)
    corrected_images = apply_correction_factors(factors[0], images)
    return corrected_images, factors


def normalize_intensity(images: np.ndarray, vmin: float = None, vmax: float = None, exclude_outliers: bool = True,
                        outlier_threshold: float = 3.5) -> np.ndarray:
    ''' Normalize image intensity by performing a linear stretch

    max of the dataset can be computed using raw pixel intensities
    or by discarding outlier pixels

    Parameters:
    images (np.ndarray): Array of image data, of shape (nframes, y, x)
    vmin (Union[float, None]): minimum data value to map to
    exclude_outliers (bool): If True, discard outliers when computing data max
    outlier_threshold (float): Threshold for outlier detection, Median Absolute Deviation

    Returns:
    np.ndarray: normalize image data
    '''
    images = images.astype('float32')
    if vmax is None:
        if exclude_outliers:
            datamax = max_exclude_outliers(images.ravel(), threshold=outlier_threshold)
        else:
            datamax = images.max()
    else:
        datamax = vmax

    if vmin is None:
        datamin = images.min()
    else:
        datamin = vmin

    print(f'    -> Stretching from [{datamin}, {datamax}] -> [0, 255]')

    norm = ((images - datamin) / (datamax - datamin)) * 255
    norm[norm > 255] = 255
    norm[norm < 0] = 0
    return norm.astype('uint8')


def register_images(images: np.ndarray, **kwargs) -> np.ndarray:
    ''' Run image registration to correct for movement

    Uses the rigid body method of pystackreg

    We reverse images in the time dimention for registration, as it seems to provide more stable results
    when the first few frames have little to no signal

    Parameters:
    images (np.ndarray): images to register
    **kwargs: additional arguments to pass to `StackReg.register_transform_stack()`

    Returns:
    np.ndarray: movement corrected image series
    '''
    progress = tqdm(total=images.shape[0], desc='Registering Images', leave=False)
    def update_progress(current_iteration: int, end_iteration: int):
        if progress.total != end_iteration:
            progress.reset(total=end_iteration)
        progress.n = current_iteration
        progress.refresh()

    stkreg = StackReg(StackReg.RIGID_BODY)
    registered = stkreg.register_transform_stack(np.flip(images, axis=0), progress_callback=update_progress, **kwargs)
    registered = np.flip(registered, axis=0)
    registered[registered > 255] = 255
    registered[registered < 0] = 0
    return registered.astype('uint8')

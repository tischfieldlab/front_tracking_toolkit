from typing import Tuple
import numpy as np
from skimage.draw import polygon2mask, polygon_perimeter


def poly_area(x: np.ndarray, y: np.ndarray) -> float:
    ''' get the area of a polygon defined by verticies x and y

    Parameters:
    x (np.ndarray): x coordinates of polygon verticies
    y (np.ndarray): y coordinates of polygon verticies

    Returns:
    float: area of the polygon
    '''
    # https://stackoverflow.com/a/30408825
    return 0.5 * np.abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))


def poly_to_px_coords(x: np.ndarray, y: np.ndarray, scale: float = 1.0) -> np.ndarray:
    ''' convert polygon from real units to pixel units given `scale`

    Parameters:
    x (np.ndarray): x coordinates of polygon verticies
    y (np.ndarray): y coordinates of polygon verticies
    scale (float): conversion factor from real units to pixels

    Returns:
    np.ndarray: polygon coordinates scaled from real units to pixel units
    '''
    x = x / scale
    y = y / scale
    return np.swapaxes(np.array([x, y]), 0, 1)


def poly_to_mask(x: np.ndarray, y: np.ndarray, shape: Tuple[int, int], scale: float = 1.0) -> np.ndarray:
    ''' Convert a polygon to a mask

    Parameters:
    x (np.ndarray): x coordinates of polygon verticies
    y (np.ndarray): y coordinates of polygon verticies
    shape (Tuple[int, int]): shape of the output mask, (x, y)
    scale (float): conversion factor from real units to pixels

    Returns:
    np.ndarray: mask
    '''
    coords = poly_to_px_coords(x, y, scale)
    msk = polygon2mask((shape[1], shape[0]), coords)
    msk = np.flipud(np.rot90(msk))
    return msk


def draw_poly(img: np.ndarray, coords: np.ndarray, color: Tuple[int, int, int] = None) -> np.ndarray:
    ''' Draw a polygon onto an image

    Parameters:
    img (np.ndarray): Image to draw upon
    coords (np.ndarray): polygon coordinates, of shape (verticies, xy)
    color (Tuple[int, int, int]): color to draw the polygon, tuple of ints in range 0-255, RGB order

    Returns:
    np.ndarray: image with polygon drawn
    '''
    msk = polygon_perimeter(coords[:, 0], coords[:, 1], shape=img.shape)
    if img.shape == 3:
        if color is None:
            color = (255, 0, 0) # Red
        # RGB image
        img[msk[0], msk[1], 0] = color[0]
        img[msk[0], msk[1], 1] = color[1]
        img[msk[0], msk[1], 2] = color[2]
    else:
        # grayscale image
        if color is None:
            color = 255 # White
        img[msk[0], msk[1]] = color
    return img

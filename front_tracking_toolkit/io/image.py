import glob
import os
import xml.etree.ElementTree as ET
from typing import Iterable, Tuple, Union, List

import cv2
import numpy as np
import tifffile
from skimage import color, img_as_ubyte


def find_images(pattern: Union[str, Iterable[str]]) -> List[str]:
    ''' Find images given a pattern

    Patten may be a glob-able string, or a list of glob-able strings. Returned results are sorted

    Parameters:
    pattern (Union[str, Iterable[str]]): string or iterable of string paths to search. glob wildcards allowed

    Returns:
    sorted list of paths
    '''
    paths = []
    if isinstance(pattern, str):
        paths.extend(glob.glob(pattern, recursive=True))
    else:
        for subpattern in pattern:
            paths.extend(glob.glob(subpattern, recursive=True))
    return sorted(paths)


def load_images(image_paths: Iterable[str]) -> np.ndarray:
    ''' Load a series of images defined by `image_paths`

    see `load_image()`

    Parameters:
    image_paths (Iterable[str]): Paths of the images to load

    Returns:
    np.ndarray: array of image data, of shape (nframes, y, x), single channel grayscale uint8 dtype
    '''
    images = []
    for path in image_paths:
        images.append(load_image(path))
    return np.array(images)


def load_image(path: str) -> np.ndarray:
    ''' Load a single tiff image from `path`

    Image returned is single channel grayscale as uint8 dtype

    Parameters:
    path (str): path to the image to read

    Returns:
    np.ndarray: array of image data, of shape (y, x)
    '''
    image = tifffile.imread(path)
    if image.shape[-1] == 3:
        image = color.rgb2gray(image)
    image = np.asarray(img_as_ubyte(image))
    return image


def save_images(images: np.ndarray, dest_dir: str) -> None:
    ''' Save an image series to disk

    Parameters:
    images (np.ndarray): Array of image data, of shape (nframes, y, x)
    dest_dir (str): destination directory for images
    '''
    os.makedirs(dest_dir, exist_ok=True)
    nimages = images.shape[0]
    fmt = '{:0' + str(int(np.ceil(np.log10(max(1, abs(nimages)+1))))) + 'd}.tif'
    for i in range(nimages):
        path = os.path.join(dest_dir, fmt.format(i+1))
        tifffile.imwrite(path, img_as_ubyte(images[i]))


def save_video(images: np.ndarray, dest_dir: str, basename: str = 'video', codec: str = 'MP4V', fps: int = 6) -> None:
    ''' Save an image series as a video, in grayscale and false color heatmap

    Parameters:
    images (np.ndarray): Array of image data, of shape (nframes, y, x)
    dest_dir (str): Destination directory for videos
    basename (str): basename for the video files
    codec (str): codec to use for writing videos
    fps (int): frame per second for the videos
    '''
    foutrcc = cv2.VideoWriter_fourcc(*codec)
    width = images.shape[2]
    height = images.shape[1]

    if codec in ['MP4V']:
        ext = 'mp4'
    elif codec in ['XVID']:
        ext = 'avi'

    os.makedirs(dest_dir, exist_ok=True)

    dest_gray = os.path.join(dest_dir, f'{basename}.gray.{ext}')
    video_gray = cv2.VideoWriter(dest_gray, foutrcc, fps, (width, height), isColor=False)

    dest_color = os.path.join(dest_dir, f'{basename}.color.{ext}')
    video_color = cv2.VideoWriter(dest_color, foutrcc, fps, (width, height), isColor=True)

    for i in range(images.shape[0]):
        video_gray.write(images[i])
        video_color.write(cv2.applyColorMap(images[i], cv2.COLORMAP_JET))
    video_gray.release()
    video_color.release()


def save_color_video(images: np.ndarray, dest_dir: str, basename: str = 'video', codec: str = 'MP4V', fps: int = 6) -> None:
    ''' Save an image series as a color video

    Parameters:
    images (np.ndarray): Array of image data, of shape (nframes, y, x, c)
    dest_dir (str): Destination directory for videos
    basename (str): basename for the video files
    codec (str): codec to use for writing videos
    fps (int): frame per second for the videos
    '''
    foutrcc = cv2.VideoWriter_fourcc(*codec)
    width = images.shape[2]
    height = images.shape[1]

    if codec in ['MP4V']:
        ext = 'mp4'
    elif codec in ['XVID']:
        ext = 'avi'

    os.makedirs(dest_dir, exist_ok=True)
    dest = os.path.join(dest_dir, '{}.{}'.format(basename, ext))
    video = cv2.VideoWriter(dest, foutrcc, fps, (width, height), isColor=True)

    for i in range(images.shape[0]):
        video.write(cv2.cvtColor(images[i], cv2.COLOR_RGB2BGR))
    video.release()


def get_scale_for_image(path: str) -> Tuple[float, str]:
    ''' Try to find scale information for an image

    Currently implemented:
      - Leica images, with metadata stored in './.Metadata/*.cal.xml

    If no scale information could be found, the tuple `(1.0, 'unknown')` is returned

    Parameters:
    path (str): path to image for which scale data should be found

    Returns:
    Tuple[float, str]: scale and units
    '''
    # try for Leica metadata
    dirname = os.path.dirname(path)
    base = os.path.basename(path)
    meta = os.path.join(dirname, '.Metadata', f'{base}.cal.xml')
    if os.path.exists(meta):
        return (read_scale_info(meta), 'mm')

    return (1.0, 'unknown')


def read_scale_info(path: str) -> float:
    ''' Read scale information from Leica metadata

    Parameters:
    path (str): path to Leica metadata xml file

    Returns:
    float: pixel size in mm
    '''
    namespace = 'http://schemas.datacontract.org/2004/07/LeicaMicrosystems.DataEntities.V3_2'
    tree = ET.parse(path)
    res = tree.find(".//{"+namespace+"}XMetresPerPixel").text
    return float(res) * 1e3

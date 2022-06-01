import cv2
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
# from front_tracking_toolkit.experiment.experiment import Experiment
from front_tracking_toolkit.geometry import poly_to_mask, poly_to_px_coords
from front_tracking_toolkit.io.image import save_color_video


def make_front_areas_movie(images: np.ndarray, fronts: pd.DataFrame, dest_dir: str, basename: str, scale: float) -> np.ndarray:
    ''' Make a movie of front areas as heatmaps overlain on images

    Parameters:
    images (np.ndarray): images to draw upon
    fronts (pd.DataFrame): Front data
    dest_dir (str): destination directory
    basename (str): basename for the movie
    scale (float): scale factor for real units to pixel units

    Returns:
    np.ndarray: movie frames, of shape (nframes, y, x, 3[rgb])
    '''
    masks = np.zeros(images.shape, dtype='bool')

    # iterate (time, front-index) groups in fronts
    for (t, _), group in fronts:
        idx = int(t/60)
        masks[idx, :, :] |= poly_to_mask(group['x'], group['y'], masks.shape[1:], scale)
    masks = np.repeat(masks[..., np.newaxis], 3, axis=3)

    out_images = np.zeros((*images.shape, 3), dtype='uint8')
    for i in range(out_images.shape[0]):
        gray_as_color = cv2.cvtColor(images[i], cv2.COLOR_GRAY2RGB)
        gray_as_heat = cv2.cvtColor(cv2.applyColorMap(images[i], cv2.COLORMAP_JET), cv2.COLOR_BGR2RGB)
        out_images[i] = np.where(masks[i], gray_as_heat, gray_as_color)

    save_color_video(out_images, dest_dir=dest_dir, basename=basename)
    return out_images



def make_front_lines_movie(images: np.ndarray, fronts: pd.DataFrame, dest_dir: str, basename: str, scale: float) -> np.ndarray:
    ''' Make a movie of front lines overlaid on images

    Parameters:
    images (np.ndarray): images to draw upon
    fronts (pd.DataFrame): Front data
    dest_dir (str): destination directory
    basename (str): basename for the movie
    scale (float): scale factor for real units to pixel units

    Returns:
    np.ndarray: movie frames, of shape (nframes, y, x, 3[rgb])
    '''

    out_images = np.zeros((images.shape[0], images.shape[1], images.shape[2], 3), dtype='uint8')
    for i in range(out_images.shape[0]):
        im = cv2.cvtColor(images[i], cv2.COLOR_GRAY2RGB)

        # iterate (time, front-index) groups in fronts
        for (t, _), group in fronts:
            if i == int(t/60):
                coords = poly_to_px_coords(group['x'], group['y'], scale)
                coords = [coords.reshape((-1, 1, 2)).astype(np.int32)]
                im = cv2.polylines(im, coords, True, (255, 0, 0), 10, cv2.LINE_AA)
        out_images[i] = im

    save_color_video(out_images, dest_dir=dest_dir, basename=basename)
    return out_images


# def plot_spatial_velocity_hexbin(exp: Experiment, fig_size=(8, 25)):
#     ''' Plot spatial velocity as a hexbin plot
#     '''
#     fig, axs = plt.subplots(len(exp.subjects), len(exp.stains), figsize=fig_size, sharex=True, sharey=True)

#     vmin = exp.data['vm'].min()
#     vmax = exp.data['vm'].max()

#     for spi, subject in enumerate(exp.data['subject'].unique()):
#         for sti, stain in enumerate(exp.config.stains):
#             subset = exp.data[(exp.data['subject'] == subject) & (exp.data['stain'] == stain)].dropna(subset=['vm'])

#             ax = axs[spi, sti]
#             ax.set_title(f'{subject} ({exp.metadata[subject]["Genotype"]}): {stain}')
#             im = ax.hexbin(subset['x'], subset['y'], C=subset['vm'], vmin=vmin, vmax=vmax)
#             ax.set_aspect('equal')

#     fig.subplots_adjust(right=0.8)
#     cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
#     fig.colorbar(im, cax=cbar_ax)

#     return fig, axs

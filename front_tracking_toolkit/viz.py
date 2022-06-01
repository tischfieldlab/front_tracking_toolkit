from typing import Iterable, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
from IPython.display import HTML
from matplotlib import animation
from matplotlib.colorbar import Colorbar
from mpl_toolkits.axes_grid1.axes_grid import ImageGrid



def visualize_frames(images: Union[np.ndarray, Iterable[np.ndarray]], vmin: float = None, vmax: float = None, cmap: str = 'viridis', interval: int = 100,
                     figsize: Tuple[float, float] = None, labels: Union[Iterable[str], str] = None) -> HTML:
    ''' Visualize frames as a matplotlib animation for jupyter notebooks

    If frames are RGB (i.e. len(images.shape) == 4), images are displayed as RGB. Otherwise, images are displayed using `cmap` and `vmin`/`vmax`.

    Parameters:
    images (Union[np.ndarray, Iterable[np.ndarray]]): one or more image series to visualize
    vmin (float): minimum value of colormap
    vmax (float): maximum value of colormap
    cmap (str): colormap to use for plotting
    interval (int): number of milliseconds between animation frames
    figsize (Tuple[float, float]): size of the figure
    labels (Iterable[str]): Labels for the images, should be same length as images

    Returns:
    HTML: matplotlib html5 video animation wrapped in IPython HTML
    '''
    if isinstance(images, np.ndarray):
        nplots = 1
        images = [images]
    else:
        nplots = len(images)

    if isinstance(labels, str):
        labels = [labels]
    elif labels is None:
        labels = []

    if figsize is None:
        figsize = (6.0 * nplots, 4.5)

    nframes = np.max([imgs.shape[0] for imgs in images])
    vmin = np.min([imgs.min() for imgs in images]) if vmin is None else vmin
    vmax = np.max([imgs.max() for imgs in images]) if vmax is None else vmax


    fig = plt.figure(figsize=figsize)
    suptitle = fig.suptitle('')
    grid = ImageGrid(
        fig=fig,
        rect=111,
        nrows_ncols=(1, nplots),
        cbar_location="right",
        cbar_mode="single",
        axes_pad=0.15,
    )

    ax_imgs = []
    needs_cbar = False

    for i, ax in enumerate(grid):
        try:
            ax.set_title(labels[i])
        except IndexError:
            pass

        if len(images[i].shape) == 3:
            ax_imgs.append(ax.imshow(images[i][0], vmin=vmin, vmax=vmax, cmap=cmap))
            needs_cbar = True

        else:
            ax_imgs.append(ax.imshow(images[i][0]))

    if needs_cbar:
        Colorbar(grid[0].cax, ax_imgs[0])

    plt.close(fig)

    def init():
        for i, ax in enumerate(ax_imgs):
            suptitle.set_text('Frame 0')
            ax.set_data(images[i][0])

    def animate(frame):
        for i, ax in enumerate(ax_imgs):
            suptitle.set_text(f'Frame {frame}')
            ax.set_data(images[i][frame])
        return grid

    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=nframes, interval=interval)
    return HTML(anim.to_html5_video())

import os
from typing import Iterable

import numpy as np
from front_tracking_toolkit.io.data import write_yaml
from front_tracking_toolkit.io.image import (find_images, load_images,
                                             save_images, save_video)
from front_tracking_toolkit.preprocess.config import PreprocessOptions
from front_tracking_toolkit.preprocess.proc import (
    correct_intensity_fluctuations, crop_images, mask_images,
    normalize_intensity, register_images)
from front_tracking_toolkit.preprocess.viz import plot_correction_factors


def preprocess_images(image_paths: Iterable[str], opts: PreprocessOptions) -> np.ndarray:
    ''' Process an image series given `opts`

    See `PreprocessOptions()`

    Parameters:
    image_paths (Iterable[str]): Image paths (in order) to be processed
    opts (PreprocessOptions): options for preprocessing pipeline
    '''
    #print(f"working on {path}")

    if not opts.dry:
        meta_dir = os.path.join(opts.dest, 'metadata')
        os.makedirs(meta_dir, exist_ok=True)

        write_yaml(os.path.join(meta_dir, 'options.yaml'), opts.dict())

    # load images
    print("  - Loading images...")
    images = load_images(find_images(image_paths))
    print(f"    -> Found and loaded {images.shape[0]} images") # pylint: disable=unsubscriptable-object

    # crop images according to options
    if opts.crop.enabled:
        print("  - Cropping images...")
        images = crop_images(images, **opts.crop.kwargs.dict())

    #mask images according to options
    if opts.mask.enabled:
        print("  - masking images...")
        images = mask_images(images, opts.mask.mask)

    # normalize image intensity (linear stretch)
    if opts.inorm.enabled:
        print("  - Normalizing image intensity")
        images = normalize_intensity(images, **opts.inorm.kwargs.dict())

    # correct for intensity fluctuations
    if opts.icorr.enabled:
        print("  - Correcting intensity fluctuations...")
        images, factors = correct_intensity_fluctuations(images)
        fig, _ = plot_correction_factors(*factors)
        if not opts.dry:
            fig.savefig(os.path.join(meta_dir, 'intensity_correction.png'))
            fig.savefig(os.path.join(meta_dir, 'intensity_correction.pdf'))
        #print(images.dtype)

    # correction for movement
    if opts.reg.enabled:
        print("  - Correcting for movement...")
        images = register_images(images, **opts.reg.kwargs.dict())
        #print(images.dtype)

    if not opts.dry:
        # save images
        print("  - Saving images...")
        save_images(images, opts.dest)

        # save a video for easier inspection
        print("  - Saving video...")
        save_video(images, opts.dest)

    return images

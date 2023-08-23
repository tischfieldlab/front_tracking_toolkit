from typing import Literal, Optional, Union

import numpy as np

from front_tracking_toolkit.util import FTBaseModel


class CropOptionsKwargs(FTBaseModel):
    ''' Kwargs for crop Options
    '''
    x1: Optional[int] = None
    y1: Optional[int] = None
    x2: Optional[int] = None
    y2: Optional[int] = None

class CropOptions(FTBaseModel):
    ''' Options for crop step
    '''
    enabled: bool = False
    kwargs: CropOptionsKwargs = CropOptionsKwargs()

class MaskOptions(FTBaseModel):
    ''' Options for masking step
    '''
    enabled: bool = False
    mask: Union[None, str, np.ndarray] = None

class INormOptionsKwargs(FTBaseModel):
    ''' Kwargs for intensity normalization Options
    '''
    exclude_outliers: bool = False
    outlier_threshold: float = 3.5
    vmin: Optional[float] = None
    vmax: Optional[float] = None

class INormOptions(FTBaseModel):
    ''' Options for intensity normalization step
    '''
    enabled: bool = True
    kwargs: INormOptionsKwargs = INormOptionsKwargs()

class ICorrOptions(FTBaseModel):
    ''' Options for intensity correction stop
    '''
    enabled: bool = False

class RegOptionsKwargs(FTBaseModel):
    ''' Kwargs for Registration Options
    '''
    reference: Literal['previous', 'first', 'mean'] = 'previous'
    n_frames: int = 1
    moving_average: int = 1

class RegOptions(FTBaseModel):
    ''' Options for image registration step
    '''
    enabled: bool = True
    kwargs: RegOptionsKwargs = RegOptionsKwargs()



class PreprocessOptions(FTBaseModel):
    ''' preprocessing options
    '''
    dry: bool = False
    dest: Optional[str] = None
    crop: CropOptions = CropOptions()
    inorm: INormOptions = INormOptions()
    icorr: ICorrOptions = ICorrOptions()
    reg: RegOptions = RegOptions()
    mask: MaskOptions = MaskOptions()


# def get_default_processing_opts() -> PreprocessOptions:
#     ''' Get a dictionary with default preprocessing options

#     Returns:
#     dict: containing default processing options
#     '''
#     return PreprocessOptions()
#     # return {
#     #     'dest': None,
#     #     'dry': False, # No saving output
#     #     'crop': {
#     #         'enabled': False,
#     #         'kwargs': {}
#     #     },
#     #     'inorm': {
#     #         'enabled': True,
#     #         'kwargs': {
#     #             'vmin': None,
#     #             'vmax': None,
#     #             'exclude_outliers': False,
#     #             'outlier_threshold': 3.5,
#     #         }
#     #     },
#     #     'icorr': {
#     #         'enabled': True
#     #     },
#     #     'reg': {
#     #         'enabled': True,
#     #         'kwargs': {
#     #             'reference': 'previous'
#     #         }
#     #     },
#     #     'mask': {
#     #         'enabled': False,
#     #         'mask': None
#     #     },
#     # }

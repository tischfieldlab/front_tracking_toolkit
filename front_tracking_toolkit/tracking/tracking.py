import inspect
import io
import math
import os
from pathlib import Path
from typing import Iterable, List, Literal, Union

import front_tracking_toolkit

import matlab
import matlab.engine
from front_tracking_toolkit.io.data import write_yaml

from front_tracking_toolkit.util import FTBaseModel


def get_matlab_paths() -> List[str]:
    ''' Get default paths that should be added to the matlab engine
    Here we get paths for front tracking libs and the main entrypoint
    '''
    mat_dir = Path(inspect.getabsfile(front_tracking_toolkit)).parent.parent.joinpath('matlab')

    return [
        mat_dir,
        os.path.join(mat_dir, 'FrontTrackingWithFlowRelease'),
        os.path.join(mat_dir, 'helpers')
    ]


class MatlabSession(object):
    ''' Context manager for a matlab session
    '''
    def __init__(self, paths: Iterable[str] = None, verbose: bool = False) -> None:
        ''' Construct a new matlab engine

        Parameters:
        paths (Iterable[str]): paths to add to matlab engine
        '''
        self.verbose = verbose
        self.engine = matlab.engine.start_matlab()
        self.paths = get_matlab_paths()
        if paths is not None:
            self.paths.extend(paths)
        for path in self.paths:
            self.add_path(path)

    def add_path(self, path: str):
        ''' Add a path to the matlab engine

        Parameters:
        path (str): path to add to matlab engine
        '''
        if self.engine is None:
            raise RuntimeError('Engine is not yet initialized!')

        if self.verbose:
            print(f'Adding "{path}" to matlab path')
        self.engine.addpath(str(path))

    def __enter__(self):
        return self.engine

    def __exit__(self, _type, value, traceback):
        self.engine.quit()


class FrontTrackOptions(FTBaseModel):
    ''' Typed class for front tracking options
    '''
    threshold: Union[int, Literal['auto']] = 50
    step: int = 1
    minarea: int = 25
    span: int = 50
    ROI: List = []
    framerange: List = [-math.inf, math.inf]
    frrate: float = 1/60
    invert: int = 0
    noisy: int = 0
    thick: int = 1
    maxSpeed: float = 1.5
    maxPx: int = 40
    bac: Union[int, Literal['first']] = 0
    scale: float = 1.0



def _prepare_front_tracking_options(options: FrontTrackOptions) -> dict:
    ''' Ensure all configuration items are the correct types
    '''
    return {
        'threshold': 'auto' if options.threshold == 'auto' else int(options.threshold),
        'step': float(options.step),
        'minarea': float(options.minarea),
        'span': float(options.span),
        'ROI': matlab.double(options.ROI),
        'framerange': matlab.double(options.framerange),
        'frrate': float(options.frrate),
        'invert': float(options.invert),
        'noisy': float(options.noisy),
        'thick': float(options.thick),
        'maxSpeed': float(options.maxSpeed),
        'maxPx': float(options.maxPx),
        'bac': 'first' if options.bac == 'first' else float(options.bac),
        'scale': float(options.scale),
    }


def run_front_tracking(images: Iterable[str], dest: str, config: FrontTrackOptions = None) -> None:
    ''' Run front tracking on images

    Parameters:
    images (Iterable[str]): paths to images to use for front tracking
    dest (str): destination directory for results
    config (dict): configuration for front tracking
    '''
    dest = os.path.abspath(dest)
    os.makedirs(dest, exist_ok=True)

    if config is None:
        config = FrontTrackOptions()
    write_yaml(os.path.join(dest, 'options.yaml'), config.dict())
    final_config = _prepare_front_tracking_options(config)

    cap = MatlabConsoleCapture()
    with MatlabSession() as engine:
        engine.ProcessSample(images, dest, final_config, nargout=0, stdout=cap.out, stderr=cap.err)

    cap.save(os.path.join(dest, 'front-tracking'))

    return cap




class MatlabConsoleCapture(object):
    ''' Utility class to help capture stdout and stderr from Matlab
    '''

    def __init__(self) -> None:
        self.out = io.StringIO()
        self.err = io.StringIO()

    def save(self, dest_base: str) -> None:
        ''' Save stderr and stdout streams to a file

        Parameters:
        dest_base (str): Destination base name (path plus basename)
        '''
        with open(f'{dest_base}.out.txt', 'w', encoding='utf8') as out_f:
            out_f.write(self.out.getvalue())

        with open(f'{dest_base}.err.txt', 'w', encoding='utf8') as err_f:
            err_f.write(self.err.getvalue())

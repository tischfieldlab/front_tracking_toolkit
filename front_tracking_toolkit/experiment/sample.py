

from enum import Enum, unique
from typing import Dict, Iterable, List, Tuple, Union

import numpy as np
import pandas as pd
from tabulate import tabulate
from front_tracking_toolkit.experiment.config import PipelineConfig
from front_tracking_toolkit.io.image import find_images, get_all_meta_for_leica_image, load_image, load_images

@unique
class ImageStage(str, Enum):
    ''' Image Stage
    '''
    RAW = 'raw'
    PREPROCESSED = 'preprocessed'


class ImageSeries(object):
    ''' Image Series
    '''

    def __init__(self, stage: ImageStage, path_spec: Union[str, Iterable[str], None] = None) -> None:
        self.stage: ImageStage = ImageStage(stage)
        self.path_spec: Union[str, Iterable[str], None] = path_spec
        self.paths: List[str] = []
        self.refresh()

    def refresh(self):
        ''' Refresh image list from path_spec '''
        self.paths.clear()
        if self.path_spec is not None:
            self.paths = find_images(self.path_spec)

    def load_image_metadata(self) -> pd.DataFrame:
        meta = []
        for i, img_path in enumerate(self.paths):
            m = get_all_meta_for_leica_image(img_path)
            m['filename'] = img_path
            m['index'] = i
            if i == 0:
                m['AcquiredDelta'] = pd.to_timedelta(0)
            else:
                m['AcquiredDelta'] = m['AcquiredDate'] - meta[i-1]['AcquiredDate']
            meta.append(m)
        return pd.DataFrame(meta)

    def __len__(self):
        return len(self.paths)


class Sample(object):
    ''' Represents a single sample within an Experiment
    '''

    def __init__(self, subject: str, stain: str) -> None:
        self._subject = str(subject)
        self._stain = str(stain)
        self.images: Dict[str, ImageSeries] = {}
        self.mask: Union[None, np.ndarray] = None
        self.scale: Tuple[float, str] = (1.0, 'unknown')
        self.metadata: dict = {}
        self.pipeline: PipelineConfig = PipelineConfig()

    @property
    def subject(self) -> str:
        ''' Get this sample's subject '''
        return self._subject

    @property
    def stain(self) -> str:
        ''' Get this sample's stain '''
        return self._stain

    def register_images(self, stage: ImageStage, path_spec: Union[str, Iterable[str]]):
        ''' Register images
        '''
        self.images[stage] = ImageSeries(stage=stage, path_spec=path_spec)

    def get_images(self, stage: ImageStage, t: Union[range, int] = None) -> List[str]:
        ''' Get a list of image file paths for a given stage and time range

        Parameters:
        stage (str): which stage of images to retrieve (raw or processed)
        t (Union[range, int, None]): if passed, restrict data to these indicies
        '''
        if t is None:
            return self.images[stage].paths
        else:
            return self.images[stage].paths[t]


    def load_images(self, stage: ImageStage = 'raw', t: Union[range, int] = None) -> np.ndarray:
        ''' Load image(s) for a given stage and time range

        Parameters:
        stage (str): which stage of images to retrieve (raw or processed)
        t (Union[range, int, None]): if passed, restrict data to these indicies
        '''
        image_paths = self.get_images(stage=stage, t=t)

        if isinstance(image_paths, str):
            return load_image(image_paths)
        else:
            return load_images(image_paths)

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}('{self.subject}', '{self.stain}')"

    def describe(self):
        ''' Describe this experiment
        '''
        data = [
            ["Attribute", "Value"],
            ["Subject", self.subject],
            ["Stain", self.stain],
            ["Has Mask", str(self.mask is not None)],
            ["Scale", f"{self.scale[0]:.3E} {self.scale[1]}"]
        ]
        for stage, imseries in self.images.items():
            data.append([f'# {stage} images', len(imseries)])

        print(self.__repr__())
        print(tabulate(data, headers='firstrow'))
        print('\n')
        print("Subject Metadata:")
        print(tabulate(self.metadata.items(), headers=["Attribute", "Value"], showindex=False))

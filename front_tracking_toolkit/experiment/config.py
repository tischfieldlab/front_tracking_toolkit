import os
from typing import List

from tabulate import tabulate

from front_tracking_toolkit.preprocess.config import PreprocessOptions
from front_tracking_toolkit.tracking.tracking import FrontTrackOptions
from front_tracking_toolkit.util import FTBaseModel


class PipelineConfig(FTBaseModel):
    ''' Configuration for Pipeline
    '''
    preprocess: PreprocessOptions = PreprocessOptions()
    tracking: FrontTrackOptions = FrontTrackOptions()


class ExperimentConfig(FTBaseModel):
    ''' Configuration data for an Experiment
    '''
    basedir: str = os.getcwd()
    tracking_pattern: str = 'tracking/{subject}/{stain}/vf.csv'
    processed_images_pattern: str = 'processed/{subject}/{stain}/*.tif'
    tracking_image_stage: str = 'processed'

    def __str__(self) -> str:
        out = 'Experiment Configuration:\n'
        out += tabulate(vars(self).items(), headers=['Property', 'Value'])
        return out


class ExperimnetDefinition(FTBaseModel):
    ''' Experiment Definition '''
    config: ExperimentConfig = ExperimentConfig()
    pipeline: PipelineConfig = PipelineConfig()
    samples: List = []

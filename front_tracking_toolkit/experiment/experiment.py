import os
from distutils.log import warn
from typing import List, Union

import numpy as np
import pandas as pd
from tabulate import tabulate
import tqdm
from front_tracking_toolkit.experiment.config import ExperimentConfig, PipelineConfig
from front_tracking_toolkit.experiment.sample import ImageStage, Sample

from front_tracking_toolkit.io.data import read_yaml
from front_tracking_toolkit.io.image import (find_images, get_scale_for_image)
from front_tracking_toolkit.preprocess.config import PreprocessOptions
from front_tracking_toolkit.preprocess.preprocess import preprocess_images
from front_tracking_toolkit.tracking.tracking import FrontTrackOptions, run_front_tracking
from front_tracking_toolkit.tracking.viz import make_front_areas_movie, make_front_lines_movie


class Experiment(object):
    ''' Represents a front tracking experiment
    '''

    def __init__(self, config: str) -> None:
        super().__init__()
        self.config_file = config
        self.config: ExperimentConfig = None
        self.pipeline: PipelineConfig = None
        self._sample_map: dict = None
        self._read_experiment_definition()

    def reload(self):
        ''' Reload this experiment
        '''
        self._read_experiment_definition()

    @property
    def samples(self) -> List[Sample]:
        ''' Get the samples in this experiment
        '''
        return list(self._sample_map.values())

    def get_sample(self, subject: str, stain: str) -> Sample:
        ''' Get the sample with subject and stain
        '''
        return self._sample_map[(subject, stain)]

    def has_sample(self, subject: str, stain: str) -> bool:
        ''' tell if this experiment has a sample with subject and stain
        '''
        return (subject, stain) in self._sample_map


    @property
    def subjects(self) -> List[str]:
        ''' Get unique subjects in this experiment, in no particular order
        '''
        return list(self.metadata.index.levels[0])


    @property
    def stains(self) -> List[str]:
        ''' Get unique stains in this experiment, in no particular order
        '''
        return list(self.metadata.index.levels[1])

    @property
    def metadata(self) -> pd.DataFrame:
        ''' Get a dataframe with sample metadata from this experiment
        '''
        metadata_items = []
        for sample in self.samples:
            metadata_items.append({
                'subject': sample.subject,
                'stain': sample.stain,
                'scale_factor': sample.scale[0],
                'scale_units': sample.scale[1],
                **sample.metadata,
            })
        df = pd.DataFrame(metadata_items)
        df = df.set_index(['subject', 'stain'])
        return df

    def get_images(self, subject: str, stain: str, stage: ImageStage = 'raw', t: Union[range, int] = None) -> List[str]:
        ''' Get a list of image file paths for a given subject, stain, stage, and timepoints
        '''
        return self.get_sample(subject, stain).get_images(t=t, stage=stage)

    def load_images(self, subject, stain, stage='raw', t: Union[range, int] = None) -> np.ndarray:
        ''' Load image(s) for a given subject and stain
        '''
        return self.get_sample(subject, stain).load_images(t=t, stage=stage)

    def get_effective_preprocess_options(self, sample: Sample, config: PreprocessOptions = None):
        ''' Get effective preprocess options given a sample and user config '''
        if config is None:
            config = {}

        cfg = PreprocessOptions()
        cfg.dest = os.path.dirname(sample.images[ImageStage.PREPROCESSED].path_spec)
        cfg.mask.mask = sample.mask
        return PreprocessOptions.merge(self.pipeline.preprocess, sample.pipeline.preprocess, cfg, config)


    def preprocess_sample(self, sample: Sample, config: PreprocessOptions = None) -> np.ndarray:
        ''' Run preprocessing on a sample
        '''

        cfg = self.get_effective_preprocess_options(sample=sample, config=config)
        processed_frames = preprocess_images(sample.get_images(stage='raw', t=None), opts=cfg)
        sample.proc_images = find_images(sample.images[ImageStage.PREPROCESSED].path_spec)

        return processed_frames


    def get_effective_tracking_options(self, sample: Sample, config: FrontTrackOptions = None):
        ''' Get effective tracking options given a sample and user config '''
        if config is None:
            config = {}

        cfg = FrontTrackOptions()
        cfg.scale = sample.scale[0]
        return FrontTrackOptions.merge(self.pipeline.tracking, cfg, sample.pipeline.tracking, config)


    def front_track_sample(self, sample: Sample, config: FrontTrackOptions = None):
        ''' Run front tracking on sample
        '''
        stage = 'preprocessed'
        cfg = self.get_effective_tracking_options(sample=sample, config=config)

        rel_dest = os.path.dirname(self.config.tracking_pattern.format(subject=sample.subject, stain=sample.stain))
        dest = os.path.join(self.config.basedir, rel_dest)

        run_front_tracking(sample.get_images(stage=stage), dest, cfg)

        fronts = self.load_tracking_results(sample)
        t_fi_grouped_fronts = fronts.groupby(['t', 'fi'])
        images = sample.load_images(stage=stage)
        lines_frames = make_front_lines_movie(images, t_fi_grouped_fronts, dest, 'lines_movie', sample.scale[0])
        areas_frames = make_front_areas_movie(images, t_fi_grouped_fronts, dest, 'areas_movie', sample.scale[0])

        return {
            'fronts': fronts,
            'lines_frames': lines_frames,
            'areas_frames': areas_frames,
        }


    def has_tracking_results(self, samples=None) -> bool:
        ''' Tell if tracking results exist '''
        if samples is None:
            samples = self.samples
        elif isinstance(samples, Sample):
            samples = [samples]


        for sample in samples:
            path = os.path.join(self.config.basedir, self.config.tracking_pattern).format(subject=sample.subject, stain=sample.stain)
            if not os.path.exists(path):
                return False
        return True


    def load_tracking_results(self, samples=None, incldue_meta: Union[None, List[str]] = None) -> pd.DataFrame:
        ''' Load tracking results
        '''
        if samples is None:
            samples = self.samples
        elif isinstance(samples, Sample):
            samples = [samples]

        dsets = []
        for sample in tqdm.tqdm(samples, desc='Loading Samples', leave=False):
            path = os.path.join(self.config.basedir, self.config.tracking_pattern).format(subject=sample.subject, stain=sample.stain)

            if os.path.exists(path):
                data = pd.read_csv(path)
                data['subject'] = sample.subject
                data['stain'] = sample.stain

                if incldue_meta is not None:
                    for meta in incldue_meta:
                        data[meta] = self.metadata.loc[(sample.subject, sample.stain), meta]

                dsets.append(data)
            else:
                warn(f'Warning: tracking results not found for subject {sample.subject} and stain {sample.stain}, expected to find results here: "{path}"')

        all_data = pd.concat(dsets, ignore_index=True)
        all_data['vm'] = np.sqrt(all_data['vx']**2 + all_data['vy']**2)
        return all_data


    def load_image_intensities(self, samples=None, stage='preprocessed', incldue_meta: Union[None, List[str]] = None) -> pd.DataFrame:
        ''' Load whole image intensities
        '''
        if samples is None:
            samples = self.samples
        elif isinstance(samples, Sample):
            samples = [samples]

        whole_img_intensities = []
        for sample in tqdm.tqdm(samples, desc='Loading Samples', leave=False):
            imgs = sample.load_images(stage=stage)
            means = np.mean(imgs, axis=(1,2))
            medians = np.median(imgs, axis=(1,2))
            sums = np.sum(imgs, axis=(1,2))
            stdev = np.std(imgs, axis=(1,2))

            if incldue_meta is not None:
                meta = {m: self.metadata.loc[(sample.subject, sample.stain), m] for m in incldue_meta}
            else:
                meta = {}

            for i in range(imgs.shape[0]):
                whole_img_intensities.append({
                    'subject': sample.subject,
                    'stain': sample.stain,
                    **meta,
                    't': i * 60,
                    'mean': means[i],
                    'median': medians[i],
                    'sum': sums[i],
                    'std': stdev[i],
                })

        return pd.DataFrame(whole_img_intensities)


    def load_image_metadata(self, samples=None, stage='raw') -> pd.DataFrame:
        if samples is None:
            samples = self.samples
        elif isinstance(samples, Sample):
            samples = [samples]

        meta = []
        for sample in tqdm.tqdm(samples, desc='Loading Samples', leave=False):
            m = sample.images[stage].load_image_metadata()
            m['subject'] = sample.subject
            m['stain'] = sample.stain
            meta.append(m)
        return pd.concat(meta, ignore_index=True)


    def _read_experiment_definition(self) -> None:
        ''' Read experiment definition from a yaml file
        '''

        config = read_yaml(self.config_file)
        self.config = ExperimentConfig(**config.get('config', {}))
        self.pipeline = PipelineConfig()
        if 'pipeline' in config:
            if 'preprocess' in config['pipeline'] and config['pipeline']['preprocess'] is not None:
                self.pipeline.preprocess = PreprocessOptions(**config['pipeline'].get('preprocess', {}))

            if 'tracking' in config['pipeline'] and config['pipeline']['tracking'] is not None:
                self.pipeline.tracking = FrontTrackOptions(**config['pipeline'].get('tracking', {}))


        all_samples = {}
        for i, sample in enumerate(config['samples']):
            try:
                sample['subject'] = str(sample['subject'])
                sample['stain'] = str(sample['stain'])
                sid = (sample['subject'], sample['stain'])
                if sid in all_samples:
                    raise ValueError(f'Error: subject "{sample["subject"]}" and stain "{sample["stain"]}" have already been parsed! duplicates are not allowed!')

                s = Sample(sample['subject'], sample['stain'])
                s.register_images(ImageStage.RAW, sample['images'])

                proc_img_base = self.config.processed_images_pattern.format(subject=sample['subject'], stain=sample['stain'])
                s.register_images(ImageStage.PREPROCESSED, os.path.join(self.config.basedir, proc_img_base))

                s.scale = tuple(sample.get('scale', get_scale_for_image(s.images[ImageStage.RAW].paths[0])))
                s.mask = sample.get('mask', None)
                s.metadata = sample.get('metadata', {})

                if 'pipeline' in sample:
                    if 'preprocess' in sample['pipeline']:
                        s.pipeline.preprocess = PreprocessOptions.merge(self.pipeline.preprocess, sample['pipeline']['preprocess'])

                    if 'tracking' in sample['pipeline']:
                        s.pipeline.preprocess = PreprocessOptions.merge(self.pipeline.preprocess, sample['pipeline']['tracking'])

                all_samples[sid] = s
            except Exception as ex:
                raise Exception(f'While reading sample #{i} (subject="{sample["subject"]}", stain="{sample["stain"]}"), an exception was encountered!') from ex

        self._sample_map = all_samples


    def describe(self):
        ''' Print information describing this Experiment
        '''
        print(self.config)
        print('\n')
        print("Experiment Metadata:")
        meta = self.metadata.reset_index()
        meta['# raw images'] = meta.apply(lambda row: len(self.get_images(row['subject'], row['stain'], stage=ImageStage.RAW)), axis=1)
        meta['# proc images'] = meta.apply(lambda row: len(self.get_images(row['subject'], row['stain'], stage=ImageStage.PREPROCESSED)), axis=1)
        meta['has tracking'] = meta.apply(lambda row: self.has_tracking_results(self.get_sample(row['subject'], row['stain'])), axis=1)
        print(tabulate(meta, headers='keys', showindex=False))


    def __str__(self) -> str:
        return self.metadata.to_string() + '\n' + self.config


    def __repr__(self) -> str:
        return f'{self.__class__.__name__}("{self.config_file}")'

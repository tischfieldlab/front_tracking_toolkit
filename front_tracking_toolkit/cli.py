import click
from front_tracking_toolkit.experiment.config import ExperimnetDefinition
from front_tracking_toolkit.io.data import write_yaml

from front_tracking_toolkit.util import click_monkey_patch_option_show_defaults


# Show click option defaults
click_monkey_patch_option_show_defaults()


@click.group()
@click.version_option()
def cli():
    ''' Toolbox for front tracking flouresence images
    '''
    pass # pylint: disable=unnecessary-pass


@cli.command('generate-experiment-config')
@click.argument('dest', default='experiment-definition.yaml', type=click.Path())
@click.option('--include-defaults', is_flag=True)
def generate_experiment_config(dest, include_defaults):
    ''' Generate an experiment definition yaml file
    '''
    write_yaml(dest, ExperimnetDefinition().dict(exclude_defaults=not include_defaults))

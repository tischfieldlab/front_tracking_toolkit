


import click
import mergedeep
from pydantic import BaseModel # pylint: disable=no-name-in-module


def click_monkey_patch_option_show_defaults():
    ''' Monkey patch click.core.Option to turn on showing default values.
    '''
    orig_init = click.core.Option.__init__
    def new_init(self, *args, **kwargs):
        ''' This version of click.core.Option.__init__ will set show default values to True
        '''
        orig_init(self, *args, **kwargs)
        self.show_default = True
    # end new_init()
    click.core.Option.__init__ = new_init # type: ignore


class FTBaseModel(BaseModel):
    ''' BaseModel with some config and extra methods
    '''

    @classmethod
    def merge(cls, *args):
        ''' Merge models
        '''
        opts = {
            'exclude_defaults': True
        }
        to_merge = [a.dict(**opts) if isinstance(a, BaseModel) else a for a in args]
        #print(to_merge)
        merged = mergedeep.merge({}, *to_merge)
        print(merged)
        return cls(**merged)

    class Config:
        ''' Configuration for FTBaseModel '''
        arbitrary_types_allowed = True
        validate_default = True
        validate_assignment = True
        extra = 'forbid'

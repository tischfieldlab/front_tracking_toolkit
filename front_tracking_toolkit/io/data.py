import ruamel.yaml as yaml


def read_yaml(yaml_file: str) -> dict:
    ''' Read a yaml file into dict object

    Parameters:
        yaml_file (str): path to yaml file

    Returns:
        return_dict (dict): dict of yaml contents
    '''
    with open(yaml_file, 'r', encoding='utf8') as in_file:
        yml = yaml.YAML(typ='safe')
        return yml.load(in_file)


def write_yaml(yaml_file: str, data: dict) -> None:
    ''' Write a dict object into a yaml file

    Parameters:
        yaml_file (str): path to yaml file
        data (dict): dict of data to write to `yaml_file`
    '''
    with open(yaml_file, 'w', encoding='utf8') as out_file:
        yml = yaml.YAML(typ='safe')
        yml.default_flow_style = False
        yml.dump(data, out_file)

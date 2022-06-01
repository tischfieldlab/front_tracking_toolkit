
from distutils.core import setup
import setuptools


setup(
    name="front-tracking-toolkit",
    version=1.0,
    license="MIT License",
    install_requires=[
        'click',
        'opencv-python',
        'numpy',
        'scikit-image',
        'seaborn',
        'matplotlib',
        'pandas',
        'statsmodels',
        'pystackreg',
        'ruamel.yaml',
        'mergedeep',
        'tabulate',
    ],
    description="Toolkit for analyzing front tracking experiments",
    packages=setuptools.find_packages(),
    include_package_data=True,
    entry_points={
        "console_scripts": [
            'front-tracking-toolkit = front_tracking_toolkit.cli:cli',
            'front-track = front_tracking_toolkit.cli:cli',
        ],
    },
    zip_safe=False,
)

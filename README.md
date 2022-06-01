

# Installation

Create a new conda environment
```
conda create -y -n front-tracking python=3.8
conda activate front-tracking
```

Install ffmpeg
```
conda install -c conda-forge ffmpeg
```

Install Matlab-Python bridge
Navigate to your matlab install and run the `setup.py` file there.
```
cd "c:\Program Files\MATLAB\R2021a\extern\engines\python"
python setup.py install
```
See https://www.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html for more details
Tested successfully with `R2021a` and `2021b`. Possibly some issues on mac using parallel-for loops


Install front-tracking-toolkit
```
pip install -e .
```
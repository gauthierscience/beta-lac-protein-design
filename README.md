# Introduction

This repository contains scripts and jupyter notebooks to generate and analyze
$\beta$-lactamase variants as well as analyze their experimental characterization 
as was done in our recent publication (see citation below):

Download this repository and navigate to the root directory
```
git clone https://github.com/gauthierscience/beta-lac-protein-design.git
cd beta-lac-protein-design
```

Notes:
- Data are available in the `analysis_scripts/data` directory. Analysis on new data must conform to the same format. Output of the scripts can be viewed as images in the jupyter notebooks.
- Analysis scripts are
all provided as jupyter notebooks in the `analysis_scripts` directory. Output of these scripts is shown as embedded in the jupyter notebook.
- There are no notable installation or run time committments for any of these scripts.
- All MATLAB scripts and data (multiple sequence alignment and model) that were used to generate the designs are available in
[potts_design_release.zip](https://github.com/gauthierscience/beta-lac-protein-design/blob/main/analysis_scripts/potts_design_release.zip). Please see `README.m` in the zip file for examples of how to execute new design generation.

## Demo / Run analysis scripts

1. **Optional** create a virtual environment for analysis:
```
python3 -m venv ~/env/betalacdesign
source ~/env/betalacdesign/bin/activate
```

2. Install dependencies
```
pip3 install --upgrade pip
pip3 install pandas
pip3 install biopython
pip3 install seaborn
pip3 install openpyxl
pip3 install scikit-learn
pip3 install https://github.com/debbiemarkslab/EVcouplings/archive/develop.zip
pip3 install jupyter
pip3 install lmfit
pip3 install "numpy<1.24"
```

3. Navigate to script directory, launch jupyter notebook and open any of the notebooks in the gui (ipynb files).
```
cd analysis_scripts
jupyter notebook
```

## Citation:

>**Simultaneous Enhancement of Multiple Functional Properties Using Evolution-informed Protein Design.**
>Benjamin Fram<sup>#</sup>,
>Yang Su<sup>\*</sup>, 
>Ian Truebridge<sup>\*</sup>,
>Adam J. Riesselman,
>John B. Ingraham,
>Alessandro Passera,
>Eve Napier,
>Nicole N. Thadani,
>Samuel Lim,
>Kristen Roberts,
>Gurleen Kaur,
>Michael A. Stiffler,
>Debora S. Marks,
>Christopher D. Bahl,
>Amir R. Khan,
>Chris Sander,
>Nicholas P. Gauthier<sup>#</sup>,
>Nature Communications, _accepted in principle_, 2024
>
> \# Correspondence should be addressed to [Benjamin Fram and Nicholas Gauthier](mailto:benjamin.fram.research@gmail.com,nicholas.gauthier.research@gmail.com)

### bioRxiv version available at https://doi.org/10.1101/2023.05.09.539914

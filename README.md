# Connectivity analysis in Alzheimer disease patients
## Installation
* prerequisites
    - python3
    - gfortran
    - liblapack
    - FSL
    - freesurfer
    - HD-BET
    - julia

```bash
git clone http::/github.com/Yutong441/ADcon
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt

cd ADcon/cdipy
python setup.py build_ext --inplace
cd -

cd ADcon/bct_for
python -m numpy.f2py -llapack -c -m fort *.f90
cd -
deactivate
```

## Analysis
Analysis is conducted in the following steps:

First, skullstrip T1 and FLAIR, register them to MNI space.
This requires specifying the dataset name, T1 and FLAIR paths regular expression
patterns and the output directory in the `data_raw/data_regex.json` file.

```bash
# specify the path to the data folder as $root
python data_raw/reg.py --data_name=EOAD_ADNI \
        --root=$root
```

Next, run the connectivity analysis
```bash
python ADcon/lesion.py \
    --root_dir=$root/ADNI/processed/EOAD_ADNI \
    --save_dir=$save_dir \
    --atlas_dir=data
```

Finally, use `data_raw/load_data.R` to summarize all the data.
Use `data_raw/quantify.R` to summarize the WMH data.
Use `data_raw/VLSM.jl` to perform VSLM analysis.

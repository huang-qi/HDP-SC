# HDP-SC

## Project Overview

**HDP-SC** is a project designed for building a coarse-grained potential for polymer. 

## Project Structure

# HDP-SC/
├── README.md
├── LICENSE
├── data/
│   ├── AA/
│   │   ├── AA.data
│   │   └── AA.lammpstrj
│   ├── CG/
│   │   ├── atom_types.npy
│   │   ├── box.npy
│   │   ├── build.json
│   │   ├── CG.data
│   │   ├── CG.lammpstrj
│   │   ├── force.npy
│   │   ├── pos.npy
│   │   ├── prior.json
│   │   └── structure.psf
│   └── DP/
│       ├── input.json
│       ├── training_data/
│       │   ├── set.000/
│       │   │   ├── box.npy
│       │   │   ├── coord.npy
│       │   │   └── force.npy
│       │   └── type.raw
│       └── validation_data/
│           ├── set.000/
│           │   ├── box.npy
│           │   ├── coord.npy
│           │   └── force.npy
│           └── type.raw
├── model/
│   └── graph.pb
└── scripts/
    ├── do_cg_map.ipynb
    ├── model_train.ipynb
    ├── prior_fit.ipynb
    ├── rdf_analysis.ipynb
    └── utils.py
    
##  Guide

1. **CG Mapping**

   Run `do_cg_map.ipynb` to generate the coarse-grained mapping.

2. **Prior Energy Fitting**

   Run `prior_fit.ipynb` to generate prior parameters.

3. **Deep Potential Training**

   Run `model_train.ipynb` to train the deep potential model.

4. **RDF Analysis**

   Run `rdf_analysis.ipynb` to perform RDF analysis.

## Dependencies

All project dependencies are listed in `requirements.txt`. The main dependencies include:

- deepmd-kit
- numpy
- scipy
- mdtraj
- moleculekit
- matplotlib
- scikit-learn

## License

This project is licensed under the Apache License. See the [LICENSE](LICENSE) file for details.

## Contact

If you have any questions or suggestions, please contact me via email: huangqi@mail.sim.ac.cn.

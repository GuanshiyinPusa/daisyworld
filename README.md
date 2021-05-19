# Interactive Daisyworld

<p align="center">
  <img src="https://raw.githubusercontent.com/mrernst/daisyworld/main/imgs/phone1.png" width="100">
  <img src="https://raw.githubusercontent.com/mrernst/daisyworld/main/imgs/phone2.png" width="100">
  <img src="https://raw.githubusercontent.com/mrernst/daisyworld/main/imgs/phone3.png" width="100">
  

Interactive daisyworld is a python implementation of the famous [daisyworld computer simulation](https://en.wikipedia.org/wiki/Daisyworld). It describes a hypothetical world orbiting a star, whose radiant energy slowly increases over time. The world has two types of daisies (black and white) as it's only life form. The population of daisies can modify the worlds albedo and thereby migitage the effects arising from the increasing radiant energy of the star. The concept was first described by [Watson and Lovelock (1983)](https://doi.org/10.1111%2Fj.1600-0889.1983.tb00031.x) [1]

[1] Watson, A.J.; J.E. Lovelock (1983). "Biological homeostasis of the global environment: the parable of Daisyworld". Tellus B. 35 (4): 286–9.

Interactive daisyworld originated during my time as a TA for the MetK: Climate Change course at Goethe-University. It came to be because I wanted to give students an opportunity to see understand the hypothetical idea by playing around with the parameters, but without necessarily understanding code.

The code is based on an Fortran implementation by [Kirsten Menking](https://serc.carleton.edu/quantskills/activities/daisyworld_lab.html) and Joel Dashnaw.
Dept. of Geology and  Geography, Vassar College, Poughkeepsie, NY  12604, July 21, 2003.

## Getting started with the repository

There are different versions of the scripts, depending on the use-case. 

*  `daisyworld.py`
*  `interactive_daisyworld.py`
*  `interactive_daisyworld.ipynb`



To simulate a single experiment, run `python3 daisyworld.py` from the main folder. This will result in a matplotlib plot with the default parameters as shown in problem set 7.

The interactive version of the script either uses [pythonista](https://omz-software.com/pythonista/), an iOS app by Ole Zorn to provide a GUI, or [jupyter-notebook](https://jupyter.org) and widgets to play around with the plant parameters.

The Pythonista version enabled me to hand around an iPad and let students work in groups to figure out how the parameters where influencing the hypothetical world.


### Prerequisites

* [numpy](http://www.numpy.org/)
* [matplotlib](https://matplotlib.org/)
* [jupyter-notebook](https://jupyter.org)
* [pythonista](https://omz-software.com/pythonista/)


### Directory structure

```bash
.
├── interactive_daisyworld.ipynb   # jupyter notebook
├── interactive_daisyworld.pyui    # GUI for pythonista
├── interactive_daisyworld.py      # script for pythonista
├── daisyworld.py                  # simple python script
├── README.md
├── LICENSE
└── imgs                           # image folder for the readme
    ├── phone1.png
    ├── phone2.png
    └── phone3.png
```

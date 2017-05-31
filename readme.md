# SORN

This code contains the necessary files to simulate the SORN model and perform the avalanche analysis as described in Del Papa et al. 2017. It is based on the previous implementation by [Hartmann et al.](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004640), which can be found [here](https://github.com/chrhartm/SORN).

## Experiments

To run simple SORN experiments, proceed in the same way as the original SORN code: navigate to the *common* folder and run `python test_single.py your_parameter_file`. The simulations store all relevant data in a `backup` folder, which you may change according to your needs.

The avalanche scripts should be run independently of the SORN experiments. They read the simulations results, perform the relevant analysis and plot the avalanche events distributions. Remember to update the simulations folder accordingly (usually named `exper_path` in the scripts).

## Files and data for Del Papa et al. 2017

The experiments and parameter files described in the paper can be found in the `delpapa` folder. The subfolder `avalanches` contains the exact code to reproduce each of the figures, as instructed below. Note that a few parameters must be change for each of the experiments (for example, in order to reproduce all the simulations for different network sizes from Fig. 2, you should run simulations for each size independently!).

* Fig. 1
Parameters and Experiments: `param_Zheng2013.py`, `experiment_Zheng2013.py`
Figure: `avalanches/Fig1.py`

* Fig. 2
Parameters and Experiments: `param_Zheng2013.py`, `experiment_Zheng2013.py`
Figure: `avalanches/Fig2.py`

* Fig. 3
Parameters and Experiments: `param_Zheng2013.py`, `experiment_Zheng2013.py`
Figure: `avalanches/Fig3.py`

* Fig. 4
Parameters and Experiments: `param_FrozenPlasticity.py`, `experiment_FrozenPlasticity.py`
Figure: `avalanches/Fig4.py`

* Fig. 5
Parameters and Experiments: `param_Zheng2013.py`, `experiment_Zheng2013.py`
Figure: `avalanches/Fig5.py`

* Fig. 6
Parameters and Experiments: `param_ExtraInput.py`, `experiment_ExtraInput.py`
Figure: `avalanches/Fig6.py`

* Fig. 7
Parameters and Experiments: `param_CountingTask.py`, `experiment_CountingTask.py` and `param_RandomTask.py`, `experiment_RandomTask.py`
Figure: `avalanches/Fig7.py`


## Dependencies

This code relies on the [powerlaw](https://pypi.python.org/pypi/powerlaw) python package to fit the power-law distributions of avalanche sizes and durations.

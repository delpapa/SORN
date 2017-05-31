# SORN

This code contains the necessary files to simulate the SORN model and perform the avalanche analysis as described in Del Papa et al. 2017. It is based on the previous implementation by [Hartmann et al.](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004640), which can be found [here](https://github.com/chrhartm/SORN).

## Experiments

To run simple SORN experiments, proceed in the same way as the original SORN code: navigate to the *common* folder and run `python test_single.py your_parameter_file`. The simulations store all relevant data in a `backup` folder, which you may change according to your needs.

The avalanche scripts should be run independently of the SORN experiments. They read the simulations results, perform the relevant analysis and plot the avalanche events distributions. Remember to update the simulations folder accordingly (usually named `exper_path` in the scripts).

## Files and data for Del Papa et al. 2017

The experiments and parameter files describe in the paper can be found in the `delpapa` folder. The subfolder `avalanches` contains the exact code to reproduce each of the figures, as instructed below.

## Dependencies

This code relies on the [powerlaw](https://pypi.python.org/pypi/powerlaw) python package to fit the power-law distributions of avalanche sizes and durations.

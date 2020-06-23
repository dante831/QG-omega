# QG-omega inversion

This repository accompanies the paper Li and O'Gorman, 2020, "Response of Vertical Velocities in Extratropical Precipitation Extremes to Climate Change". It includes all data used to plot figures in the main paper, and some of the figures in the supporting information. For a copy of the paper, please see this [link](https://pog.mit.edu/src/li_omega_equation_precipitation_extremes_2020.pdf), and the supporting information can be seen at this [link](https://pog.mit.edu/src/li_omega_equation_precipitation_extremes_supp_2020.pdf). 

## Prerequisites

You will need at least Matlab version 2016b for the code to function properly. For the best result, please update to the newest version of Matlab, which is 2019b by the time when this note is written. 

## The main_plot routine

For the full analysis, open a matlab session and directly run main_plot.m in the 'QG-omega' folder: 
```
main_plot
```

## Inverting QG-omega equation using an example event

Routine 'inversion_example.m' in the 'source' folder performs the QG-omega inversion using the strongly implicit method originally due to Stone, 1968. This Matlab routine first reads in an example event in NetCDF format ('data/example_event.nc '), adjusts the computation box according to the data availability, calculates the terms on the right-hand side of the QG-omega equation, adds smoothing to the fields where necessary, and performs the numeric inversion. 

'inversion_example.m' calls 'synoptic_map.m' to plot Fig. 1 in the paper. 'synoptic_map.m' uses 'cbrewer'[[1]](#1) for its color schemes. 

To run the code, simply open a matlab session in the root folder ('QG-omega') and type in
```
inversion_example
```

## Plotting figures

'figure_2.m', 'figure_4.m', 'figure_5.m', 'figure_7.m', 'figure_8.m', and 'plot_terms.m' are the respective plotting routines for the figures in our main paper and the supplementary material. 

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Reference
<a id="1">[1]</a> Charles (2020). cbrewer : colorbrewer schemes for Matlab (https://www.mathworks.com/matlabcentral/fileexchange/34087-cbrewer-colorbrewer-schemes-for-matlab), MATLAB Central File Exchange. Retrieved May 25, 2020.

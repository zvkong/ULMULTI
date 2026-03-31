Files
- functions_indicator_nongroup_revised.R
- functions_indicator_group_revised.R
- empirical_indicator_nongroup_revised.R
- empirical_indicator_group_revised.R

Run
1. Put IL21.csv in the working directory.
2. Run empirical_indicator_nongroup_revised.R for the unit-level binomial version.
3. Run empirical_indicator_group_revised.R for the grouped binomial version.

Notes
- The scripts are self-contained and load the required packages directly.
- The grouped version keeps the Gaussian block at the unit level and collapses only the binomial block.
- These files follow the uploaded reference implementation, including tau_2 in the binomial shared effect.

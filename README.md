# ctc-cytometry-analyis-v1
This is software developed in Matlab 2015a

The main function is analyzeExperimentalPerformance.m
This function will segment the data into training-testing subsets

Depending on bLoad=0,1, it will either load previously saved regressions or perform regressions and save them for each of the training testing subsets. Perforing regressions with parfor and 8 polled threads takes about 1.5 days. Using loaded regressions will take ~20 minutes to complte. 

Functions performing regressiosn are selectAndRegress2.m and feature_selection.m . 


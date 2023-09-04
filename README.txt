# Mixture Hidden Markov Model for Patient Clustering

This repository contains R and Stan code for implementing a Mixture Hidden Markov Model for Patient Clustering. The model is designed to analyze sequential data with hidden states and mixture components.

## Overview

Our model is designed for clustering patient trajectories.

## Repository Contents

casey_master.csv: This contains the dataset that can be used to test and run the model.
Preprocsing.R: The files for preprocessing the dataset. 
Preprocessed.MHMM.Data.Rdata. Data (Rdata file) with prepcoessed data. 
StanCode/: The Stan code for implementing the Mixture HMM is located here.
Analysis.R: The interface to run the Stan code for the Mixture HMM model.

## Getting Started

Clone this repository to your local machine.

Ensure you have the required dependencies installed, including R, Stan, and any necessary R packages.

## Dependencies

R (>= 3.6.0)
Stan (>= 2.18.1)
RStan package (for interfacing R and Stan)

Additional R packages as specified in the R scripts

## License

This project is licensed under the MIT License.

Feel free to use, modify, and distribute this code as needed.
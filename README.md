# integrated
Integrated models of population dynamics in R

Copyright &copy; 2018, Jian Yen

*****

## License details
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

*****

## Overview
Integrated model to estimate demographic vital rates from data on individual growth trajectories and population size distributions.

Maintainer: Jian Yen (jdl.yen@gmail.com)

Updated: 12 March 2018

## About
This project uses integrated Bayesian models to estimate demographic vital rates from multiple data sources. Models are fitted in R using the [greta R package](https://github.com/greta-dev). The current models use a stage-, age-, or un-structured population matrix model with eight possible classes in the example given here.

These models use two types of data: population abundance surveys with individual size data and individual growth trajectories. Simulated data are provided here, intended to match broadly the characteristics of a large-bodied, predatory freshwater fish species.

## Usage
The main.R script contains the full model run, calling on helper functions split into several different scripts. The primary helper function is estimate_mpm.R in the mpm_est.R script. Model runs with 1000 iterations take up to an hour on a MacBook Air.

A helper script, install_packages.R, is called within main.R and will install all required packages from GitHub. Issues with installing the greta package can be resolved by following instructions at https://github.com/greta-dev.


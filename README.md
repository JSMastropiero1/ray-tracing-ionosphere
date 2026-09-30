# Ray Tracing of Radio Waves in the Ionosphere

A numerical physics project exploring the propagation of radio-frequency waves through the ionosphere using simplified electron-density and refractive-index models.

The project was originally developed in 2022 as part of our Master's-levels for the subject Computational and Numerical Methods. It combines physical modeling, numerical ray tracing, parameter sweeps, and visualization in Python.

## Project overview

The aim of the project is to study how radio-wave trajectories depend on the properties of the ionospheric plasma. The code models a spatially varying electron density, converts it into a refractive-index profile, and propagates rays through the resulting medium.

The repository explores, among other cases:

- day and night ionospheric electron-density profiles;
- linear and quadratic approximations for the density distribution;
- radio-wave trajectories for different launch angles and frequencies;
- refraction and reflection as the refractive index changes with altitude;
- propagation range as a function of frequency and incidence angle;
- perturbations of the electron-density profile;
- Cartesian and polar-coordinate implementations of the ray trajectory.


The numerical implementation uses a discretized refractive-index profile and applies geometrical-optics/Snell-law considerations as the ray propagates through the inhomogeneous medium.

## Repository contents

The project evolved iteratively, so the repository contains both core routines and exploratory notebooks developed while testing different numerical approaches.

Subdirectory "notebooks" shows a cleaner version of the developed code, while subdirectory "archive" keeps the original files.

Several similarly named notebooks are retained intentionally because they document the development and testing process rather than representing separate finished software releases.


## Project status

This repository preserves the original exploratory character of the 2022 project. It is intended primarily as a record of the numerical modeling work and thwe development process, rather than as a production-ready ray-tracing package.

The code may therefore contain duplicated exploratory notebooks, hard-coded model parameters, and implementation choices that would be refactored in a modern software project. These have largely been retained to preserve the original work and its development history.


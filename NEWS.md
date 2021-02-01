# JMbayes2 0.1.4

## Major
* Added two vignettes (available on the [website of the package](https://drizopoulos.github.io/JMbayes2/)) to showcase joint models with competing risks and joint models with non-Gaussian longitudinal outcomes.


# JMbayes2 0.1.3

## Minor
* Used `parallel::clusterSetRNGStream()` in `jm_fit()` for distributing the seed in the workers.
* Changed the default position of the knots for the B-spline approximation of the log baseline hazard.


# JMbayes2 0.1.2

## Minor
* Changed calls to `floor()` in the C++ code.


# JMbayes2 0.1.0

## General
* First version of the package.


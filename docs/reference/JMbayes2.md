# Extended Joint Models for Longitudinal and Time-to-Event Data

Fit joint models for longitudinal and time-to-event data under the
Bayesian approach. Multiple longitudinal outcomes of mixed type
(continuous/categorical) and multiple event times (competing risks and
multi-state processes) are accommodated.

## Details

|          |            |
|----------|------------|
| Package: | JMbayes2   |
| Type:    | Package    |
| Version: | 0.6-0      |
| Date:    | 2026-01-30 |
| License: | GPL (\>=3) |

This package fits joint models for longitudinal and time-to-event data.
It can accommodate multiple longitudinal outcomes of different type
(e.g., continuous, dichotomous, ordinal, counts), and assuming different
distributions, i.e., Gaussian, Student's-t, Gamma, Beta, unit Lindley,
censored Normal, Binomial, Poisson, Negative Binomial, and
Beta-Binomial. For the event time process, right, left and interval
censored data can be handled, while competing risks and multi-sate
processes are also covered.

JMbayes2 fits joint models using Markov chain Monte Carlo algorithms
implemented in C++. The package also offers several utility functions
that can extract useful information from fitted joint models. The most
important of those are included in the **See also** Section below.

## Author

Dimitris Rizopoulos, Grigorios Papageorgiou, Pedro Miranda-Afonso

Maintainer: Dimitris Rizopoulos \<d.rizopoulos@erasmusmc.nl\>

## References

Rizopoulos, D. (2012). Joint Models for Longitudinal and Time-to-Event
Data With Applications in R. Boca Raton: Chapman & Hall/CRC.

## See also

[`jm`](https://drizopoulos.github.io/JMbayes2/reference/jm.md),
[`methods.jm`](https://drizopoulos.github.io/JMbayes2/reference/methods.md),
[`coda_methods.jm`](https://drizopoulos.github.io/JMbayes2/reference/coda_methods.md)

# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

# Tutorial to use SciPy library with Beta distribution:
# https://towardsdatascience.com/probability-distributions-with-pythons-scipy-3da89bf60565
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time
import math
from scipy import stats

from scipy.stats import (
    norm, beta, expon, gamma, genextreme, logistic, lognorm, triang, uniform, fatiguelife,            
    gengamma, gennorm, dweibull, dgamma, gumbel_r, powernorm, rayleigh, weibull_max, weibull_min, 
    laplace, alpha, genexpon, bradford, betaprime, burr, fisk, genpareto, hypsecant, 
    halfnorm, halflogistic, invgauss, invgamma, levy, loglaplace, loggamma, maxwell, 
    mielke, ncx2, ncf, nct, nakagami, pareto, lomax, powerlognorm, powerlaw, rice, 
    semicircular, trapezoid, rice, invweibull, foldnorm, foldcauchy, cosine, exponpow, 
    exponweib, wald, wrapcauchy, truncexpon, truncnorm, t, rdist
    )

a, b = 2, 8

x = beta.rvs(a, b, size=1000)

# Histogram
fig, ax = plt.subplots(1, 1)
ax.hist(x, density=True, histtype='stepfilled', alpha=0.2)
ax.legend(loc='best', frameon=False)
plt.show()

# Store parameters in rv
rv = beta(a,b)
moments = rv.stats("mvsk")
_ = [print(f'{v:.3f}') for v in moments]


# plot the pdf
x = np.linspace(rv.ppf(0.01),
                rv.ppf(0.99), 100)

fig, ax = plt.subplots(1, 1)
ax.plot(x, rv.pdf(x)/sum(rv.pdf(x)), 'r-', lw=5, alpha=0.6, label='beta pdf')"""



# https://vitalflux.com/beta-distribution-explained-with-python-examples/
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import beta
#
# Set the shape paremeters
#
a, b = 2, 6
#
# Generate the value between
#
x = np.linspace(beta.ppf(0.01, a, b),beta.ppf(0.99, a, b), 100)
#
# Plot the beta distribution
#
plt.figure(figsize=(7,7))
#plt.xlim(0.7, 1)
const_norm=sum(beta.pdf(x, a, b))
plt.plot(x, beta.pdf(x, a, b)/const_norm, 'r-')
plt.title('Beta Distribution', fontsize='15')
plt.xlabel('Values of Random Variable X (0, 1)', fontsize='15')
plt.ylabel('Probability', fontsize='15')
plt.show()

"""

# https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.skewnorm.html

import numpy as np
from scipy.stats import skewnorm
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 1)

a = 20

mean, var, skew, kurt = skewnorm.stats(a, moments='mvsk')


x = np.linspace(skewnorm.ppf(0.001, a), skewnorm.ppf(0.999, a), 100)

ax.plot(x, skewnorm.pdf(x, a),'r-', lw=5, alpha=0.6, label='skewnorm pdf')
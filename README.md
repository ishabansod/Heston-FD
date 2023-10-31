# Finite difference methods for stochastic volatility modelling

## Volatility smile modelling
The 1973 paper of F. Black and M. Scholes revolutionized the theory of derivative pricing, paved the way to the concept of risk-neutral valuation and brought forth the famous Black-Scholes formula. A drawback of their framework is that it cannot explain a phenomenon that is known as the volatility-smile or volatility-skew.
Under the Black-Scholes model, assets are modelled as geometric Brownian Motions and are characterized by two main parameters: a constant drift and volatility. The latter cannot be directly measured and is therefore often implied by inverting the option-valuation formula and inserting option-prices observed in the market.
Should the Black-Scholes framework be realistic, then the implied volatility would be independent of the strike price and maturity of an option contract. In reality however, the implied volatility as a function of the strike price exhibits a convex (a smile) or even asymmetric (a skew) shape and changes with contract duration yielding a non-constant volatility surface, see Figure 1. The Heston model [[2]](#2) is a popular alternative to Black-Scholes, which addresses this short-coming by modelling the volatility as a stochastic process.

## Project description
For this project you will familiarize yourself with the Heston stochastic volatility model [[2]](#2). A modern outlook on the computational aspects of this model can be found in [[1]](#1).
The first objective is to implement a European option pricer under the Heston model. There exist several versions of the semi-analytical valuation formula, which typically require numerical integration. You may choose any programming language to set-up your implementations.
The main objective is to construct a finite-difference (FD) framework to price options under the Heston model.
For this you will need to derive and discretise the Heston PDE. A nice introduction to finite-difference methods can be found in [[3]](#3). Experiment with the degrees of freedom that FD methods offer. Test the accuracy, convergence and stability by comparing to the semi-analytical pricer.

## References
<a id="1">[1]</a> 
Jim Gatheral. The volatility surface: a practitioner’s guide. John Wiley & Sons, 2011.

<a id="2">[2]</a> 
Steven L Heston. “A closed-form solution for options with stochastic volatility with applications to bond
and currency options”. In: The review of financial studies 6.2 (1993), pp. 327–343.

<a id="3">[3]</a> 
R¨udiger Seydel. Tools for computational finance. Vol. 3. Springer, 2006.
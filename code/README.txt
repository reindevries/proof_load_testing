This C++ code belongs to the PhD thesis of Rein de Vries titled "Structural reliability updating 
through proof load testing – A Bayesian methodology applied to reinforced concrete road bridges and 
viaducts". The data is organised following the structure of the dissertation, which can be found via
the TUD repository: https://repository.tudelft.nl/

The dissertation followed from the PhD project "Probabilistic substantiation of proof load testing" 
executed at the Concrete Structures department of the Delft University of Technology (TUD) from 1 
September 2020 to 1 September 2025. The project number is CS3B09, and the Directorate-General for 
Public Works and Water Management (Rijkswaterstaat, RWS) has financed the work. This project was a 
collaboration between TUD, TNO and RWS.

This code requires the following compiler, libraries and editor:
• GCC C++ 10.2 or later, 64-bit version installed in %USERPROFILE%\msys64\mingw64
• Boost C++ library 1.65 or later, located in %USERPROFILE%\cpp_include (headers only)
• Eigen C++ library 3.4 or later, located in %USERPROFILE%\cpp_include
• Visual Studio Code 1.103.2 or later, with extensions: C/C++, C/C++ Themes, F5 Anything

Description of files:
• .vscode/c_cpp_properties.json
  .vscode/launch.json
  .vscode/tasks.json
  These are the project configuration files for Visual Studio Code. The project is set up such that
  when running (pressing F5) the currently displayed source file is compiled and executed. Debug and
  Release configurations are provided, with the executables following from the latter are
  significantly faster. This speed is absolutely required for the Monte Carlo simulations, which
  would otherwise take a lot of time.

• rdv/distributions.hpp
  rdv/erf_inv.hpp
  ...
  Together, these files contain the header-only "rdv" probabilistic library that is used in the rest 
  of the code, containing the implementation of distributions, mathematical functions, etc.

• chapter_2.cpp
  Contains the Monte Carlo simulation procedure to produce graphs displaying the evolution of
  reliability over time. Importance sampling is used on variables f_y, C_0Q and θ_E to increase the
  accuracy and ultimately reduce the run time. The output is displayed in Figures 2.2 and 2.3.

• chapter_3.cpp
  This code iteratively calculates the target proof load factors for CC2 and CC3 in bending and
  shear. These are the factors provided in data file Chapter 3 Factors CC2 and CC3.csv and are
  plotted in Figure 3.3.

• chapter_4_prior_sensitivity.cpp
  Using this code, the sensitivity to various choices for the prior distribution of the resistance
  is explored. The distribution is changed by uncommenting the relevant rv_R_hat definition. The
  result of the calculation is the posterior distribution of the resistance, assuming survival of
  the proof load effect. The results of these calculations are displayed in Figures 4.2 and 4.3.

• chapter_4_time-dependent.cpp
  This code is quite similar to chapter_2.cpp as it also performs a time-dependent reliability
  analysis. However, in this case a low-informative prior distribution is assumed for the resistance
  with its mean value directly calculated from the mean values of the load variables. The output of
  the analysis is provided in Figure 4.4.

• chapter_5_case_study_1.cpp
  In this code the reliability analysis of case study 1, where laboratory data is used to predict
  the resistance, is performed. Several target load levels are used to simulate the incremental
  load application. For each step a Markov-chain Monte Carlo (MCMC) simulation is performed by which
  the reliability during and after a successful load step is calculated. The output of this code is
  provided in Table 5.1.

• chapter_5_case_study_2_sectional_analysis.cpp
  Since for the second case study no laboratory data is available, a numerical simulation using
  Latin Hypercube Sampling (LHS) is performed. The output of this code has resulted in the data of
  Chapter 5 Sectional analysis.csv and is displayed in Figure 5.6.

• chapter_5_case_study_2_reliability_updating.cpp
  This code is similar to chapter_5_case_study_1.cpp, but in this case the resistance is predicted
  using the numerically simulated data. The resulting reliability indices are provided in Table 5.2.

• chapter_6_load_correlation_decay.cpp
  In this code a Monte Carlo simulation is performed to determine the decay of the correlation
  between load effects when the reference period increases. Gumbel distributed random variables with
  given initial correlation and block sizes are defined using vectors. The output is provided in
  Figure 6.5.
  
• chapter_6_case_study_1_bending.cpp
  chapter_6_case_study_1_shear.cpp
  These codes calculate transfer factors for case study 1, considering bending and shear. The flag
  'correlated_Q' switches between the correlated and uncorrelated situations. Using interpolation, 
  the target load corresponding to beta = 4 may be determined. The ratio between the load required 
  in each configuration and the n = N case results in the transfer factors. The results for bending 
  and shear are provided in Tables 6.3 and 6.4.

• chapter_6_case_study_2.cpp
  Using this code, the reliability of the continuous bridge with two spans can be calculated, given
  two load testing strategies and two load test vehicle configurations (tipper truck and semi-low
  trailer). If the axle loads are set to zero, the code can also be used to calculate the in-service
  proven strength from one year of traffic loading. The calculated reliability indices are provided 
  in Table 6.5.

Contacts:
• Rein de Vries, rein.devries MONKEY TAIL tno.nl
• Eva Lantsogt, e.o.l.lantsoght MONKEY TAIL tudelft.nl

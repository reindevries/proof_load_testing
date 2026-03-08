This dataset belongs to the PhD thesis of Rein de Vries titled "Structural reliability updating 
through proof load testing – A Bayesian methodology applied to reinforced concrete road bridges and 
viaducts". The data is organised following the structure of the dissertation, which can be found via
the TUD repository: https://repository.tudelft.nl/

The dissertation followed from the PhD project "Probabilistic substantiation of proof load testing" 
executed at the Concrete Structures department of the Delft University of Technology (TUD) from 1 
September 2020 to 1 September 2025. The project number is CS3B09, and the Directorate-General for 
Public Works and Water Management (Rijkswaterstaat, RWS) has financed the work. This project was a 
collaboration between TUD, TNO and RWS.

Description of files:
• Chapter 3 Fits bending and shear SI.csv
  Chapter 3 Fits bending and shear US.csv
  These files contain the Gumbel distribution fits, where column "x" refers to the span length in m
  (SI) and in ft (US). The remaining columns provide the mean value for each road and the 
  corresponding coefficient of variation (COV). The distribution parameters referring to shear force 
  are indicated with suffix "_Shear". The mean values are provided in MNm (SI) and in kip·ft (US).
  The mean and coefficient of variation plots are collected in Figure 3.2.

• Chapter 3 Factors CC2 and CC3.csv
  This file contains the target proof load factors that have been calculated using the lower-bound 
  method in Chapter 3. The first column "x" refers to the span length in m. In the remainder of 
  the columns containing the factors, CC refers to the consequence class, M indicates moment, and V 
  indicates shear. The proof load factors are plotted in Figure 3.3.

• Chapter 5 Beam test data H121.csv
  Chapter 5 Beam test data H401.csv
  Chapter 5 Beam test data H403.csv
  Chapter 5 Beam test data H404.csv
  Chapter 5 Beam test data H602.csv
  These files contain the data that followed from the DIC measurements on the lab tests of the
  beams tested in the TUD laboratory. The first row contains the applied load (P) in kN, and 
  each corresponding column contains the calculated nominal crack width in mm using a virtual LVDT
  centred around the location that is indicated using the fraction M/(Vd) in the first column.

• Chapter 5 Beam test data (incl post-processing).xlsx
  In this file, the same data are contained as described above, including their post-processing
  calculating the maximum that has occurred during all previous load steps and all considered
  positions, and utilising a weighing factor based on their locations. Note that the weighted
  crack width is still expressed in mm, but it does not directly correspond to measured crack widths
  any more. The results are also plotted in Figure 5.3.

• Chapter 5 Sectional analysis.csv
  This file contains the section analysis results from the LHS procedure in which the mechanical 
  properties of a cross-section are varied. The first column contains the reinforcement steel 
  strains (no unit) and the second row provides the yield moment of each realisation. In total 100 
  virtual experiments were performed, resulting in the curves indicated with M_i where i is the 
  experiment or realisation number. The moments are provided in kNm.

• Chapter 5 Sectional analysis (incl post-processing).xlsx
  Just like with the laboratory experiments, a second file is provided that contains the data and
  also includes the post-processing. No additional weighing is performed; the ratio is directly
  calculated as the current load effect divided by the yield moment.

• Chapter 6 Repeat tests with shear reinforcement.csv
  Chapter 6 Repeat tests without shear reinforcement.csv
  These files contain the test results that were selected from the large ACI-DAfStb shear da‍ta‍ba‍se,
  separated for beams with and without shear reinforcement. The mean compressive strength is 
  provided in MPa, and the value of the test results is in kN. The coefficient of variation (COV) 
  was calculated from the test results. A copy of the complete shear database may be requested via:
  https://dafstb.de/aci-dafstb.html

As described in Chapter 3 of the dissertation, weigh-in-motion (WIM) data have been used to 
characterise the load effect from trucks statistically. This dataset may be obtained by contacting 
Martin Nijgh of RWS, specifying the intended usage of the data.

Contacts:
• Rein de Vries, rein.devries MONKEY TAIL tno.nl
• Eva Lantsogt, e.o.l.lantsoght MONKEY TAIL tudelft.nl
• Martin Nijgh, martin.nijgh MONKEY TAIL rws.nl

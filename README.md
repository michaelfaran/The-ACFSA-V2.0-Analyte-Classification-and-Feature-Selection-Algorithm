<p align="center"> <img src="TOC_Image.png" > </p>
About:

This repository contains the code framework for ACFSA V2.0: Analyte Classification and Feature Selection Algorithm.

It serves as an open source with a graphical user interface (GUI) for anyone interested in choosing a minimal sensor set out of analyte-sensor screening data, that preserves both accuracy and minimal sensor number- a key purpose in experiment design.

The algorithm, with some key examples, appears in the upcoming paper "Rational Selection of Minimal Sensor Arrays for Analyte Fingerprinting" by Michael Faran, Minyeong Yoon, Soo-Yeon Cho, and Gili Bisker [1]. The paper is now under revision.

The paper introduces ACFSA V2.0 as a practical tool for experimental design, with the following merits:

Optimized sensor sets – Takes a full screening table and returns the smallest sensor subset that preserves high classification accuracy.

Cost and complexity reduction – Fewer sensors mean lower experimental cost, simpler hardware, and faster move from screening to deployable panels.

Robust to real lab data – Designed for noisy or limited datasets and does not rely on neural networks that may overfit.

Transparent decision-making – Visual decision maps show why chosen sensors work and where misclassifications may occur.

Platform-agnostic – Works with any cross-reactive array and with both continuous and discrete readouts.

No coding required – GUI workflow: load data, set accuracy targets, and export publication-ready figures.

Statistical stability features introduced in the revised version

ACFSA V2.0 includes several statistical stabilization mechanisms introduced during peer-review revision:

Adaptive PCA dimensionality.
When the first principal component explains 95% or more of the total variance, the algorithm automatically performs classification in the PC1 space only, instead of the PC1–PC2 plane. This prevents unstable covariance estimation when the dataset is effectively one-dimensional.

Regularized QDA classifier.
When the QDA classifier is selected, covariance matrices are stabilized using Ledoit–Wolf shrinkage toward the pooled covariance matrix, followed by ridge eigenvalue regularization, improving robustness for small sample sizes and highly correlated sensor arrays.

Consistent statistical treatment.
The same covariance regularization logic is applied both during classification and during the feature-selection step to maintain numerical stability throughout the sensor elimination process.

Proven peformance:

Builds on prior ACFSA – Extends the original ACFSA framework [2].

Validated on multiple datasets – Demonstrated on:

Metal-ion SWCNT dataset [1]

Artificial datasets based on [2]

DNA-SWCNT sensors for sweat-related analytes [3]

DNA-SWCNT sensors for urine analytes [4]

Supports alternative activation modes – Benchmarked across different ACFSA V2.0 activation settings.

This code was written by Michael Faran, 16/11/2025, in MATLAB.
For any questions or inquiries, please email michaelfaran[at]gmail.com.
The code will be continually updated upon request.

Software Requirements:

MATLAB: 24.1 (R2024a) and later versions.

MATLAB Mapping Toolbox 24.1

MATLAB Statistics and Machine Learning Toolbox 24.1

The code uses Excel: Version 2508 (Build 19127.20264), but will likely work on other versions.

The code was tested on Windows: Windows 11, Version 22H2 (OS Build 22621.3880), but will likely work on other versions or other operating systems.

Installation:

Ensure the software requirements above are verified.

Clone or download the repository to a chosen folder.

Activation:

Activate the main script "Main.m" from the chosen repository folder, and the GUI main screen should appear.

Put the inputs of the number of sensors, analytes, and the number of samples per analyte. The code requires at least 3 samples per analyte.

Then, click Create a template & open in Excel- a template to input each sensor measurement per the same sample appears based on the input provided.

Fill the template values with sensor responses after post-processing if needed (For SWCNTs fluorescence datasets, normalized responses were used before). Suppose sensor measurements are taken from different samples of the same analyte. In that case, it is possible to concatenate them into a single analyte measurement row, with the risk of inducing spurious artificial correlations and losing realistic in-sample sensor reading correlations.

Save the Excel file and load it using the GUI "import from XLSX", validate, and approve your choices.

Click "Next: Configure & RUN ACFSA V2.0 ".

Fill up the different fields according to Configure & Run inputs below, activate as default, and check [1] for more details.

Click "Run ACFSA"-All other figures will be closed now, including the GUI activation window, and the algorithm will run. Do not press newly created figures when MATLAB runs, as it might distort the activation output figures.

Enjoy the new selected minimal sensor data set, appearing in the repo/results/CONFIG folder name (see "output.txt" for main output and others in "supporting output")

Activation notes:

We suggest first activate the ACFSA V2.0 on the default dataset found in "[Repo local address]/Open Source/examples/measurement_data_default_data_set.xlsx". If any errors happen during this run, please reach out to michaelfaran [at] gmail.com for help.

If you encounter an error during the run of your own data set, please see the current limitations of the scheme below. If this does not solve the issue, please contact michaelfaran[at]gmail.com for assistance.

Scheme Current Limitations:

The code assumes the same number of measurements per analyte, and at least 3 measurements for each analyte are required.

At least two analytes are required.

The code assumes the input response is already after post-processing. That might mean different implications for different sensor sets. As an example, a SWCNT normalized fluorescence intensity response can be obtained by integrating the emission spectrum and normalizing by the integral measured at zero analyte concentration

The code is limited to plotting a maximum of six analytes in the data set, due to graphical constraints.

Some graphical mismatches in the output files can arise due to long analyte names.

ACFSA V2.0 assumes that each analyte measurement is drawn from Gaussian statistics, and classification error is calculated accordingly.

When the variance explained by the first principal component exceeds 95%, the classifier automatically switches to a 1-dimensional PC1 classification, avoiding unstable covariance estimation when the data are effectively one-dimensional.

References:

[1] Faran, Michael, et al. “The ACFSA V2.0: Analyte Classification and Feature Selection Algorithm.” Manuscript submitted (2025).
[2] Petresky, Gabriel, et al. "Metal-Ion Optical Fingerprinting Sensor Selection via an Analyte Classification and Feature Selection Algorithm." Analytical Chemistry 97.16 (2025): 8821-8832.
[3] Lee, Yeon Soo, et al. "Spatiotemporal molecular tracing of ultralow-volume biofluids via a soft skin-adaptive optical monolithic patch sensor." Nature Communications 16.1 (2025): 3272.‏
[4] Yoon, Minyeong, et al. "Enzyme-free optical detection of uric acid using corona phase molecular recognition in near-infrared fluorescent single-walled carbon nanotubes." Nanoscale 17.17 (2025): 10652-10662.‏

[5] Ledoit, Olivier, and Michael Wolf. "A Well-Conditioned Estimator for Large-Dimensional Covariance Matrices." Journal of Multivariate Analysis 88.2 (2004): 365–411.

[6] Friedman, Jerome H. "Regularized Discriminant Analysis." Journal of the American Statistical Association 84.405 (1989): 165–175.

The code uses adjusted versions of the following MATLAB codes:
plot_elipse.m, taken from:
[7] https://www.mathworks.com/matlabcentral/fileexchange/116610-plot-ellipse-on-scattered-2d-data?s_tid=prof_contriblnk

[8] Ohad Gal (2025). fit_ellipse (https://www.mathworks.com/matlabcentral/fileexchange/3215-fit_ellipse
), MATLAB Central File Exchange. Retrieved November 15, 2025.
[9] Jakob Sievers (2025). VoronoiLimit(varargin) (https://www.mathworks.com/matlabcentral/fileexchange/34428-voronoilimit-varargin
), MATLAB Central File Exchange. Retrieved November 15, 2025.
[10] qqffssxx (2025). Rand and Adjusted Rand Index Calculator for Cluster Analysis (https://www.mathworks.com/matlabcentral/fileexchange/130779-rand-and-adjusted-rand-index-calculator-for-cluster-analysis
), MATLAB Central File Exchange. Retrieved November 15, 2025.
&lt;p align="center"&gt;
  &lt;img src="TOC_Image.png" &gt;
&lt;/p&gt;

## About:
This repository contains the code framework for ACFSA V2.0: Analyte Classification and Feature Selection Algorithm.

It serves as an open source with a graphical user interface (GUI) for anyone interested in choosing a minimal sensor set out of analyte-sensor screening data, that preserves both accuracy and minimal sensor number- a key purpose in experiment design.

The algorithm, with some key examples, appears in the upcoming paper  "The ACFSA V2.0: Analyte Classification and Feature Selection Algorithm" by Michael Faran, Minyeong Yoon, Soo-Yeon Cho, and Gili Bisker [1]. The paper is now under revision.

The paper introduces ACFSA V2.0 as a practical tool for experimental design, with the following merits:

* **Optimized sensor sets** – Takes a full screening table and returns the smallest sensor subset that preserves high classification accuracy.
* **Cost and complexity reduction** – Fewer sensors mean lower experimental cost, simpler hardware, and faster move from screening to deployable panels.
* **Robust to real lab data** – Designed for noisy or limited datasets and does *not* rely on neural networks that may overfit.
* **Transparent decision-making** – Visual decision maps show why chosen sensors work and where misclassifications may occur.
* **Platform-agnostic** – Works with any cross-reactive array and with both continuous and discrete readouts.
* **No coding required** – GUI workflow: load data, set accuracy targets, and export publication-ready figures.

Proven peformance:

- **Builds on prior ACFSA** – Extends the original ACFSA framework [2].  
- **Validated on multiple datasets** – Demonstrated on:  
  - Metal-ion SWCNT dataset [1]  
  - Artificial datasets based on [2]  
  - DNA-SWCNT sensors for sweat-related analytes [3]  
  - DNA-SWCNT sensors for urine analytes [4]  
- **Supports alternative activation modes** – Benchmarked across different ACFSA V2.0 activation settings.

This code was written by Michael Faran, 16/11/2025, in MATLAB. 
For any questions or inquiries, please email michaelfaran[at]gmail.com. 
The code will be continually updated upon request.

## Software Requirements:
- MATLAB: 24.1  (R2024a) and later versions. 
- MATLAB Mapping Toolbox  24.1 
- MATLAB Statistics and Machine Learning Toolbox  24.1
- The code uses Excel: Version 2508 (Build 19127.20264), but will likely work on other versions.
- The code was tested on Windows: Windows 11, Version 22H2 (OS Build 22621.3880), but will likely work on other versions or other operating systems.

## Installation:
1. Ensure the software requirements above are verified.
2. Clone or download the repository to a chosen folder.

## Activation:
1. Activate the main script "Main.m" from the chosen repository folder, and the GUI main screen should appear.
2. Put the inputs of the number of sensors, analytes, and the number of samples per analyte. The code requires at least 3 samples per analyte.
3. Then, click Create a template & open in Excel- a template to input each sensor measurement per the same sample appears based on the input provided. 
4. Fill the template values with sensor responses after post-processing if needed (For SWCNTs fluorescence datasets, normalized responses were used before). Suppose sensor measurements are taken from different samples of the same analyte. In that case, it is possible to concatenate them into a single analyte measurement row, with the risk of inducing spurious artificial correlations and losing realistic in-sample sensor reading correlations.
5. Save the Excel file and load it using the GUI "import from XLSX", validate, and approve your choices.
6. Click "Next: Configure & RUN ACFSA V2.0 ". 
7. Fill up the different fields according to Configure & Run inputs below, activate as default, and check [1] for more details.
8. Click "Run ACFSA"-All other figures will be closed now, including the GUI activation window, and the algorithm will run. Do not press newly created figures when MATLAB runs, as it might distort the activation output figures. 
9. Enjoy the new selected minimal sensor data set, appearing in the repo/results/CONFIG folder name (see "output.txt" for main output and others in "supporting output")

## Activation notes:
1. We suggest first activate the ACFSA V2.0 on the default dataset found in "[Repo local address]/Open Source/examples/measurement_data_default_data_set.xlsx". If any errors happen during this run, please reach out to michaelfaran [at] gmail.com for help.
2. If you encounter an error during the run of your own data set, please see the current limitations of the scheme below. If this does not solve the issue, please contact michaelfaran[at]gmail.com for assistance.

## Scheme Current Limitations:
1. The code assumes the same number of measurements per analyte, and at least 3 measurements for each analyte are required.
2. At least two analytes are required.
3. The code assumes the input response is already after post-processing. That might mean different implications for different sensor sets. As an example, a SWCNT normalized fluorescence intensity response can be obtained by integrating the emission spectrum and normalizing by the integral measured at zero analyte concentration
4. The code is limited to plotting a maximum of six analytes in the data set, due to graphical constraints. 
5. Some graphical mismatches in the output files can arise due to long analyte names. 
6. ACFSA V2.0 assumes that each analyte measurement is drawn from Gaussian statistics, and classification error is calculated accordingly. 

## User Inputs:

### Step 1 – Template & Import (data/layout inputs)
- Analytes (A) — integer ≥ 1 (default 3)
- Sensors (S) — integer ≥ 1 (default 20)
- Measurements per analyte (M) — integer ≥ 1 (default 3)
- Import XLSX file — chosen via file picker, where the measurements data are inserted (produces DAT, analyte_name_vec, sensor_name_vec, n_measure)

### Step 2 – Configure & Run (algorithmic inputs)
- Dataset title (supertitle) — text (default “ACFSA V2.0 Run”)
- Inflate data to cover statistical uncertainty? — 0/1 (default 1)
- Artificial dataset: STD multiplier — ≥ 0 (default 0; 0 = none)
- Classifier — 0 = QDA, 1 = Voronoi (default 0)
- Weighted feature selection? — 0/1 (default 0)
- Allowed classification error per sensor (%) — Working point 1 — [0,100] (default 3)
- Allowed classification error per sensor (%) — Working point 2 — [0,100] (default 0.5)

### Packed vector passed downstream
- Dataset title 
- input_vec / inputs_mat (1×6) = [inflateFlag, artificialSTDmultiplier, classifierFlag, weightedFSFlag, WP1_pct, WP2_pct]  
→ i.e. [1, 0, 0, 0, 3, 0.5] with defaults.

## Code Outputs

### Output Token:
This run saves results under:
results/CONFIG/
where CONFIG is constructed from your supertitle and the flags below.

- inflateFlag (index 1): Inflate Gaussian STD for uncertainty  
Token: _I (if 1), _NI (if 0)
- artificialSTDmultiplier (index 2): Synthetic stress-test multiplier  
Token: _Nbuff (if 0), otherwise _&lt;val&gt;_buff where &lt;val&gt; is formatted with %.15g and then “.” → “_”.  
Examples: 1.4 → _1_4_buff, 0.001 → _0_001_buff
- classifierFlag (index 3): Decision boundaries  
Token: _QDA (if 0), _Vor (if 1)
- standardFSFlag (index 4): χ²/feature-selection flavor  
Token: _WFS (if 0), _SFS (if 1)
- WP1_pct (index 5): Allowed classification error per sensor (%) — working point 1  
Token: (none) — affects plots/thresholds only
- WP2_pct (index 6): Allowed classification error per sensor (%) — working point 2, this input defines the output PCA and decision lines last presented iteration.

### Outputs:
1. The CONFIG folder name  
CONFIG =
  "Config " + supertitle
  + Inflate_name          (% from inflateFlag → _I / _NI)
  + Buff_name             (% from artificialSTDmultiplier → _Nbuff or _&lt;val&gt;_buff)
  + Stop_name             (% now empty string)
  + chi_name              (% from standardFSFlag → _WFS or _SFS)
  + chi_name2             (% now empty string)
  + Bond_name             (% from classifierFlag → _QDA or _Vor)

   Example with (new) typical defaults  
   Assuming your GUI defaults are: inflateFlag=1, artificialSTDmultiplier=0, classifierFlag=0, standardFSFlag=0, WP1=3, WP2=0.5, and supertitle="ACFSA v2 Run"  
   → *CONFIG = Config ACFSA v2 Run_I_Nbuff_WFS_QDA*

2. Files created inside results/CONFIG/  
(Replace CONFIG below with the actual string.)

#### Main Output
The minimal sensor set for each working point, with its corresponding mean classification error, is printed in "output.txt".  
Example for the default dataset:

*ACFSA V2.0 run summary  
Config: Config ACFSA v2 Run_I_Nbuff_SFS_QDA  
Interpretation: Two working points are shown as vertical lines in Config ACFSA v2 Run_I_Nbuff_SFS_QDA_Summary_Fig.png (dashed = WP1, solid = WP2).  
[WP1 – Alternative] lambda_1 = 0.500%  
  Minimal sensor set size: 2  
  Mean classification error at WP1: 0.0645%  
  Remaining sensors (2): (6,5)_NOX_SWCNT-Gly, (6,5)_OX_SWCNT-Gly  
[WP2 – Default]     lambda_2 = 3.000%  
  Minimal sensor set size: 1  
  Mean classification error at WP2: 2.1347%  
  Remaining sensors (1): (6,5)_NOX_SWCNT-Gly*

#### Supporting Output
- CONFIG_Summary_Fig.(fig|png): This figure outputs the key results of the ACFSA v2 scheme. It shows the first iteration and default working-point results, PCA and decision-boundary plots, alongside the working-point selection panel. If the default working point (WP2) leaves only one sensor, the PCA and decision-boundary plots instead display the two-sensor case.
- Surviving_Sensors_CONFIG_&lt;N&gt;.mat: This .mat file contains the data of the per-iteration state:  
  - names_vec: The remaining sensor names.  
  - sensor_original_number_removed: The sensor number, out of the initial sensor number labels, that was eliminated during iteration N.  
  - ARI: The adjusted Rand index of the classifier versus the ground truth labels.  
  - updated_num_sensors: The current remaining number of sensors.

- First_and_Last_Iterations_PCA.(fig|png)This figure shows the first-iteration and default working-point results of the PCA plot, where WP2 influences which N is shown.
- &lt;Classifier&gt;_Classifier_First_and_Last_Iteration.(fig|png)This figure shows the first-iteration and default working-point results of the decision lines plot, where WP2 influences which N is shown.  
  &lt;Classifier&gt; is:  
  – QDA if classifierFlag = 0  
  – Vor if classifierFlag = 1
- PCA of CONFIG &lt;N&gt; Sensors.(fig|png): This figure shows the PCA plot of iteration N during elimination, e.g., N = 2…#sensors.  
- Surviving_Sensors_&lt;N&gt;_QDA_Classifier.(fig|png) if classifierFlag = 0: This figure shows the QDA classifier plot of iteration N during elimination.
- Surviving_Sensors_&lt;N&gt;_Vor_Classifier.(fig|png) if classifierFlag = 1: This figure shows the Vor classifier plot of iteration N during elimination.
- &lt;supertitle&gt;.mat: Run-level summary saved in results/CONFIG/. Contains:  
  - config: The full configuration tag used for this run’s results folder.  
  - all_sensors_length: Initial number of sensors at the start of the elimination loop, N.  
  - sensor_idx: Baseline index reference for the original sensor order, i.e., [1, 2, …, N]. Helpful in mapping any later indices back to the initial ordering.  
  - avg_dist_vec: Average inter-class centroid separation at each elimination step (computed in the PCA space used for plotting).  
  - Ordering: entry j corresponds to the configuration with (N − j + 1) remaining sensors. Thus:  
    - avg_dist_vec(1) → all N sensors,  
    - avg_dist_vec(2) → N−1 sensors,  
    - …,  
    - avg_dist_vec(N) → 1 sensor.  
  - Vektor_ARI: The adjusted Rand index of the classifier versus the ground truth labels vs. each iteration, in reverse ordering, hence Vektor_ARI(j) corresponds to (N − j + 1) remaining sensors, with the last element for the 1-sensor case.  
  - Gaussian_error_val_mat: Per-class Gaussian misclassification estimates across the elimination path. The column j hold each class error values for the configuration with (N − j + 1) remaining sensors (last column = 1 sensor).  
    Notes: values reflect the decision rule in use (_QDA or _Vor) and the inflation/buffer settings; near-zero values indicate very low estimated error for that class at that step.  
  - supertitle: The dataset/run title you provided (also used as the .mat filename).  
- eliminated_one.mat: cell array logging removed sensors in order.

## References:
[1] Faran, Michael, et al. “The ACFSA V2.0: Analyte Classification and Feature Selection Algorithm.” Manuscript submitted (2025).  
[2] Petresky, Gabriel, et al. "Metal-Ion Optical Fingerprinting Sensor Selection via an Analyte Classification and Feature Selection Algorithm." Analytical Chemistry 97.16 (2025): 8821-8832.  
[3] Lee, Yeon Soo, et al. "Spatiotemporal molecular tracing of ultralow-volume biofluids via a soft skin-adaptive optical monolithic patch sensor." Nature Communications 16.1 (2025): 3272.‏  
[4] Yoon, Minyeong, et al. "Enzyme-free optical detection of uric acid using corona phase molecular recognition in near-infrared fluorescent single-walled carbon nanotubes." Nanoscale 17.17 (2025): 10652-10662.‏  

The code uses adjusted versions of the following MATLAB codes:  
plot_elipse.m, taken from:  
[5] https://www.mathworks.com/matlabcentral/fileexchange/116610-plot-ellipse-on-scattered-2d-data?s_tid=prof_contriblnk  
[6] Ohad Gal (2025). fit_ellipse (https://www.mathworks.com/matlabcentral/fileexchange/3215-fit_ellipse), MATLAB Central File Exchange. Retrieved November 15, 2025.  
[7] Jakob Sievers (2025). VoronoiLimit(varargin) (https://www.mathworks.com/matlabcentral/fileexchange/34428-voronoilimit-varargin), MATLAB Central File Exchange. Retrieved November 15, 2025.  
[8] qqffssxx (2025). Rand and Adjusted Rand Index Calculator for Cluster Analysis (https://www.mathworks.com/matlabcentral/fileexchange/130779-rand-and-adjusted-rand-index-calculator-for-cluster-analysis), MATLAB Central File Exchange. Retrieved November 15, 2025.  

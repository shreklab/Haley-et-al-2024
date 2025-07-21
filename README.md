# haley_et_al_2024
> Analysis package for [Haley et al. (2024)](https://doi.org/10.7554/eLife.103191.2)

This code was written and used to analyze the behavior of *C. elegans* exploring environments with small, dilute patches of bacteria. This package is comprehensive and includes all code used for this paper. Some files are only useful for this particular set of experiments and require adaption for other uses.

# Table of Contents
1. [Software](#software)
2. [File list](#file-list)
   - [Compile experiment metadata](#foragingInfo)
     - [Get file metadata](#metadata)
     - [Get arena properties](#arena)
     - [Get lawn properties](#lawn)
     - [Image registration](#registration)
     - [Create grids](#grid)
   - [Analyze worm behavior](#behavior)
     - [Analyze worm movement and identify encounters](#analysis)
     - [Plot behavior](#plot)
   - [Classify worm behavior](#label)
   - [Estimate bacterial density](#OP50-GFP)
   - [Model the exploitation decision](#model)
   - [Statistical analysis](#statistics)
   - [Figures](#figures)
   - [Figure supplements](#supplements)
3. [Links](#links)
4. [Contact](#contact)

---
## Software <a name="software"></a>

I used the following software for analysis and figure production:
- [MATLAB 2024a](https://www.mathworks.com/products/matlab.html)
  - Image Processing Toolbox
  - Statistics and Machine Learning Toolbox
  - [Bioformats Image Toolbox](https://github.com/Biofrontiers-ALMC/bioformats-matlab)
  - [Violinplot-Matlab](https://github.com/bastibe/Violinplot-Matlab)
- [WormLab](https://www.mbfbioscience.com/products/wormlab/)
- Adobe Illustrator 2024
- Adobe Photoshop 2024

---
## File list <a name="file-list"></a>

### Compile experiment metadata <a name="foragingInfo"></a>
- `getForagingInfo.m` : compiles all of the metadata for a set of experiments outlined in an **infoFile** (`.xls`) and saves that info as a table in a `.mat` file. This metadata includes information about the experiment (e.g. time stamps, conditions, camera, temperature), the video (e.g. pixels, frame rate), the behavioral arena (e.g. diameter, mask, orientation), and the lawns (e.g. mask, size, spacing).

  - #### Get file metadata <a name="metadata"></a>

    - `getExperimentInfo.m` : retrieves metadata from a `.xls` file.
    - `getVideoInfo.m` : retrieves metadata from a `.avi` video file.

  - #### Get arena properties <a name="arena"></a>

    - `getAcetateArena.m` : takes an image from a video containing an acetate arena, locates the arena, and calculates information about the video.
    - `getAcetateOrientation.m` : takes an image from a video containing both an acetate arena and a reference dot, locates the reference dot, and calculates information about the orientation of the arena.

  - #### Get lawn properties <a name="lawn"></a>

    - `getInvisibleLawns.m` : takes the experimental video and the contrast video (if applicable) and finds/estimates the location of bacterial lawns in the experiment.
      - `getLawnsTemplate.m` : creates an estimate of the locations of bacterial lawns in the behavioral arena given information about size and spacing of a single lawn or an isometric grid of lawns.
      - `getLawnsThreshold.m` : creates an estimate of the locations of bacterial lawns in the behavioral arena given a frame from the video (ideally without worms) and information about the expected position of the lawns using image processing and morphological operations. `getLawnsThreshold.m` expects that the image used is not necessarily from the experimental video itself and enables translation, rotation, and scaling of the lawn mask to match the experimental video.
      - `getLawnsFilter.m` : creates an estimate of the locations of bacterial lawns in the behavioral arena given a VideoReader object (ideally for a video without worms) and information about the expected position of the lawns using image processing and morphological operations. `getLawnsFilter.m` expects a video where a dark object (e.g. black cardstock) is passed between the lightsource and the experimental plate, creating a "darkfield" image of plate. The cardstock scatters the light creating increased contrast of the bacteral lawns, enabling automatic identification of the lawn edges. `getLawnsFilter.m` expects that the video used contains an arena that is not necessarily in the same position and orientation as the experimental video itself and enables translation, rotation, and scaling of the lawn mask to match the experimental video.
      - `getLawnsCircle.m` : creates an estimate of the locations of bacterial lawns in the behavioral arena given a VideoReader object (ideally for a video without worms) and information about the expected position of the lawns using image processing and morphological operations. `getLawnsCircle.m` expects a video where a dark object (e.g. black cardstock) is passed between the lightsource and the experimental plate, creating a "darkfield" image of plate. The cardstock scatters the light creating increased contrast of the bacteral lawns, enabling automatic identification of the lawn edges. `getLawnsCircle.m` uses [`imfindcircles`](https://www.mathworks.com/help/images/ref/imfindcircles.html) to identify the putative lawns on every frame and creates a weighted average image. `getLawnsCircle.m` expects that the video used contains an arena that is not necessarily in the same position and orientation as the experimental video itself and enables translation, rotation, and scaling of the lawn mask to match the experimental video.
      - `getLawnProperties.m` :  takes a binary mask of regions of interest (e.g., bacterial lawns), labels each distinct region, and computes their geometric properties.
    - `getManualLawns.m` : updates lawn locations with any manual corrections that were made in Adobe Photoshop (when needed).

  - #### Image registration <a name="registration"></a>

    - `registerImage.m` : takes in two images, detects SURF features (using [`detectSURFFeatures`](https://www.mathworks.com/help/vision/ref/detectsurffeatures.html)) of each, estimates a 2D rigid transformation (using [`estgeotform2d`](https://www.mathworks.com/help/vision/ref/estgeotform2d.html)), and performs the transformation (using [`imwarp`](https://www.mathworks.com/help/images/ref/imwarp.html)). 
    - `transformLabeledImage.m` : takes in an image and its properties as a structure and performs given image transformations (i.e. registrations); performed on the grids (hexagon, diamond, triangle) to register them to arena. This function use simple geometric tranformations (e.g., [`imrotate`](https://www.mathworks.com/help/images/ref/imrotate.html)) and is intended for use with images that have been labeled using [`bwlabel`](https://www.mathworks.com/help/images/ref/bwlabel.html). If using a video frame, use `registerImage.m` for more accurate image registration.

  - #### Create grids <a name="grid"></a>

    - `createHexagonalGrid.m` : creates a grid of isometric hexagons and corresponding triangles (6 per hexagon) given information about their size and spacing.
    - `createDiamondGrid.m` : creates a grid of diamonds and corresponding triangles (2 per diamond) given information about their size and spacing.

### Analyze worm behavior <a name="behavior"></a>
- `analyzeForaging.m` : uses the body segment location of animal combined with the location of bacterial patches/lawns to identify encounters and properties of those encounters.

   -  #### Analyze worm movement and identify encounters <a name="analysis"></a>
    
      - `analyzeWormLabTracks.m` : takes body segment data exported from WormLab and computes numerous metrics related to the worm's position, velocity, and turning behavior as well as it's location relative to the arena and lawns.
      - `offsetTime.m` : corrects the wormNum identifiers and time values for experiments where the recording spans multiple `.avi` files.
      - `defineEncounter.m` : uses high-resolution body segment locations of animals exploring environments with a single small patch of bacteria to define an encounter event for use when only the mid-point of the animal can be reliably tracked.
      - `analyzeEncounters.m` : uses the location of animal(s) and lawn(s) to identify encounters and properties of those encounters.
   
    - #### Plot behavior <a name="plot"></a>
      - `plotTracks.m` : writes an image to file showing the tracks of an animal overlaid onto an image of the lawns and arena. Tracks are colored based on a given metric (e.g., time, velocity, path angle).
      - `createVideo.m` : writes a video to file showing either the tracks of an animal or the original video downsampled with scale bar and time stamp showing.

### Classify worm behavior <a name="label"></a>
- `labelEncounters.m` : uses the body segment location of an animal combined with the location of bacterial patches/lawns to identify encounters and properties of those encounters.
   - `estimatePosterior.m` : predicts the conditional probabilities of sensing for encounters where only minimum velocity on-patch was a reliable metric. Marginalizes the conditional probabilities over the other the other two metrics and numerically integrates using an adaptive quadrature method over the product of the QDA estimated conditional probability distribution and the kernel density estimate of the joint probability distribution.

### Estimate bacterial density <a name="OP50-GFP"></a>
- `analyzeGFP.m` : gets all the meta data for a set of experiments, identifies brightfield images of acetate templates, identifies background images of empty plates, and extracts and analyzes the fluorescence intensity profile for each patch.
   - #### Get file metadata <a name="metadata2"></a>
      - `getPlateInfo.m` : retrieves metadata from a `.xls` file.
      - `getMetaDataCZI.m` : retrieves metadata from `.czi` microscopy image files.
   - #### Analyze fluorescence <a name="fluorescence"></a>
      - `getFluorescenceBackground.m` : estimates the image of the microscopy light source's diffusion pattern on an empty agar plate by fitting an elliptic paraboloid OR by smoothing the image using an averaging filter.
      - `analyzeLawnProfiles.m` : normalizes a fluorescent image of bacterial patches, detects each patch, and then extracts and analyzes their radial intensity profiles to quantify properties like peak fluorescence.
- `plotForagingGFP.m` : loads a saved workspace from the `analyzeGFP.m` analysis pipeline, performs statistical analysis, and generates plots and source data related to the fluorescence profiles of bacterial lawns.
- `borderAmplitude.csv`: contains the parameters for linear models describing the change in bacterial density over time for various experimental conditions. The models were generated by `plotForagingGFP.m`, which fits a linear regression to the "border amplitude"—a measure of fluorescence intensity at the patch edge—as a function of growth time. The `slope` and `intercept` values from this file are used in other scripts to estimate the effective bacterial density of a patch at any given time during an experiment.

### Model the exploitation decision <a name="model"></a>
- `modelExploit.m` : fits and evaluates generalized linear models (GLMs) to predict a worm's decision to exploit a food patch.

### Statistical analyses <a name="statistics"></a>
- `benjaminiHochberg.m` : corrects for multiple comparisons by controlling the false discovery rate (FDR) using the Benjamini-Hochberg procedure.
- `silvermansTest.m` : computes a bandwidth using Silverman's rule of thumb and iteratively finds the "critical bandwidth" where the number of modes (peaks) in the KDE switches from multimodal to unimodal.
- `gaussianContours.m` : calculates the coordinates for plotting ellipses that represent the standard deviation contours of a 2D Gaussian distribution.
- `writeSourceData.m` : compiles and writes source data, including a header, multiple data tables with titles, and accompanying notes, into a single, formatted `.xlsx` file.
- `permutePatchLocation.m` : generates a null model for foraging behavior by randomizing the spatial locations of food patches. For each experimental plate, the script creates thousands of new layouts by individually rotating and translating the original patches to new, non-overlapping locations within the arena. It then replays the worms' actual recorded trajectories on these shuffled layouts to calculate a null distribution of patch occupancy time, which can be used to assess whether the observed behavior is significantly different from random chance.

### Figures <a name="figures"></a>
- `fig1.m` : generates figures showing fundamental foraging metrics, including patch occupancy, encounter duration, and velocity changes upon patch entry.
- `fig2.m` : generates figures related to behavioral classification, including velocity profiles, time on patch, and encounter type probabilities over time.
- `fig3_mini.m` and `fig3_single.m` : generate figures for experiments with arrays of small patches and single large patches, respectively, focusing on relative bacterial density and the time/number of encounters before exploitation.
- `fig4.m`, `fig4_matching.m`, and `fig4_sensory.m` : generate figures for the GLM analysis, including model coefficients and predicted vs. observed probabilities of exploitation for different experimental paradigms.

### Figure supplements <a name="supplements"></a>
- `figS1.m` : assembles graphics illustrating the experimental and video processing workflow.
- `figS2.m` : generates figures explaining the encounter definition methodology, using head and midpoint tracking data.
- `figS3.m` : creates plots for the GMM used to classify encounters, including the effect of censored data.
- `figS4.m` : illustrates the patch location permutation process for generating the null model.
- `figS5.m` : shows the quantification of deceleration upon patch entry.
- `figS6.m` : details the GFP analysis pipeline for quantifying bacterial density, from raw images to extracted fluorescence profiles.
- `figS7_S8.m` : displays example images and border amplitude measurements for various patch conditions in the GFP experiments.
- `figS9.m` : generates example worm traces for the single large patch experiments.
- `figS10.m` : shows on-patch probability over time for all conditions compared to a permuted null model.
- `figS11_S12.m` : srovides detailed analysis of the GMM and QDA classifiers, including Silverman's test for bimodality and 3D visualization of the QDA decision boundary.
- `figS13.m` : generates example worm traces for the large single patch concentration experiments.
- `figS14.m` : generates example worm traces for the small patch array ("mini") experiments.
- `figS15.m` : compares the performance (Log-Likelihood, AIC, BIC) of different GLM combinations and visualizes model coefficients.
- `figS16.m` : plots the probability of the first exploitation event versus encounter number for all conditions and models.
- `figS17.m` : generates figures for the food-deprivation experiments, including encounter timelines and GLM coefficients.
- `figS18.m` : generates figures for the matching law experiments, including example traces and comparisons of GLM models with different time-dependent terms.
- `figS19.m` : generates figures for the sensory mutant experiments, including ridge regression optimization and GLM coefficients.
- `VideoSupp.m` : generates supplementary videos, including raw behavior, contrast videos for lawn detection, and tracking overlays.
- `reportMetaData.m` : summarizes metadata for all experiments in the project.

---
## Links <a name="links"></a>

- Paper: [doi.org/10.7554/eLife.103191.2](https://doi.org/10.7554/eLife.103191.2)
- Git Repository: [github.com/jesshaley/haley_et_al_2024](https://github.com/jesshaley/haley_et_al_2024)
- Related projects: [github.com/shreklab](https://github.com/shreklab)

---
## Contact <a name="contact"></a>

If you are having issues or have questions about the code, please contact jess.allison.haley at gmail.com.
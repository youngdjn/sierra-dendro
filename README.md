# Sierra dendro project

Analysis of ~800 tree cores collected across latitudinal and elevational gradients in the Sierra Nevada

## Data files for analysis

The analysis-ready data files are located in: [data/compiled-for-analysis/](https://github.com/youngdjn/sierra-dendro/tree/master/data/compiled-for-analysis)

* **plots.csv**: Plot-level data, including coordinates and normal climate. One row per plot.\
* **trees.csv**: Tree-level data, including species, DBH, plot, and voronoi area. One row per tree.\
* **years.csv**: Year-level data, including ring widths and annual climate. One row per year per tree.\


#### Sample size summary tables

*Also see "Quality filtering", below, for more details on these tables*

* **[Sample size by plot and species](https://raw.githack.com/youngdjn/sierra-dendro/master/data/dendro/sample_size_summaries/plot_sample_size.html)**\ (The rightmost columns show the number of trees with >15, >30, and >50 rings for each plot-species combination)

* **[Number of rings and truncation years for each core](https://raw.githack.com/youngdjn/sierra-dendro/master/data/dendro/sample_size_summaries/core_ring_counts.html)** 


#### Details on the ring-width data columns (in years.csv)

* **rwi**: Ring-width index, computed by dividing the raw width by a 30-year spline\
* **raw_width**: The raw width in mm\
* **ba**: The basal area of the tree including that year's growth. Only computed for trees measured to the pith.\
* **bai**: The basal area increment (the additional basal area grown) that year. Only computed for trees measured to the pith.\
* **bai.ba**: Basal increment divided by previous year's basal area. Equivalent to relative growth rate of basal area. Only computed for trees measured to the pith.\
* **ba.prev**: Previous year's BA. Only computed for trees measured to the pith.\
* **ba.ext, bai.ext, bai.ba.ext, ba.extprev**: Equivalent to the above metrics, but computed based on the tree's measured DBH, using a regression to relate external DBH (which includes bark) to internal DBH (which doesn't). Only computed for trees where the core extends to the bark.\
* **ba.comb, bai.comb, bai.ba.comb, ba.prev.comb**: Equivalent to the above metrics, but they are a composite the uses internal BA when it exists and external BA when it doesn't. Still does not cover all trees, because some cores did not extend to the pith AND the core did not extend to the bark.

#### Details on the climate/weather data columns (in years.csv)

* Climate variables without a number suffix are current-year
* Climate variables with a number suffix X reflect the weather X years ago.
* Climate variables with a z reflect the z-score (number of standard deviations from the long-term mean)
* More climate variables can be extracted relatively easily. I still have to port the scripts to the git repository though.



## Quality filtering of tree ring data

Tree ring data are filtered (i.e., ring widths are excluded from the data files) in the following ways:

* If the crossdating technician reported multiple equally-good ways to get good alignment, or no good ways, all ring widths are excluded\
* If the crossdating technician reported adding or removing a ring that was not justified based on the core image (but substantially improved alignment), ring widths from that year and older are excluded (note that this filter is easy to turn off).\
* Using a moving window, the core is compared against its cluster-level reference chronology. If the correlation drops below a certain threshold, the ring widths older than that point are excluded, **unless** the correlation improves again within X years. Improved is defined as being above the correlation threshold for at least Y years. At the time I wrote this, X was 15, Y was 10, the window width was 12, and the correlation threshold was 0.1. All these parameters can easily be tuned. The reference chronology is computed excluding cores that had a ring addition or deletion that was not justified based on the image. The correlation between the focal and reference is evaluated after normalizing both by computing each ring's proportional difference from the average of the two rings adjacent to it.

Filtering (along with many other tasks) is done by the script compile_data_for_analysis.R (in scripts/analysis-prep/).

Here is a table to visualize when and why the cores in each plot are being truncated:\
**[Number of rings and truncation years for each core](https://raw.githack.com/youngdjn/sierra-dendro/master/data/dendro/sample_size_summaries/core_ring_counts.html)**

The truncation reasons are:

* **Unjustified mod**: there was a ring added or removed that did not make sense based on the image
* **Ab**: there was an aberration in the core that obscured the rings
* **Poor align throughout**: The technician could not find a good way to get alignment
* **Multiple ways**: There were multiple equally good ways to get good alignment
* **Alignment drops**: The chronology was truncated by the moving-window quality filter
If there is more than one code, the first one happened more recently.

Here is a table to visualize the sample sizes per plot. The rightmost columns show the number of trees with >15, >30, and >50 rings for each plot-species combination:\
**[Sample size by plot and species](https://raw.githack.com/youngdjn/sierra-dendro/master/data/dendro/sample_size_summaries/plot_sample_size.html)**

Note that we still have to decide whether it is important to do a more local-level crossdating.


## Compiling data for analysis

This has already been done (and the output of the most recent run is in data/compiled-for-analysis/). This only has to be done if there are changes to the underlying raw data (e.g., plot data or ring width measurements) or if the tree ring quality filtering parameters change.

Note that compiling analysis-ready data requires some large files that are too big to sync to GitHub. These go in the "data/non-synced" folder of the repo (which you have to create after you clone the repo). They are stored [on Box here](https://ucdavis.box.com/s/3j6dnkzjuyhi3vrarfye81apbbf2ezzi).

The above text documents the data that is already prepared for analysis. Here are details on the scripts that process the raw data (also in the repo) into the analysis-ready data. Scripts are in the folder  Data prep scripts should be run in the following order:

**Process field-based plot data** (scripts in scripts/plot-tree/)

1. **tree_loc_calc.R**: Compute tree locations based on distance and bearing from reference-point coordinates
2. **extract_tree_elevation.R**: Extract elevation from DEM for each tree
3. **compute_tree_cluter.R**: Determine the cluter (NL, SL, NH, SH) that each tree belongs to.
4. **voronoi_calc.R**: Extract the voronoi polygon area for each tree.

**Process climate/weather data** (scripts in scripts/climate-prep)

5. *Section not yet completed. Scripts and climate source data files must be ported to the git repo. The scripts are already in the repo, but paths etc need to be updated, and the climate data they depend on must be added. The main script (which runs all the sub-scripts is dendro_climate_compile.R). This step produces the large year-by-tree climate data file tree_clim_full.csv in "data/non-synced/climate-extracted/"*

**Compile all ring width, plot, and climate/weather data for analysis** (script in scripts/analysis-prep/)

6. **compile_data_for_analysis.R**: Reads in ring-width data, performs quality filtering, and reads in and merges with plot and climate data. This script produces the final files for analysis in data/compiled-for-analysis/. This is the script that contains the parameters to tweak for the moving-window quality filter.


## Raw data files

*This section is under construction. It currently only lists selected major data input files/folders.*

Subfolders of the "data/" folder:

* dendro/coorecorder-measurements: CooRecorder .pos files
* dendro/crossdating-records: Technician records of crossdating (e.g., when there was no clear way to get good alignment)
* dendro/reference-chronologies: Auto-generated reference chronologies.
* dendro/sample-size-summaries: Summaries of the sample sizes of the quality-filtered tree ring data
* plot-and-tree/field-data: Plot-data from the field (e.g., tree DBH, competitors)
* plot-and-tree/processed: Processed plot data (e.g., tree locations, voronoi area, processed tree climate data)
  * trees_loc.csv, in addition to the tree locations, this file contains most other measured plot data and is the main file to use to get plot-based data for each tree. This is not technically raw data because it is compiled based on data in the raw plot data using the scripts described above.
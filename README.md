
<!-- README.md is generated from README.Rmd. Please edit that file -->

# shewhart4hrms

<!-- badges: start -->
<!-- badges: end -->

A package for automated processing of LC-HRMS data to record intensity,
retention time, mass accuracy and peak width of user defined internal
standards. The data is displayed in the form of Shewhart control charts
and density distributions along with user defined warning limits for
quick visual indication of exceedances.

## Installation of shewhart4hrms

R version 4.4.2 and `ntsworkflow` are required (see
github.com/bafg-bund/ntsworkflow).

First follow the installation instructions for `ntsworkflow`, then run

``` r
renv::install("bafg-bund/shewhart4hrms")
```

## Usage

A local directory tree to hold the raw data, settings and processing
results is created first. The directory can be given any name.

``` r
shewhart4hrms::newDirTree("~/linsensuppe")
```

In this directory there are sub-directories `mzXML-data-files`,
`results` and `settings`.

The subdirectory `settings` contains two `is-table-*.csv` files which
lists the internal standards you are using (separated for positive and
negative ionization). Example tables are provided. You must provide a
name, formula, retention time and adduct for each internal standard. See
the example files for the formating. For example, deuterium atoms are
indicated with \[space\]2H.

| name               | formula          |   rt | adduct   |
|:-------------------|:-----------------|-----:|:---------|
| Bezafibrat-d4      | C19H16 2H4ClNO4  | 10.9 | \[M+H\]+ |
| Olmesartan_acid-d6 | C24H20 2H6N6O3   |  7.3 | \[M+H\]+ |
| Iopromide-d3       | C18H21 2H3I3N3O8 |  4.7 | \[M+H\]+ |

Example is-table-pos.csv file

The `settings/processingSettings.yml` file contains all settings for
measurement file processing. This is a Yaml formatted text file best
viewed with Rstudio or notepad++. For more information on Yaml see
[yaml.org](https://yaml.org/). Some default values used on a Sciex 6600
are already filled in but these will most likely be very different
depending on the mass spectrometer. Especially the `area_threshold_*`
parameter may be as high as 1e6 depending on the mass spectrometer used.

Example processingSettings.yml file

``` yaml
# Settings file for shewhart4hrms ##############################################
# Processing settings ##########################################################

# mz_tol and eic_extraction width in Da
# rt_tol in min
# area_threshold in counts

# pos 
file_pattern_pos: .*pos.*
mz_tol_pos: 0.005
rt_tol_pos: 0.5
area_threshold_pos: 1
eic_extraction_width_pos: 0.04

# neg
file_pattern_neg: .*neg.*
mz_tol_neg: 0.005
rt_tol_neg: 0.5
area_threshold_neg: 1
eic_extraction_width_neg: 0.04
```

In the `settings/warningLevels.csv` file the user can specify where
warning levels are to be drawn on the plots. The name of the internal
standard is given under name_is and must match the name given in the
`is-table-*.csv`. The parameter name must match the parameter name
displayed in the app. Missing values are left empty.

| name_is           | parameter      | polarity | upper_warning | lower_warning |
|:------------------|:---------------|:---------|--------------:|--------------:|
| DummyIS1          | intensity      | pos      |         600.0 |           400 |
| DummyIS3HighInten | intensity      | pos      |               |         40000 |
| DummyIS1          | area           | pos      |        1500.0 |          1800 |
| DummyIS1          | delta_mz_mDa   | pos      |           1.0 |               |
| DummyIS1          | delta_rt_min   | pos      |           0.3 |               |
| DummyIS1          | peak_width_min | pos      |           0.2 |               |

Example warningLevels.csv

You must first copy the mzXML files of your measurements into the
`mzXML-data-files` sub-directory. For conversion to mzXML use
[Proteowizard](https://proteowizard.sourceforge.io). Each file must be
an LC-ESI-HRMS measurement of the defined internal standards (see
“settings”). The file name must contain either “pos” or “neg” to
indicate the polarity.

Once you have your first mzXML files copied into `mzXML-data-files`
(both pos and neg files are needed) and an initial guess at the
appropriate integration parameters has been set, the processing can be
initiated.

Start the app by running:

``` r
shewhart4hrms::viewShewhart("~/linsensuppe")
```

The shiny interface should open in your browser.

Click “Refresh” to start the processing. The results of the peak
integration, if successful, will be stored in the `results` directory
and the plots should be updated. To restart the processing delete the
`results_*.csv` file.

New files can be copied into the `mzXML-data-files` at any time and
clicking “Refresh” in the shiny window will process the new files.

To start the shiny interface directly in Windows you can use a `.bat`
file. An example file is provided in the directory you created, and
should work as-is regardless of where it is placed on the system.

## License

Copyright 2025 Bundesanstalt für Gewässserkunde (German Federal
Institute of Hydrology) Package and source code available under GPL 3
(or greater).

shewhart4hrms is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

shewhart4hrms is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
Public License for more details.

You should have received a copy of the GNU General Public License along
with shewhart4hrms. If not, see [gnu.org](https://www.gnu.org/licenses).

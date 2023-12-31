# shewhart4hrms

A package for automated processing of LC-HRMS data to record intensity, retention time, mass accuracy and peak width of user defined internal standards. The data is displayed in the form of Shewhart control charts and density distributions along with user defined warning limits for quick visual indication of exceedances.


## License
Copyright 2020-2022 Bundesanstalt für Gewässserkunde (German Federal Institute of Hydrology)
Package and source code available under GPL 3 (or greater). 

shewhart4hrms is free software: you can redistribute it and/or modify it under the 
terms of the GNU General Public License as published by the Free Software 
Foundation, either version 3 of the License, or (at your option) any 
later version.
 
shewhart4hrms is distributed in the hope that it will be useful, but WITHOUT ANY 
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 
You should have received a copy of the GNU General Public License along 
with shewhart4hrms. If not, see [gnu.org](https://www.gnu.org/licenses).


## Installation of shewhart4hrms

R version 4.2.1 with the following cran packages: `Rcpp`, `RcppArmadillo`, `shiny`, `tidyverse` and the Bioconductor package `xcms` is required.

Start R and run the following (an Internet connection is required):

```
install.packages(c(
  "Rcpp",
  "RcppArmadillo",
  "shiny",
  "tidyverse",
  "cowplot",
  "BiocManager"
))
BiocManager::install("xcms")
```

On Windows, the compiled package archive `shewhart4hrms_0.1.zip` can then be installed through the R-GUI menu "Packages" => "Install packages from local files..." or through the RStudio menu "Tools" => "Install Packages..." => "Install from" => "Package archive file". The compiled archive can be found on the BfG FTP server: `ftp.bafg.de/pub/REFERATE/g2/quality-control-shewhart4hrms/`

Alternatively, you can compile the package for your system using the `devtools` package and the source code. Only Windows 10 has been tested so far.

## Usage

A local directory tree to hold the raw data, settings and processing results is created first. The directory can be given any name. 

```
shewhart4hrms::create_dir("~/linsensuppe")
```

In this directory there are sub-directories `mzXML-data-files`, `results` and `settings`.

The subdirectory `settings` contains two `is-table-*.csv` files with the tables (positive and negative ionisation) of internal standards you are using. Example tables are provided. You must provide a name, formula, retention time and adduct for each internal standard. See the example files for the formating. For example, deuterium atoms are indicated with [space]2H. The `settings.yml` file contains all settings for peak integration, processing and warning limits. This is a Yaml formatted text file best viewed with Rstudio or notepad++. For more information on Yaml see [yaml.org](https://yaml.org/). Some default values used on a Sciex 6600 are already filled in but these will most likely be very different depending on the mass spectrometer. Especially the `area_threshold_*` parameter may be as high as 1e6 depending on the mass spectrometer used.

You must first copy the mzXML files of your measurements into the `mzXML-data-files` sub-directory. For conversion to mzXML use [Proteowizard](https://proteowizard.sourceforge.io). Each file must be an LC-ESI-HRMS measurement of the defined internal standards (see "settings"). The file name must contain either "pos" or "neg" to indicate the polarity.

Once you have your first mzXML files copied into `mzXML-data-files` (both pos and neg files are needed) and an initial guess at the appropriate integration parameters has been set, the automated processing can be initiated.

```
shewhart4hrms::initiate("~/linsensuppe")
```

The results of the peak integration, if successful, will be stored in the `results` sub-directory. To view the results run:
```
shewhart4hrms::viewshewhart("~/linsensuppe")
```
The shiny interface should open in your browser showing the results from one data-file in each polarity (one point on the time-series for each category).

To process the remaining files click the "Refresh" button. If processing is successful, the plots should be updated. 

New files can be copied into the `mzXML-data-files` at any time and clicking "Refresh" in the shiny window will process the new files.

To start the shiny interface directly in Windows you can use a `.bat` file. An example bat file is provided in the directory you created, and should work as-is regardless of where it is placed on the system.



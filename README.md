# disnap

`disnap` is an R package for analysis of disease snapshot data.


## Installation

To install the github version of `disnap`:

```r
devtools::install_github("mc30/disnap")
```


## Principal purpose

This package is used to analyse disease snapshot data, i.e. data on infection status of individual units (locations/entities) taken at discrete points in time, e.g. annual data on disease occurrence from farms located in a particular region.

Main functionality is to provide functions to compute several statistical properties of a data set and associated data (e.g. properties of units, links between units)

This package treats the temporal component as a discrete value and is not suitable to perform continuous analysis of the temporal dimension.


## Input data

Disease snapshot data are represented by a matrix with [number of units] rows and [Number of timepoints] columns, where each element is a disease status: 0 (negative), positive number (infected, can be used as a strain), NA (unit is not present or status is not defined).


## Further development

The following features will be added in next versions:
* Analysis of network data
* Performance improvement for large data sets


## Feedback

If you have any comments, suggestions or a portion of criticism, please contact me at <mikhail.churakov@gmail.com>.

Pull requests and general comments are welcome.

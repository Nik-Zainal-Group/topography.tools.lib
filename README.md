# topography.tools.lib R package


## Table of content

- [Introduction to the package](#intro)
- [Versions](#version)
- [Installation](#installation)
- [Testing the package](#test)
- [Package documentation](#docs)
- [Functions provided by the package](#functions)

<a name="intro"></a>

## Introduction to the package

```topography.tools.lib``` is an R package for the analysis of genomic features.
This package is designed to provide low level functions that can be used to
describe, manipulate and analyse generic genomic features. There are two types
of genomic features considered: positions and bed regions. Each position or region
is defined by a genomic location: chromosome and position for positions, and 
chromosome, start and end for bed regions. They also require a unique id, and, 
optionally, a class attribute, which allows to aggregate results by class.

Functions available in the package include calculating the intersection between 
two sets of genomic features and their statistical significance, as well as testing
a set of bed regions for overlaps and aggregating overlapping regions. Functions
for calculating intersections can also be used for annotating a set of positions
or bed regions with other genomic features.

<a name="version"></a>

## Versions

1.0.0

- first version contains basic functions for positions and bed regions intersections and their statistical significance

<a name="installation"></a>

## Installation

Download the `topography.tools.lib` repository:

```
git clone https://github.com/Nik-Zainal-Group/topography.tools.lib.git
cd topography.tools.lib
```

You can install ```topography.tools.lib``` by entering the R environment from the main
directory and typing:

```
install.packages("devtools")
devtools::install()
```

<a name="test"></a>

## Testing the package

You can test the package by entering the package main directory and
typing from the R environment:

```
devtools::test()
```

<a name="docs"/>

## Package documentation

**DOCUMENTATION:** Documentation for each of the functions below is
provided as R documentation, and it is installed along with the R
package. The documentation should give detailed explanation of the input
data required, such as a list of data frame columns and their
explanation. To access the documentation you can use the ```?function```
syntax in R, for each of the functions below. For example, in R or
RStudio, type ```?sortChroms```.

<a name="functions"></a>

## Functions provided by the package

Function for basic manipulation of positions and bed regions:

- **```sortPositions```**: order positions according to chromosome name
and base pair position
- **```sortBedRegions```**: order bed regions according to chromosome name
and base pair position
- **```checkBedRegionsOverlap```**: given a set of bed regions, determine
whether there is any overlap and if so return in which chromosomes the overlaps
were observed
- **```breakDownOverlappingBedRegions```**: given a set of bed regions that
contain overlaps, break the regions so that the resulting set of bed regions
do not overlap any more. If a ```signal``` column is specified, then signal of
overlapping segments can be aggregated either as a sum or average
- **```assignBedRegionsToNonOverlappingSets```**: given a set of bed regions,
group the regions into sets such that each set contains non-overlapping regions
- **```getIMD```**: calculate positions inter-mutational distance
- **```getIRD```**: calculate bed regions inter-region distance, for non-overlapping
bed regions

Function for intersecting positions and bed regions:

- **```intersectPositionsAndBedRegions```**: calculate the intersection of a set of
positions and a set of bed regions
- **```intersectBed```**: calculate the intersection between two 
sets of bed regions

Function for calculating the statistical significance of positions and bed regions
intersections:

- **```correlatePositionsAndBedRegions```**: calculate the intersection of a set of
positions and a set of bed regions and determine the statistical significance
- **```correlateBed```**: calculate the intersection between two 
sets of bed regions and determine the statistical significance
- **```multipleCorrelations```**: test the overlap between one genomic feature and a
list of genomic features

Function for random resampling positions or bed regions:

- **```resamplePositions```**: random resampling of a set of positions preserving the
chromosome
- **```resampleBedRegions```**: random resampling of a set of regions preserving the
chromosome

Function for annotating positions or bed regions:

- **```annotateBedWithGenes```**: annotate a set of bed regions with genes



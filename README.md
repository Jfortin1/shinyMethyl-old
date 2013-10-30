shinyMethyl
===========

Interactive visualization of the 450k methylation data

Some changes and bugs fixed: 

- Within the predictGender panel, it is now possible to download a .csv file containing not only the predicted gender, but also the actual gender if it was provided in the phenotype data, it will be also included in the file. The function also handles missing data. Thanks to Brent Pedersen
- With the update of minfi in Bioconductor 2.13, IlluminaHumanMethylation450kannotation.ilmn.v1.2 is no longer supported. ShinyMethyl is now updated to be annotation free (shinyMethyl does not really need an annotation). Thanks to Kathleen Fisch for pointing that out.

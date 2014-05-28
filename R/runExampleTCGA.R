runExampleTCGA <- function(){
	summary.tcga.raw <- get(data("summary.tcga.raw", package = "shinyMethylData"))
	summary.tcga.norm <- get(data("summary.tcga.norm", package = "shinyMethylData"))
	runShinyMethyl(summary.tcga.raw, summary.tcga.norm)
}
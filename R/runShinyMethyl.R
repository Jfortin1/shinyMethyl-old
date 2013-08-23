runShinyMethyl <-
function(extractedData){
 
     directory <- system.file(package="shinyMethyl","shinyMethyl")
     require(shiny)
     shinyMethylData <<- extractedData
 
 
     runApp(directory)
 
}

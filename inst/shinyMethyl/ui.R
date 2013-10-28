#### Created by Jean-Philippe Fortin
#### Aug 19 2013


library(RColorBrewer)


############################################################


betaQuantiles <- shinyMethylData$betaQuantiles
mQuantiles <- shinyMethylData$mQuantiles
methQuantiles <- shinyMethylData$methQuantiles
unmethQuantiles <- shinyMethylData$unmethQuantiles
greenControls <- shinyMethylData$greenControls
redControls <- shinyMethylData$redControls
# oobControls <- shinyMethylData$oobControls
XYMedians <- shinyMethylData$XYMedians
sampleDistance <- shinyMethylData$sampleDistance
covariates <- shinyMethylData$pd
covariates <- covariates[match(colnames(shinyMethylData$betaQuantiles$II),rownames(covariates)),]

	sampleNames <- colnames(betaQuantiles[[1]])
	slideNames     <- substr(sampleNames,1,10)
	arrayNames    <- substr(sampleNames,12,17)
	plateNames <- substr(sampleNames,1,7)
	 
	
	designInfo <- data.frame(sampleNames = sampleNames,
												slideNames  = slideNames,
												arrayNames = arrayNames,
												plateNames = plateNames)
	
	#### Ordering:
	designInfo <- designInfo[order(plateNames,slideNames,arrayNames), ]
	order <- match(designInfo$sampleNames, 
				                colnames(betaQuantiles[[1]]))


    sampleNames <- designInfo$sampleNames
	slideNames     <- designInfo$slideNames
	arrayNames    <- designInfo$arrayNames
    plateNames    <- designInfo$plateNames
	
	plate <- as.numeric(as.factor(plateNames))
	slides <- unique(slideNames)
	names(slideNames) <- slideNames
	names(slides) <- slides

	for (i in 1:3){
		betaQuantiles[[i]]         <- betaQuantiles[[i]][,order]
		mQuantiles[[i]]             <- mQuantiles[[i]][,order]
		methQuantiles[[i]]        <- methQuantiles[[i]][,order]
		unmethQuantiles[[i]]    <- unmethQuantiles[[i]][,order]
	}
	
	for (i in 1:12){
		greenControls[[i]]    <- greenControls[[i]][,order]
		redControls[[i]]        <- redControls[[i]][,order]
	}

   	XYMedians$medianXU = XYMedians$medianXU[order]
	XYMedians$medianXM = XYMedians$medianXM[order]
	XYMedians$medianYU = XYMedians$medianYU[order]
	XYMedians$medianYM = XYMedians$medianYM[order]
	
   covariates <- covariates[order,]
   
    controlNames <- names(greenControls)





shinyUI(pageWithSidebar(
	
	

  	
###########################  ---  Header ------------ ##############
  	



 
headerPanel(HTML("<p style=\"color:#000000;font-family:\"Times New Roman\",Georgia,Serif\">
         shiny<span style=\"color:#E56717\">M</span><span style=\"color:#000000\">ethyl</span></p>")),






  	
###########################  ---  Sidebar ---------- ###############

	
  sidebarPanel(


	wellPanel(
	  HTML("<p><span style=\"color:#336666;font-size:16px\">
			      Quality control exploration</span></p>"),

	selectInput("mOrBeta", "Methylation measure:",
	       list("Beta-value","M-value"),
	           multiple=FALSE, 
	               selected="Beta-value"),
	    selectInput("probeType", "Choose a probe type:", 
			                	choices = c("I Green","I Red","II"),
			                	selected="II"),
	    selectInput("controlType", "Choose a control type:", 
			                	choices = controlNames,selected=controlNames[1]),
			                	
			                	
	    if ("Sample_Plate" %in% colnames(covariates)){
					choices <- colnames(covariates)
					selectInput("phenotype", "Choose a phenotype:",
					       choices,
					           multiple=FALSE, 
					               selected ="Sample_Plate")
		} else {
				 choices <- colnames(covariates)
					selectInput("phenotype", "Choose a phenotype:",
					       choices,
					           multiple=FALSE)},          	
			                	                	
			                	
		selectInput("slides", "Choose slides to explore:",slides,multiple=TRUE),
			             
	
	checkboxInput("mean","Average density by slide"),

	checkboxInput("solidLine","Solid lines", value=TRUE)
	),
	
	
	
	    ##########  -- Sliders       ---######       
		wellPanel(

		    		 sliderInput(inputId = "bandwidth",
		                  label = "Bandwidth Beta-value",
		                  min = 0, max = 0.05, step = 0.001, value = 0.02),
		             sliderInput(inputId = "bandwidth2",
		                  label = "Bandwidth M-value",
		                  min = 0, max = 0.5, step = 0.001, value = 0.35)
	                        
		 )
		


  ),



  mainPanel(
  	tabsetPanel(

  	
###########################  ---  Home   ------------ 
  	
  		
  				tabPanel("Home",
  		HTML("<br>
  		<p style=\"width:500px;text-align:justify\"><span style=\"color:#000000;font-size:16px\">
  		Welcome to
  		 <span style=\"font-weight:bold\">shiny</span><span style=\"color:#E56717;font-weight:bold\">M</span><span style=\"font-weight:bold\">ethyl</span>, 
  		 the interactive R-package for 
  		exploration of methylation data based on shiny. The current version
  		is designed for the Illumina HumanMethylation 450k Array. 
 <br><br>For more information, please visit <span style=\"font-style:italic\">shinyMethyl.com</span><br><br><br></p>")),
  		


  	
###########################  ---  Quality control --- 
  	


	  	tabPanel("Quality control",


				##########  --- First plots ---#######
				
				wellPanel(	
				HTML('<table border=0 width="100%"><tr bgcolor="#f5f5f5"><td>'),

  				div(style="width:100%;max-width:600px;",
						plotOutput("internalControls", clickId = "controlsHover")
			             ), 
				HTML('</td><td>'),
			    plotOutput("rawDensities"),
			    HTML('</td></tr></table>')
			    ),
			    
			    

			    ##########  --- Second plots ---######
			    wellPanel(	
				HTML('<table border=0 width="100%"><tr bgcolor="#f5f5f5"><td>'),

  				div(style="width:80%;max-width:600px;",
						plotOutput("medianChannels")
			             ), 
				HTML('</td><td>'),

						
						
						#plotOutput("normalizedDensities")
			          #    ), 
			
			    HTML('</td></tr></table>')
			    )
  
		),
		

######################   ----   Array Design  -------- 


tabPanel("Array Design",



HTML("<br>
  		<p style=\"width:800px;text-align:justify\"><span style=\"color:#000000;font-size:16px\">
  		The Illumina 450k samples are divided into slides. A slide contains 12 samples (6 by 2 grid) an a plate contains 8 slides (96 samples). The plot below shows the allocation of the samples to the plates and the coloring allows the user to judge if the design is well-balanced for different phenotype covariates.
  		</span></p><br><br>") ,
  		
  		
				wellPanel(
			  

			       div(style="max-height:800px;",
						plotOutput("arrayDesign")
			    ),
			    
			     div(style="max-height:800px;",
						plotOutput("arrayDesignLegend")
			    )
				)

	

		),
		

######################   ----   Gender clustering  --------  



		tabPanel("Gender clustering",
		
	

   
 HTML("<br>
  		<p style=\"width:800px;text-align:justify\"><span style=\"color:#000000;font-size:16px\">
  		By comparing the median total intensity of the Y-chromosome-mapped probes to the median total intensity of the X-chromosome-mapped probes, where the total intensity is the sum of the methylated and unmethylated signals, it is possible to predict the gender of the sample by looking at the two distinct clusters of intensities. See the minfi function <span style=\"font-style:italic\">getSex()</span>.  		</span></p><br><br>") ,


		
		plotOutput("genderClustering", clickId = "genderCutoff"),
		

  sidebarPanel(

   
    HTML("
<p style=\"color:#000000;font-size:17px\">A. Click on the plot to choose the cutoff</span></p>
"),   
                
        


    
       HTML("
<p style=\"color:#000000;font-size:17px\">B. Save the predicted gender in a csv file:</span></p>
"),   
   downloadLink("downloadClusters","predictedGender.csv")

  ),
  
  sidebarPanel(

             
      HTML("
<p style=\"color:#000000;font-size:17px\">List of samples whose predicted gender do not agree on the gender provided in the phenotype data:</span></p>
"),

  verbatimTextOutput(outputId = "diffPrint")
                
        


  )

    
		),


######################   ----   PCA plot  --------  



		tabPanel("PC Analysis",
		
		plotOutput("pcaPlot"),
		
HTML("
<p style=\"color:#000000;font-size:17px\">A. Choose two principal components to visualize: </span></p>
"),   
                

    selectInput("pc1", "PC in x:",
	       seq(1,ncol(betaQuantiles[[1]]),1),
	           multiple=FALSE, 
	               selected=1),
    selectInput("pc2", "PC in y:",
	       seq(1,ncol(betaQuantiles[[1]]),1),
	           multiple=FALSE, 
	               selected=2),	               
	               
	               
  HTML("
<p style=\"color:#000000;font-size:17px\">B. Choose a principal component to explore below:</span></p>
"),   

 selectInput("pcToExplore", "PC:",
	       seq(1,ncol(betaQuantiles[[1]]),1),
	           multiple=FALSE, 
	               selected=1),
 HTML("
<p style=\"color:#000000;font-size:17px\">C. Choose a covariate to regress against the chosen PC:</span></p>
"),   


if (exists("covariates")){
	choices <- colnames(covariates)
	selectInput("covToRegress", "Covariate:",
	       choices,
	           multiple= FALSE)
} else {
 selectInput("covToRegress", "Covariate:",
	       list("Batch"),
	           multiple= FALSE, 
	               selected="Batch")},
	               
	               

 # ),
  


  
         
      HTML("
<p style=\"color:#000000;font-size:17px\">Association with the PC:</span></p>
"), 
                
  verbatimTextOutput(outputId = "modelPrint")
  



		
		

		),



######################   ----   About  --------  

############################################################

		tabPanel("About",
		
		
		
		HTML("<br>
  		<p style=\"width:500px;text-align:justify\"><span style=\"color:#000000;font-size:16px\">
  		<span style = \"font-weight:bold\">About</span><br><br>
  		<span style=\"font-style:italic\">shinyMethyl</span> is a complementary tool to <span style=\"font-style:italic\">minfi</span> for visualizing methylation data from the Illumina 450k array. The application is still in development and will be soon be compatible with methylation data from the Illumina 27k array. 
  		<br><br>
          <span style = \"font-weight:bold\">Acknowledgements</span><br><br>
           The <span style=\"font-style:italic\">shinyMethyl</span> application is based on the package <span style=\"font-style:italic\">minfi</span>, and a large part of the source code is inspired by the work done by the <span style=\"font-style:italic\">minfi</span>'s authors. The gender clustering panel and the lower QC plot of the quality control panel are based on algorithms available in the development version of <span style=\"font-style:italic\">minfi</span> at the Bioconductor project.
          <br><br> 
          <span style=\"font-style:italic\">shinyMethyl</span> is currently developed at the Johns Hopkins Department of Biostatistics, under the supervision of Kasper D. Hansen. Many thanks to Elizabeth M. Sweeney and John Muschelli for their help and for providing precious feedbacks. 
          </br></br>
          <span style = \"font-weight:bold\">Author:</span>  Jean-Philippe Fortin
          <span style=\"color:#FFFFFF\">aa</span>   
          <span style=\"font-style:italic\">  jfortinbiostats.com</span><br><br>
             <span style = \"font-weight:bold\">Webpage:</span>  <span style=\"font-style:italic\">shinymethyl.com</span><br><br>
          
  		</span></p>")
  		


    
		)
########################################################################################		
	
	)

  )
  
  
))

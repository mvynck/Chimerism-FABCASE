library(rhandsontable)
library(shinydashboard)
library(shinyBS)
library(shiny)
library(ggplot2)
library(rmarkdown)

numSamples<-20
numInf<-3
nSets = 1
ui<-shinyUI(dashboardPage(
    skin="red",
	dashboardHeader(title = "FABCASE"),
	dashboardSidebar(
	sidebarMenu(
    	menuItem("Calculate", tabName = "data", icon = icon("line-chart")),
    	menuItem("Help", tabName = "help", icon = icon("question-circle"),
    		#menuSubItem("Citation", tabName = "cite"),
   	 		menuSubItem("FAQ", tabName = "faq"))
    )

),
dashboardBody(
	tabItems(
    tabItem(tabName = "data",
		fluidRow(
			box(
				#HTML("<hr>"),
				title="Configuration",
				width=8,
				height=450,
				sliderInput("nsample", "Number of markers screened:", 1, 150, numSamples),
				sliderInput("ninf", "Minimum number of informative markers:", 1, 10, 3),
				numericInput("nallele", "Number of alleles measured:", 100, min = 10, max = 10000),
				radioButtons("eligible", "Marker types to consider:",
				             c("Informative (donor homozygous and not identical to recipient)" = "i",
				               "Informative and potentially informative (donor and recipient not identical)" = "i-ii"))
			),
			box(
				title="Example data",
				width=4,
				height=450,
    			actionButton("exampleLoad", "Load example dataset"),
				
			)

	   ),
	   fluidRow(
			box(style='overflow-x: scroll;overflow-y: scroll;',
				title="Data",
				width=3,
				footer = "Note: you can copy-paste values which is especially useful for larger number of markers.",
				rHandsontableOutput("hot")
			),
			box(style='overflow-x: scroll;overflow-y: scroll;',
			  title="Results",
			  width=9,
			  tableOutput('tableInf')
			)
		)
	),
    tabItem(tabName = "faq",
    	box(
			h4("Where can I find the underlying methodology?"),
			p("The full text of the paper is available at ...."),
 	    h4("What does informative and potentially informative mean?"),
		 	p("Markers that are heterozygous in the donor and where the genotype is
  		 	different from the genotype of the recipient are termed 'potentially
  		 	informative' (also: 'type-II'). Such markers often suffer biases and are less precise,
  		 	rendering them less useful than markers that are homozygous in the donor
  		 	and where the genotype is not identical to that of the recipient
  		 	('informative markers', 'type-I')."),tags$a(href="https://www.sciencedirect.com/science/article/pii/S1525157821001719", "See Vynck et al. (2021) Journal of Molecular
  		 	Diagnostics for examples", target="_blank"),
		 	h4("Questions? Feel free to contact me at Matthijs.Vynck@UGent.be.")
      	)
	)
  )
)
)
)

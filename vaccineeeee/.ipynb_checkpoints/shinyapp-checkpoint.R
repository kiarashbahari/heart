#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(tidyverse)

#user interface
ui <- fluidPage(pageWithSidebar( 
  
  headerPanel = headerPanel("BIOL336: Haploid model of selection"),
  
  sidebarPanel(
    
    HTML("<p style='font-size:14px'><B>Frequency of allele A over time (p).</B>"),

    sliderInput(inputId = "wA", label = "wA", value = 1.1, 
                min = 0, max = 2, step = 0.01),
    sliderInput(inputId = "wa", label = "wa", value = 1, 
                min = 0, max = 2, step = 0.01),
    
    sliderInput(inputId = "p_0", label = "Initial p", value = 0.1, 
                min = 0, max = 1, step = 0.01),
    
    sliderInput(inputId = "gen", label = "Number of generations", value = 100, 
                min = 5, max = 1000, step = 5),
  
    HTML("<p style='font-size:12px'>Alternatively, <A href='https://shiney.zoology.ubc.ca/otto/HaploidSelectionInput'>click here</a> if you would like to input numbers directly.<P>"), 

    HTML("<p style='font-size:12px'>Take me to:
         <UL><LI style='font-size:12px'><A href='https://shiney.zoology.ubc.ca/otto/HaploidSelection/'>Haploid selection</a> (<A href='https://shiney.zoology.ubc.ca/otto/HaploidMutationSelection/'>with mutation</a>)
         <LI style='font-size:12px'><A href='https://shiney.zoology.ubc.ca/otto/DiploidSelection/'>Diploid selection</a> (<A href='https://shiney.zoology.ubc.ca/otto/DiploidMutationSelection/'>with mutation</a>)
         <LI style='font-size:12px'><A href='https://shiney.zoology.ubc.ca/otto/DiploidDriftSelection/'>Diploid selection with drift</a>
         </UL><P>"),
    HTML("<p style='font-size:8px'>Download R <A href='https://www.zoology.ubc.ca/~otto/Research/ShineyAppsForPopGen.zip'>source code</a>.<BR>
          Modified from: Copyright (c) 2017 Silas Tittes, MIT License, <A href='https://github.com/silastittes/shiny_popgen'>https://github.com/silastittes/shiny_popgen</a></p>")
    
  ), 
  
  mainPanel =  mainPanel(
    plotOutput(outputId = 'viz'),
    
    selectInput("select", label = "Plot options", 
                choices = list("Mean fitness by p" = 1, "Change in p by p" = 2,
                               "Allele frequency over time" = 3, "Genotypes over time" = 4, "Mean fitness over time" = 5), selected = 3)
    
    
  )
))

#back end code and response to user input
server <- function(input, output){
  
  output$viz <- renderPlot({
    
    p <- seq(0, 1, length.out = 1000)
    #parameters
    wA = input$wA
    wa = input$wa
    p_0 = input$p_0
    gen = (input$gen)+1
    t = seq(0, gen, length.out = 1000)
    
    #if(input$w_plot){
    if(input$select == 1){
      W <- p*wA + (1-p)*wa
      data.frame(p=p, W=W) %>%
        ggplot(aes(x = p, y = W)) +
        geom_line(color="forestgreen", size=2) +
        xlab("Mean Fitness") +
        xlab("Allele frequency, p") +
#        xlim(0, 1) +
        scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) + 
        scale_y_continuous(expand = c(0, 0), limits = c(min(wa,wA)-0.001, max(wa,wA)+0.001)) + 
#        theme(panel.background =  
#                element_rect(fill =  rgb(30, 144, 255, 25, 
#                                         maxColorValue = 255)),
#              text = element_text(size=16, family= "Times"))
        theme_classic()+theme(text = element_text(size=20))
      
      #} else if(input$delta_plot){
    } else if(input$select == 2){
      W <- p*wA + (1-p)*wa
      delta_p <- p*wA/W - p
      setlim <- (max(wa,wA)-min(wa,wA))/(wA + wa)
      setlimround <- round(10*setlim)/10+0.1
      data.frame(p=p, delta_p=delta_p) %>%
        ggplot(aes(x = p, y = delta_p)) +
        geom_line(color="firebrick", size=2) +
        geom_hline(yintercept = 0, lty = 2) +
        ylab(expression(paste(Delta,p))) +
        xlab("Allele frequency, p") +
#        xlim(0, 1) +
        scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) + 
        scale_y_continuous(expand = c(0, 0), limits = c(-setlimround,setlimround)) + 
#        theme(panel.background =  
#                element_rect(fill =  rgb(30, 144, 255, 25, 
#                                         maxColorValue = 255)),
#              text = element_text(size=16, family= "Times"))
        theme_classic()+theme(text = element_text(size=20))
      
      #} else if(input$time_plot){
    } else if(input$select == 3){
      
      p_t <- rep(NA, gen)
      p_t[1] <- p_0
      for(i in 2:gen){
        W <- p_t[i-1]*wA + (1-p_t[i-1])*wa
        p_t[i] <- (p_t[i-1] * wA) / W
      }
      
      data.frame(t = 1:gen, p_t = p_t) %>%
        ggplot(aes(x = t-1, y = p_t)) +
        geom_point(color="firebrick", size=2) +
        xlab("Generation") +
        ylab("Allele frequency, p") +
        annotate("text", x = 0.9*gen, y = 0.14, label = "Final p")+
        annotate("text", x = 0.9*gen, y = 0.1, label = round(p_t[gen],5))+
#        ylim(0, 1) + 
        scale_x_continuous(expand = c(0, 0), limits = c(0, input$gen+2)) + 
        scale_y_continuous(expand = c(0, 0), limits = c(0, 1.01)) + 
#        theme(panel.background =  
#                element_rect(fill =  rgb(30, 144, 255, 25, 
#                                         maxColorValue = 255)),
#              text = element_text(size=16, family= "Times"))
        theme_classic()+theme(text = element_text(size=20))
    } else if(input$select == 4){
      
      p_t <- rep(NA, gen)
      p_t[1] <- p_0
      for(i in 2:gen){
        W <- p_t[i-1]*wA + (1-p_t[i-1])*wa
        p_t[i] <- (p_t[i-1] * wA) / W
      }
      
      data.frame(t = 1:gen, p_t = p_t) %>%
        ggplot(aes(x = t-1)) +
        scale_color_manual(values = c("genotype A" = "firebrick", "genotype a" = "gold")) +
        geom_point(aes(y = (1-p_t), color="genotype a"), size=1) +
        geom_point(aes(y = p_t, color="genotype A"), size=1) +
        labs(x= "Generation",y="Allele frequency, p",
             color = "Legend") +
        scale_x_continuous(expand = c(0, 0), limits = c(0, input$gen+2)) + 
        scale_y_continuous(expand = c(0, 0), limits = c(0, 1.01)) + 
        #        theme(panel.background =  
        #                element_rect(fill =  rgb(30, 144, 255, 25, 
        #                                         maxColorValue = 255)),
        #              text = element_text(size=16, family= "Times"))
        theme_classic()+theme(text = element_text(size=20))
    } else if(input$select == 5){
      p_t <- rep(NA, gen)
      W_t <- rep(NA, gen)
      p_t[1] <- p_0
      W_t[1] <- p_0*wA + (1-p_0)*wa
      for(i in 2:gen){
        W_t[i] <- p_t[i-1]*wA + (1-p_t[i-1])*wa
        p_t[i] <- (p_t[i-1] * wA) / W_t[i]
      }
      
      data.frame(t = 1:gen, W_t = W_t) %>%
        ggplot(aes(x = t-1, y = W_t)) +
        geom_point(color="forestgreen", size=2) +
        xlab("Generation") +
        ylab("Mean fitness") +
        #        ylim(0, 1) + 
        scale_x_continuous(expand = c(0, 0), limits = c(0, input$gen+2)) + 
        scale_y_continuous(expand = c(0, 0), limits = c(min(wa,wA)-0.001, max(wa,wA)+0.001)) + 
        #        theme(panel.background =  
        #                element_rect(fill =  rgb(30, 144, 255, 25, 
        #                                         maxColorValue = 255)),
        #              text = element_text(size=16, family= "Times"))
        theme_classic()+theme(text = element_text(size=20))
    }
  })
}
# Run the application 
shinyApp(ui = ui, server = server)
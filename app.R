## Author: Mabel Chang
## mchang15@bu.edu
## BF591
## Final Project

library(shiny)
library(DT)
library(tidyverse)
library(colourpicker)
library(ggplot2)
library(reshape2)
library(pheatmap)
library('fgsea')

#define UI for application
ui <- fluidPage(
  titlePanel("BF591 Final Project - Huntington’s Disease"),
  
  #create panel of tabs
  tabsetPanel(
####Samples tab
    tabPanel("Samples",
             sidebarLayout(
               sidebarPanel(
                 #input for the sample tab
                  fileInput(
                    "samples_tab_file",
                    label = "Load metadata csv file",
                    accept = c("text/csv",
                                "text/comma-separated-values,text/plain",
                                ".csv"))
                ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Summary", 
                            textOutput("samples_summary_dim"),
                            DTOutput("samples_summary_table")),
                   tabPanel("Table", 
                            DTOutput("samples_table_table")),
                   tabPanel("Plots", 
                            sidebarLayout(
                              sidebarPanel(
                                radioButtons("samples_plots_yaxis", "y-variable",
                                  choices = c("age_of_death", "AvgSpotLen", "Bases", "Bytes", "MBases", "MBytes", "mrna.seq_reads","pmi","RIN","age_of_onset","cag","duration","h.v_cortical_score","h.v_striatal_score","vonsattel_grade"),
                                  selected = "age_of_death"),
                                submitButton("Submit", width = "100%")
                              ),
                              mainPanel(
                                plotOutput("samples_plots_histogram")
                              )
                            )
                    )
                 )
               )
             )
    ),
####Counts tab
    tabPanel("Counts", 
             sidebarLayout(
               sidebarPanel(
                 #input for counts tab
                 fileInput(
                   "counts_tab_file",
                   label = "Load normalized counts csv file",
                   accept = c("text/csv",
                              "text/comma-separated-values,text/plain",
                              ".csv")),
                 #slider variance for counts tab
                 sliderInput("counts_variance_slider", "include genes with at least X percentile of variance:",
                             min = 0, max = 100, value = 80, step = 1),
                 #slider non-zero for counts tab
                 sliderInput("counts_nonzero_slider", "include genes with at least X samples that are non-zero:",
                             min = 0, max = 69, value = 60, step = 1),
                 submitButton("Submit", width = "100%")
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Summary", 
                            DTOutput("counts_summary_table")),
                   tabPanel("Diagnostics", 
                            plotOutput("counts_diagnostic_variance_plot"),
                            plotOutput("counts_diagnostic_zero_plot")),
                   tabPanel("Heatmap", 
                            plotOutput("counts_heatmap")),
                   tabPanel("PCA",
                            sidebarLayout(
                              sidebarPanel(
                                sliderInput("counts_PCA_xaxis_slider", "X-axis Principle Component",
                                            min = 1, max = 69, value = 1, step = 1),
                                sliderInput("counts_PCA_yaxis_slider", "Y-axis Principle Component",
                                            min = 1, max = 69, value = 2, step = 1),
                                submitButton("Submit", width = "100%")
                              ),
                              mainPanel(
                                plotOutput("counts_PCA_plot")
                              )
                            )
                          )
                 )
               )
             )
    ),
####DE tab
    tabPanel("DE",
             sidebarLayout(
               sidebarPanel(
                 #input for DE tab
                 fileInput(
                   "de_tab_file",
                   label = "Load normalized counts csv file",
                   accept = c("text/csv",
                              "text/comma-separated-values,text/plain",
                              ".csv"))
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Input",
                            DTOutput("de_input_table")),
                   tabPanel("Result",
                            sidebarLayout(
                              sidebarPanel(
                                radioButtons("de_result_xaxis", "Column for the x-axis",
                                  choices = c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"),
                                  selected = "log2FoldChange"),
                                radioButtons("de_result_yaxis", "Column for the y-axis",
                                  choices = c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"),
                                  selected = "padj"),
                                colourInput("de_result_base_color", "Base point color",
                                  value = "black", closeOnClick = T),
                                colourInput("de_result_highlight_color", "Highlight point color",
                                  value = "grey", closeOnClick = T),
                                sliderInput("de_result_magnitude_slider","Magnitude of the p adjusted coloring:",
                                  min = -300, max = 0, value = -150, step = 1),
                                submitButton("Submit", width = "100%")
                              ),
                              mainPanel(
                                tabsetPanel(
                                  tabPanel("Plot",
                                           plotOutput("de_result_plot_volcano")),
                                  tabPanel("Table",
                                           DTOutput("de_result_table"))
                                )
                              )
                            ))
                 )
               )
             )),
    tabPanel("GSEA", 
             sidebarLayout(
               sidebarPanel(
                 #input for GSEA tab
                 fileInput("gsea_tab_file",
                   label = "Load preprocessed csv file",
                   accept = c("text/csv",
                              "text/comma-separated-values,text/plain",
                              ".csv"))
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Barplot",
                            sidebarLayout(
                              sidebarPanel(
                                sliderInput("gsea_barplot_slider","Number of top pathways (pos and neg) to plot:",
                                            min = 1, max = 50, value = 10, step = 1),
                                submitButton("Plot", width = "100%")
                              ),
                              mainPanel(
                                plotOutput("gsea_barplot_plot")
                              )
                            )),
                   tabPanel("Table",
                            sidebarLayout(
                              sidebarPanel(
                                sliderInput("gsea_table_pvalue_slider", "Adjusted p-value threshold:",
                                            min = 0, max = 1, value = 0.05, step = 0.00005),
                                radioButtons("gsea_table_filter_button", "Pathways with +/- NES",
                                             choices = c("All", "Positive", "Negative"),
                                             selected = "All"),
                                submitButton("Submit", width = "100%"),
                                downloadButton("gsea_table_download_button", "Download")
                              ),
                              mainPanel(
                                DTOutput("gsea_table_table")
                              )
                            )),
                   tabPanel("Scatterplot",
                            sidebarLayout(
                              sidebarPanel(
                                sliderInput("gsea_scatterplot_pvalue_slider", "Adjusted p-value threshold:",
                                            min=0, max = 1, value = 0.05, step = 0.00005),
                                submitButton("Submit", width = "100%")
                              ),
                              mainPanel(
                                plotOutput("gsea_scatterplot_plot")
                              )
                            ))
                )
              )
            )
          )
  )
)


server <- function(input, output, session) {
  #increase the size of datasets that can be used to up to 30mb
  options(shiny.maxRequestSize=30*1024^2)
  
###Functions
 ##Summary tab
  gen_samples_summary_table <- function(data) {
    #get column names and data types
    column_info <- data.frame(
      Column_Name = names(data),
      Type = sapply(data, class)
    )
    #get mean of each column
    means <- sapply(data, mean, na.rm = TRUE)
    #combine information into a summary table
    summary_table <- cbind(column_info, Mean = means)
    #remove rownames
    rownames(summary_table) <- NULL
    return(summary_table)
  }

 ##Counts tab
  gen_counts_summary_table <- function(data, varslider, zeroslider) {
    #dimensions of sample
    number_of_columns <- ncol(data) - 1
    number_of_rows <- nrow(data)
    #add a column for the number of non-zero entries
    data$NonZeroCount <- apply(data[, -1], 1, function(row) sum(row != 0))
    #add a column for the row variance
    data$RowVariance <- apply(data[, -c(1, ncol(data))], 1, var)
    #filter the data
    data <- data %>%
      filter(RowVariance >= quantile(RowVariance, prob = as.numeric(varslider)/100),
             NonZeroCount >= as.numeric(zeroslider))
    #get number of passing rows
    number_of_passing_rows <- nrow(data)
    #find percentages for output
    formatted_pass_percentage <- sprintf("%.1f%%", (number_of_passing_rows/number_of_rows) * 100)
    formatted_nonpass_percentage <- sprintf("%.1f%%", ((number_of_rows-number_of_passing_rows)/number_of_rows) * 100)
    #output table
    output_table <- data.frame(
      num_samples = c("Number of samples", number_of_columns),
      num_genes = c("Number of genes", number_of_rows),
      num_passing = c("Number of passing genes", number_of_passing_rows),
      per_passing = c("Percent of passing genes", formatted_pass_percentage),
      non_passing = c("Number of non-passing genes", number_of_rows-number_of_passing_rows),
      per_non_passing = c("Percent of non-passing genes", formatted_nonpass_percentage)
    )
    output_table <- t(output_table)
    rownames(output_table) <- NULL
    colnames(output_table) <- c("metric", "value")
    return(output_table)
  }
  
  melt_data_heatmap <- function(data, varslider, zeroslider) {
    #add a column for the number of non-zero entries
    data$NonZeroCount <- apply(data[, -1], 1, function(row) sum(row != 0))
    #add a column for the row variance
    data$RowVariance <- apply(data[, -c(1, ncol(data))], 1, var)
    #filter the data
    data <- data %>%
      filter(RowVariance >= quantile(RowVariance, prob = varslider/100),
             NonZeroCount >= zeroslider)
    #remove the last two columns
    data <- data[, -c((ncol(data) - 1):ncol(data))]
    #set row names
    rownames(data) <- data$X
    data <- data[, -1]
    #convert to matrix
    data <- as.matrix(data)
    #melt the data
    data_melt <- melt(data)
    return(data_melt)
  }
  
  gen_counts_pca_results <- function(data, varfilter, zerofilter){
    #add a column for the number of non-zero entries
    data$NonZeroCount <- apply(data[, -1], 1, function(row) sum(row != 0))
    #add a column for the row variance
    data$RowVariance <- apply(data[, -c(1, ncol(data))], 1, var)
    #filter the data
    data <- data %>%
      filter(RowVariance >= quantile(RowVariance, prob = varfilter/100),
             NonZeroCount >= zerofilter)
    #remove last two columns
    data <- data[, -c((ncol(data) - 1):ncol(data))]
    #set gene names as row names
    rownames(data) <- data$X
    data <- data[, -1]
    #perform PCA
    pca_result <- prcomp(data, scale. = TRUE)
    return(pca_result$x)
  }
  
 ##DE tab
  gen_de_result_plot_volcano <- function(data, xaxis, yaxis, magnitude, basecolor, highlightcolor) {
    plot <- 
      ggplot(data, aes_string(x=xaxis, y=paste0("-log10(", yaxis, ")"))) +
      geom_point(aes(color= ifelse(-log10(get(yaxis)) > (-as.numeric(magnitude)), "sig", "not_sig"))) +
      scale_color_manual(values = c("sig" = highlightcolor, "not_sig" = basecolor)) +
      theme_bw() +
      labs(x=xaxis, 
           y=paste0('-log10(', yaxis, ')'), 
           color=paste0(yaxis,'< 10^(', magnitude, ')')) +
      theme(legend.position = "bottom")
    return(plot)
  }
  
  gen_de_result_table <- function(data, magnitude) {
    #filter data frame based on p-adjusted values
    filtered_data <- data[data$padj < 10^magnitude, ]
    #remove all rows with NA
    filtered_data <- na.omit(filtered_data)
    # #format p-value and p-adjusted value columns
    # filtered_data$pvalue <- formatC(filtered_data$pvalue, format = "e", digits = 3)
    # filtered_data$padj <- formatC(filtered_data$padj, format = "e", digits = 3)
    return(filtered_data)
  }

 ##GSEA Tab
  gen_gsea_barplot <- function(fgsea_results, num_paths){
    #filter the top positive
    top_positive <- fgsea_results %>%
                    filter(NES > 0) %>%
                    arrange(desc(NES)) %>%
                    mutate(pos_neg = "positive") %>%
                    head(num_paths)
    #filter the top negative
    top_negative <- fgsea_results %>%
                    filter(NES < 0) %>%
                    arrange(NES) %>%
                    mutate(pos_neg = "negative") %>%
                    head(num_paths)
    #combine the positive and negative
    top_pathways_combined <- rbind(top_positive, top_negative)
    #create a ggplot bar chart
    top_pathways_plot <- ggplot(top_pathways_combined, aes(x = NES, y = reorder(pathway, NES), fill = NES)) +
      geom_bar(stat = "identity") +
      labs(
        title = "fgsea results for MSigDB C2 Curated Gene Sets",
        x = "Normalized Enrichment Score (NES)") +
      theme_minimal() +
      theme(axis.text.y = element_text(angle = 45, hjust = 1))
    return(top_pathways_plot)
  }
  
  
  
  
  
  
###Elements inputs
 ##Summary tab
  #load in data for samples tab
  samples_load_data <- reactive({
    #make sure a file is uploaded
    req(input$samples_tab_file)
    #read in file
    df <- read.csv(input$samples_tab_file$datapath)
    #remove first column
    df <- df[, -1]
  })
  
  #filter selected column for histogram
  samples_plots_yaxis_column <- reactive({
    #make sure something is selected
    req(input$samples_plots_yaxis)
    #filter dataframe to only the chosen column
    data <- samples_load_data()[, input$samples_plots_yaxis]
    data <- na.omit(data)
  })
  
 ##Counts Tab
  #load in data for counts tab
  counts_load_data <- reactive({
    #make sure a file is uploaded
    req(input$counts_tab_file)
    #read in file
    df <- read.csv(input$counts_tab_file$datapath)
    #remove first column
    df <- df[, -1]
  })
  
 ##DE Tab
  #load in data for DE tab
  de_load_data <- reactive({
    #make sure a file is uploaded
    req(input$de_tab_file)
    #read in file
    df <- read.csv(input$de_tab_file$datapath)
    #remove last two columns
    df <- df[, -1]
  })
  
 ##GSEA tab
  #load in data for GSEA tab
  gsea_load_data <- reactive({
    #make sure a file is uploaded
    req(input$gsea_tab_file)
    #read in file
    df <- read.csv(input$gsea_tab_file$datapath)
    # #remove last two columns
    # df <- df[, -c(ncol(df) - 1, ncol(df))]
    # #create a named vector with GeneSymbol as names and log2FoldChange as values
    # rnks <- setNames(df$log2FoldChange, df$symbol.x)
    # #run fgsea
    # c2_pathways <- gmtPathways("data/geneset.gmt")
    # fgsea_results <- fgsea(c2_pathways, 
    #                        rnks, 
    #                        minSize=15, 
    #                        maxSize=500) %>% 
    #                   as_tibble()
  })
  
  #reactive slider for filtering table
  gsea_table_pvalue_filter <- reactive({
    req(input$gsea_table_pvalue_slider)
    data <- gsea_load_data()
    #filter based on adjusted p-value
    data <- data[data$padj <= input$gsea_table_pvalue_slider, ]
    nes_type <- input$gsea_table_filter_button
    #filter based on NES type
    if (nes_type == "Positive") {
      data <- data[data$NES > 0, ]
    } else if (nes_type == "Negative") {
      data <- data[data$NES < 0, ]
    } else if (nes_type == "All") {
      data <- data
    }
  })
  
  gsea_scatterplot_catgeorize <- reactive({
    req(input$gsea_scatterplot_pvalue_slider)
    data <- gsea_load_data()
    threshold <- input$gsea_scatterplot_pvalue_slider
    data %>%
      mutate(Color = ifelse(padj > threshold, "Above", "Below"))
  })
  
###Elements output
 ##Summary tab
  #generate dimensions
  output$samples_summary_dim <- renderText({
    #get the dimensions of the metadata
    dimensions <- dim(samples_load_data())
    #print the number of rows and columns
    paste("Number of rows:", dimensions[1], " Number of columns:", dimensions[2])
  })
  
  #generate summary table
  output$samples_summary_table <- renderDT({
    datatable(gen_samples_summary_table(samples_load_data()))
  })
  
  #generate sortable table
  output$samples_table_table <- renderDT({
    datatable(samples_load_data(), options = list(ordering = TRUE))
  })
  
  #generate histogram
  output$samples_plots_histogram <- renderPlot({
    hist(samples_plots_yaxis_column(), main = paste("Histogram of", input$samples_plots_yaxis),
         xlab= "Counts", ylab = input$samples_plots_yaxis, col = "grey", border = "white")
  })

 ##Counts tab
  #generate summary of variance table
  output$counts_summary_table <- renderDT({
    datatable(gen_counts_summary_table(counts_load_data(), input$counts_variance_slider, input$counts_nonzero_slider))
  })
  
  #generate diagnostic table for variance
  output$counts_diagnostic_variance_plot <- renderPlot({
    #get the data
    data <- counts_load_data()
    #add a column for the row variance
    data$RowVariance <- apply(data[, -1], 1, var)
    #calculate median sample values
    data$Median <- apply(data[, -c(1, ncol(data))], 1, median)
    #create a color vector based on percentile threshold
    color_vector <- ifelse(data$RowVariance >= (quantile(data$RowVariance, probs = input$counts_variance_slider/100, na.rm = TRUE)), "black", "grey")
    #create plot
    plot(
      rank(data$Median), log(data$RowVariance),
      col = color_vector,
      pch = 16,
      xlab = "Rank of Median Sample Values",
      ylab = "Log of Variances",
      main = "Scatter Plot of Log Variances vs. Ranked Median Sample Values",
    )
    #add a legend
    legend("topright", legend = c("Variances above threshold", "Variances below threshold"), fill = c("black", "grey"))
  })
  
  #generate diagnostic table for zeros
  output$counts_diagnostic_zero_plot <- renderPlot({
    #get the data
    data <- counts_load_data()
    #add a column for the row variance
    data$ZeroCount <- apply(data[, -1], 1, function(row) sum(row == 0))
    #calculate median sample values
    data$Median <- apply(data[, -c(1, ncol(data))], 1, median)
    #create a color vector based on percentile threshold
    color_vector <- ifelse(data$ZeroCount <= (69-input$counts_nonzero_slider), "black", "grey")
    #create plot
    plot(
      rank(data$Median), data$ZeroCount,
      col = color_vector,
      xlab = "Rank of Median Sample Values",
      ylab = "Number of Zeros",
      main = "Scatter Plot of Number of Zeros vs. Ranked Median Sample Values",
    )
    #add a legend
    legend("topright", legend = c("Number of non-Zeros above threshold", "Number of non-Zeros below threshold"), fill = c("black", "grey"))
  })
  
  #generate the counts heatmap
  output$counts_heatmap <- renderPlot({
    # #plot the heatmap
    ggplot(melt_data_heatmap(counts_load_data(),input$counts_variance_slider, input$counts_nonzero_slider), aes(x = Var1, y = Var2, fill = log(value))) +
      geom_tile() +
      scale_fill_gradient(low = "blue", high = "red") +
      labs(title = "Heatmap of Log-transformed Counts", x = "Genes", y = "Samples") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 1))
    # pheatmap(melt_data_heatmap(counts_load_data(),input$counts_variance_slider, input$counts_nonzero_slider))
  })
  
  #generate the counts PCA
  output$counts_PCA_plot <- renderPlot({
    # plot the principal components
    plot(gen_counts_pca_results(counts_load_data(), input$counts_variance_slider, input$counts_nonzero_slider)[, input$counts_PCA_xaxis_slider],
         gen_counts_pca_results(counts_load_data(), input$counts_variance_slider, input$counts_nonzero_slider)[, input$counts_PCA_yaxis_slider],
         pch = 16, col = "black", main = "PCA Plot",
         xlab = paste("Principal Component", input$counts_PCA_xaxis_slider),
         ylab = paste("Principal Component", input$counts_PCA_yaxis_slider))
  })
  
 ##DE tab
  output$de_input_table <- renderDT({
    datatable(de_load_data(), options = list(ordering = TRUE))
  })

  output$de_result_plot_volcano <- renderPlot({
    p <-gen_de_result_plot_volcano(de_load_data(), input$de_result_xaxis, input$de_result_yaxis, 
                                   input$de_result_magnitude_slider, 
                                   input$de_result_base_color, input$de_result_highlight_color)
    return(p)
  }, height = 700)
  
  output$de_result_table <- renderDT(
    datatable(gen_de_result_table(de_load_data(), input$de_result_magnitude_slider))
  )
  
 ##GSEA tab
  output$gsea_barplot_plot <- renderPlot({
    gen_gsea_barplot(gsea_load_data(), as.numeric(input$gsea_barplot_slider))
  })
  
  output$gsea_table_table <- renderDT({
    datatable(gsea_table_pvalue_filter())
  })
  
  output$gsea_table_download_button <- downloadHandler(
    filename = paste0(input$gsea_tab_file, ".csv", sep = "\t"),
    content = function(file) {
      write.csv(gsea_table_pvalue_filter(), file)
    }
  )
  
  # observe({
  #   if (input$gsea_table_download_button >0){
  #     filename <- paste0(Sys.Date(), "gsea_download.csv")
  #     write.csv(gsea_table_pvalue_filter(), filename)
  #   }
  # })
  
  output$gsea_scatterplot_plot <- renderPlot({
    ggplot(gsea_scatterplot_catgeorize(), aes(x = NES, y = -log10(padj), color=Color)) +
      geom_point() +
      scale_color_manual(values = c("Above" = "black", "Below" = "grey")) +
      labs(title = "Scatterplot of NES vs -log10 Adjusted P-Value",
           x = "NES", y = "-log10(Adjusted P-Value)") +
      theme_minimal()
  })
}

# Run the application
shinyApp(ui = ui, server = server)

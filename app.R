library(shiny)
library(shinyjs)
library(tidyverse)

options(shiny.maxRequestSize = 500 * 1024^2)




ui <- fluidPage(
  titlePanel("Interactive Gene Counts Viewer"),
  helpText("(The uploaded RDS must be a data frame with rownames as gene symbols and columns as sample names.)"),
  sidebarLayout(
    sidebarPanel(
      fileInput("counts_file", "Upload RDS of Counts", accept = ".rds"),
      uiOutput("gene_selector"),
      uiOutput("condition_selector"),
      checkboxInput("toggle_paired", "Paired"),
      checkboxInput("show_debug", "Show Debug Output", value = FALSE),
      uiOutput("person_identifier_selector"),
      conditionalPanel(
        condition = "input.show_debug == true",
        verbatimTextOutput("debug_df")
      ),
      actionButton("plot_gene", "Plot Gene by Condition")
    ),
    
    
    mainPanel(
      verbatimTextOutput("col_example"),
      verbatimTextOutput("split_example"),
      verbatimTextOutput("condition_part"),
      verbatimTextOutput("identifier_part"),
      plotOutput("gene_plot")
    )
  )
)

server <- function(input, output, session) {
  
  normCounts_data <- reactiveVal(NULL)
  split_parts <- reactiveVal(NULL)
  
  # Step 1: Upload file and split first column
  observeEvent(input$counts_file, {
    req(input$counts_file)
    
    counts <- readRDS(input$counts_file$datapath)
    counts <- as.data.frame(counts)
    colnames(counts) <- gsub("\\.", "_", colnames(counts))
    colnames(counts) <- gsub("-", "_", colnames(counts))
    
    normCounts_data(counts)
    
    # Take first column name as example
    first_col <- colnames(counts)[1]
    
    output$col_example <- renderPrint({
      cat("Example column name:\n", first_col)
    })
    
    parts <- strsplit(first_col, "_")[[1]]
    split_parts(parts)
    
    output$split_example <- renderPrint({
      cat("Split parts \n")
      print(parts)
    })
    output$gene_selector <- renderUI({
      req(normCounts_data())
      selectInput(
        "gene",
        "Select gene symbol:",
        choices = rownames(normCounts_data()),
        selected = rownames(normCounts_data())[1]
      )
    })
  })
  
  # Step 2: Dynamically select which part(s) represent the condition
  output$condition_selector <- renderUI({
    req(split_parts())
    
    selectInput(
      "condition_indices",
      "Select which split part(s) represent the condition:",
      choices = seq_along(split_parts()),
      selected = 1,
      multiple = TRUE
    )
  })
  
  output$person_identifier_selector <- renderUI({
    req(split_parts(), input$toggle_paired)
    
    if (!isTRUE(input$toggle_paired)) return(NULL)
    
    selectInput(
      "person_identifier_indices",
      "Select which split part(s) uniquely idenitfy person:",
      choices = seq_along(split_parts()),
      selected = 1,
      multiple = TRUE
    )
  })
  
  output$identifier_part <- renderPrint({
    req(split_parts(), input$person_identifier_indices)
    selected_parts <- split_parts()[as.numeric(input$person_identifier_indices)]
    combined_condition <- paste(selected_parts, collapse = "_")
    cat("Selected Person Identifier (recombined):\n", combined_condition)
    
  })
  
  # Step 3: Display selected condition part
  output$condition_part <- renderPrint({
    req(split_parts(), input$condition_indices)
    selected_parts <- split_parts()[as.numeric(input$condition_indices)]
    combined_condition <- paste(selected_parts, collapse = "_")
    cat("Selected condition (recombined):\n", combined_condition)
    
  })
  

  # Step 4: Plot gene counts by condition using base R
  plot_data <- reactive({
    req(normCounts_data(), input$gene, input$condition_indices)
    
    counts <- normCounts_data()
    gene_counts <- counts[input$gene, , drop = FALSE]
    
    split_colnames <- strsplit(colnames(gene_counts), "_")
    
    data.frame(
      Sample = colnames(gene_counts),
      Count = as.numeric(gene_counts[1, ]),
      Condition = sapply(split_colnames, function(x) paste(x[as.numeric(input$condition_indices)], collapse = "_")),
      Person = sapply(split_colnames, function(x) paste(x[as.numeric(input$person_identifier_indices)], collapse = "_"))
      
      
    )
  })
  
  output$debug_df <- renderPrint({
    df = plot_data()
    wide_df = wide_data()
    wide_df_filtered = wide_data_filtered()
    long_df = long_data()
    
    cat("=== LONG FORMAT (plot_data: for unpaired ttest/anova) ===\n")
    print(head(df))
    
    cat("\n=== WIDE FORMAT (paired ttest) ===\n")
    print(head(wide_df))
    
    cat("\n=== WIDE FORMAT FILTERED (paired ttest) ===\n")
    print(head(wide_df_filtered))

    cat("\n=== LONG FORMAT FILTERED (repeated measures anova) ===\n")
    print(head(long_df))
  })
  
  
  
  
  wide_data <- reactive({
    df <- plot_data()  # your long format
    conds <- unique(df$Condition)
    
    reshape(
      df[,c("Person","Condition","Count")],
      idvar   = "Person",
      timevar = "Condition",
      direction = "wide"
    )
  })
  
  wide_data_filtered = reactive({
    df <- plot_data()
    wide_df <- wide_data()
    conds <- unique(df$Condition)

    conds <- unique(df$Condition)
    samples_per_group <- table(df$Condition)

    cond_cols <- paste0("Count.", conds)
    cond_cols <- cond_cols[cond_cols %in% colnames(wide_df)]
    wide_data_filtered <- wide_df[complete.cases(wide_df[, cond_cols]), ]

    wide_data_filtered

  })
  

  
  
  long_data = reactive({
    df = plot_data()
    wide_df_filtered <- wide_data_filtered() 
    
    
    conds <- unique(df$Condition)
    
    if(nrow(wide_df_filtered) == 0) {
      # return an empty long data frame with correct columns
      return(data.frame(Person=character(0),
                        Condition=factor(levels=conds),
                        Count=numeric(0)))
    }


    
    long_data <- reshape(
      wide_df_filtered,
      varying = paste0("Count.", conds),
      v.names = "Count",
      timevar = "Condition",
      times = conds,
      idvar = "Person",
      direction = "long"
    )
    
    rownames(long_data) <- NULL
    long_data$Condition <- factor(long_data$Condition)
    
    
    long_data
  })
  
  
  

  
  output$gene_plot <- renderPlot({
    df <- plot_data() 
    wide_df_filtered = wide_data_filtered()
    long_df = long_data()
    

    
  
    
    samples_per_group <- table(df$Condition)
    
    conds <- unique(df$Condition)
    cond_cols <- paste0("Count.", conds)
    cond_cols <- cond_cols[cond_cols %in% colnames(wide_df_filtered)]

    
    
    msg = ""
    if(length(conds) == 2) {
      pval <- t.test(Count ~ Condition, data = df)$p.value
      msg = "; Unpaired t.test"
    }
    if(length(conds) > 2) {
      pval <- summary(aov(Count ~ Condition, data = df))[[1]][["Pr(>F)"]][1]
      msg = "; Unpaired anova"
    }
    if(length(conds) < 2) {
      pval = NA
      msg = "; Your Selection Results in only 1 Group"
    }
    if(!all(samples_per_group >= 2)) {
      pval = NA
      msg = "; One or More Groups have only 1 Sample"
    }
    if(length(cond_cols) == 2 && isTRUE(input$toggle_paired)){
  if(nrow(wide_df_filtered) < 2){
    pval <- NA
    msg <- "; Not enough samples"
  } else {
    pval <- t.test(wide_df_filtered[[cond_cols[1]]], 
                   wide_df_filtered[[cond_cols[2]]], 
                   paired = TRUE)$p.value
    msg <- "; Paired t.test"
  }
}

    if(length(cond_cols) > 2 && isTRUE(input$toggle_paired)){
      
      # Safety checks first
      if (nrow(long_df) == 0) {
        pval <- NA
        msg <- "; No samples for selected Person/Condition (Your Person identifier is likely in the condition selection)"
      } else if (length(unique(long_df$Person)) < 3) {
        pval <- NA
        msg <- "; Need at least 3 individuals"
      } else if (length(unique(long_df$Condition)) < 2) {
        pval <- NA
        msg <- "; Need at least 2 conditions"
      } else {
        # safe to run repeated measures ANOVA
        long_df$Person <- factor(long_df$Person, levels = unique(long_df$Person))
        long_df$Condition <- factor(long_df$Condition, levels = unique(long_df$Condition))
        long_df$Count <- as.numeric(long_df$Count)
        
        fit <- aov(
          Count ~ Condition + Error(Person / Condition),
          data = long_df
        )
        sm <- summary(fit)
        within_table <- sm[["Error: Within"]][[1]]
        pval <- within_table["Condition", "Pr(>F)"]
        msg <- "; Paired anova"
      }
    }
    
    
    

    
    
    
    # Determine y-axis limits
    ylim_vals <- c(
      min(df$Count, na.rm = TRUE),
      max(df$Count, na.rm = TRUE)
    )
    
    boxplot(Count ~ Condition, data = df,
            col = "skyblue",
            main = paste("Counts for", input$gene),
            ylim = ylim_vals,
            las = 2)
    
    stripchart(Count ~ Condition, data = df,
               vertical = TRUE,
               add = TRUE,
               pch = 16,
               col = "darkblue")
    
    if (isTRUE(input$toggle_paired)) {
      conditions <- levels(factor(df$Condition))
      
      # Get x-axis positions of the boxes
      x_pos <- setNames(seq_along(conditions), conditions)
      
      # Loop over each Person
      for (p in unique(df$Person)) {
        person_data <- df[df$Person == p, ]
        
        # Only draw line if the person has data for all conditions
        if (all(conditions %in% person_data$Condition)) {
          y_vals <- person_data$Count[match(conditions, person_data$Condition)]
          lines(x = x_pos, y = y_vals, col = "red", lty = 2)
        }
      }
    }
    
    mtext(paste("p-value =", signif(pval, 3), msg), side = 3, line = 0.5)
    
  })
  
  
  
  
  
}

shinyApp(ui, server)




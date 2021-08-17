# Load libraries
# Shiny app & layout
library(shiny)
library(bslib)
library(shinycssloaders)
library(shinyWidgets)
library(stringi)

# Error checking
library(randomForest)
library(e1071)

# PCA plot
library(DESeq2)
library(preprocessCore)

# Heatmap
#library(biomaRt)
library(dplyr)
library(pheatmap)
library(RColorBrewer)

# Boxplot
library(ggplot2)
library(reshape)
library(tidyr)
library(forcats)

## Global settings
# Upload size to 128 MB
options(
    shiny.maxRequestSize = 128 * 1024 ^ 2,
    #shiny.error = browser,
    shiny.fullstacktrace = F,
    repos = BiocManager::repositories()
)
Sys.setlocale("LC_NUMERIC", "C")

## Global functions
# Function 'not in',

'%ni%' <- Negate('%in%')

## Shiny user interface
ui <- fluidPage(
    title = "HLCompR",
    theme = bs_theme(version = 4, bootswatch = "lux"),
    tags$head(tags$style(
        HTML(
            "div.tabbable ul.nav { margin-bottom:10px; }
                  .shiny-download-link { float:right }
                 h2 { text-transform: none;} "
        )
    )),
    titlePanel("HLCompR"),
    sidebarLayout(
        sidebarPanel(
            helpText(
                "This application allows you to compare your sequencing data to the training dataset from our paper."
            ),
            a("> Article", href = "https://www.nature.com", target = "_blank"),
            br(),
            br(),
            helpText(
                "Please carefully read the instructions on how to process your sequencing data for analysis."
            ),
            a("> Instructions", href = "https://github.com/iardisasmita/HLCompR", target =
                  "_blank"),
            br(),
            br(),
            helpText("By using this app you agree with the Terms of Usage."),
            actionLink(inputId = "terms_of_usage", label = "> Read Terms of Usage"),
            br(),
            br(),
            helpText(
                "Upload your own coldata and countdata to add them to our analyses, or start by using only our training datasets."
            ),
            hr(),
            helpText(""),
            prettySwitch(
                inputId = "use_training_data_only",
                label = "Use only training dataset",
                value = FALSE
            ),
            hr(),
            uiOutput("user_file_input") %>% withSpinner(color =
                                                            "lightgrey")
        ),
        mainPanel(uiOutput("tabs"))
    )
)

## Server side functionality
server <- function(input, output, session) {
    # First, load training datasets: both coldata and countdata
    train_Data <-
        read.csv("train_count_data.csv", header = TRUE) #read count table
    train_Data$ensembl_gene_name <-
        train_Data$ensmbl_gene_name # Typo in Ibrahim's first column
    train_Data$ensmbl_gene_name <- NULL
    train_colData <-
        read.csv("train_new_colData.csv", header = TRUE) #coldata table
    
    #DEBUG_Data <- read.table(file="user_count_data.csv", sep=",", header = T)
    
    observeEvent(input$terms_of_usage, {
        showModal(modalDialog(
            h3("Terms of Usage"),
            h4("Utrecht Unversity Shiny Server"),
            tags$ol(
                tags$li(
                    "Purpose of the service “utrecht-university.shinyapps.io” is to provide a digital place for trying out, evaluating and/or comparing methods developed by researchers of Utrecht University for the scientific community worldwide. The app and its contents may not be preserved in such a way that it can be cited or can be referenced to."
                ),
                tags$li(
                    "The web application is provided ‘as is’ and ‘as available’ and is without any warranty. Your use of this web application is solely at your own risk."
                ),
                tags$li(
                    "You must ensure that you are lawfully entitled and have full authority to upload data in the web application. The file data must not contain any data which can raise issues relating to abuse, confidentiality, privacy, data protection, licensing, and/or intellectual property. You shall not upload data with any confidential or proprietary information that you desire or are required to keep secret."
                ),
                tags$li("By using this app you agree to be bound by the above terms.")
            ),
            easyClose = TRUE,
            footer = NULL
        ))
    })
    
    train_colData_filtered <- reactive({
        # Check if filters are set - if so, assign variable
        filter_culture_type <-
            ifelse(input$filter_culture_type != "No preference",
                   input$filter_culture_type,
                   NA)
        filter_personalized_medicine <-
            ifelse(
                input$filter_personalized_medicine != "No preference",
                input$filter_personalized_medicine,
                NA
            )
        filter_expandability <-
            ifelse(
                input$filter_expandability != "No preference",
                input$filter_expandability,
                NA
            )
        filter_culture_duration <-
            ifelse(
                input$filter_culture_duration != "No preference",
                input$filter_culture_duration,
                NA
            )
        filter_cell_source <-
            ifelse(input$filter_cell_source != "No preference",
                   input$filter_cell_source,
                   NA)
        filter_viral_transduction <-
            ifelse(
                input$filter_viral_transduction != "No preference",
                input$filter_viral_transduction,
                NA
            )
        
        # Load colData into filter variable
        filter <- train_colData
        
        # Apply filters -> always select Control + selected item
        if (!is.na(filter_culture_type))
            filter <-
            filter %>% dplyr::filter(Culture.Type == filter_culture_type |
                                         Culture.Type == "Control")
        if (!is.na(filter_personalized_medicine))
            filter <-
            filter %>% dplyr::filter(
                Personalized.Medicine == filter_personalized_medicine |
                    Personalized.Medicine == "Control"
            )
        if (!is.na(filter_expandability))
            filter <-
            filter %>% dplyr::filter(Expandability == filter_expandability |
                                         Expandability == "Control")
        if (!is.na(filter_culture_duration))
            filter <-
            filter %>% dplyr::filter(Culture.Duration == filter_culture_duration |
                                         Culture.Duration == "Control")
        if (!is.na(filter_cell_source))
            filter <-
            filter %>% dplyr::filter(Cell.Source == filter_cell_source |
                                         Cell.Source == "Control")
        if (!is.na(filter_viral_transduction))
            filter <-
            filter %>% dplyr::filter(
                Viral.Transduction == filter_viral_transduction |
                    Viral.Transduction == "Control"
            )
        
        return(filter)
    })
    
    train_Data_filtered <- reactive({
        # Countdata columns need to match coldata Samples....
        # Get filtered coldata
        train_colD <- train_colData_filtered()
        
        select_cols <- c("ensembl_gene_name", train_colD$id)
        
        # Only return ensembl_gene_name + filtered columns
        return(train_Data[select_cols])
    })
    
    # User input: coldata; return NULL if not uploaded
    user_colData <- reactive({
        file <- input$file_ext_coldata
        if (is.null(file)) {
            return()
        }
        read.table(
            file = file$datapath,
            sep = input$col_sep,
            header = T
        )
    })
    
    # User input: countdata; return NULL if not uploaded
    user_Data <- reactive({
        file <- input$file_ext_countdata
        if (is.null(file)) {
            return()
        }
        cont <-
            read.table(
                file = file$datapath,
                sep = input$count_sep,
                header = T
            )
        # Update first column name to 'ensembl_gene_name'
        colnames(cont)[1] <- "ensembl_gene_name"
        return(cont)
    })
    
    override_count <- reactive({
        input$override
    })
    
    parse_genes <- function(input, nl2br = F) {
        if (nl2br == T) {
            # Paste list to linebreak list
            list <- paste(input, collapse="\n")
            return(list)
        } 
        else {
            # Split inputted list by newlines
            list <- strsplit(input, "\n")
            # Create vector element --> vector now in [[1]]
            list <- as.vector(list)
            # Remove all empty lines and NA's
            list <- stri_remove_empty(list[[1]], na_empty = TRUE)
            return(list)
        }
    }
    
    geneset_genes_filtered <- function() {
        genes <- genes_in_merged_data()
        submitted_genes <-
            parse_genes(input = input$geneset_genes)
        not_in_dataset <- setdiff(submitted_genes, genes)
        # If no genes not in dataset, return all submitted genes
        if (length(not_in_dataset) == 0) {
            return(submitted_genes)
        } else {
            # Don't return anything
            return()
        }
    }
    
    geneset_sample_genes <- function() {
        return(
            "AGT
AGXT
AHSG
ALB
ALDOB
AMBP
APCS
APOA1
APOA2
APOC1
APOC2
APOC3
APOE
APOH
ATF5
C1R
C3
CRP
CYP2E1
FGA
FGB
FGG
FGL1
FTL
GC
HP
HPD
HPX
IFITM3
IGFBP4
ITIH4
MT1G
MT2A
MTATP6P1
ORM1
ORM2
PEBP1
RBP4
SAA1
SAA2
SAA4
SERPINA1
SERPINA3
SERPINC1
SERPINF2
SERPING1
TF
TTR
VTN"
        )
    }
    
    load_geneset <- function(name, unique = FALSE) {
        if (name == "GO Biological Process") {
            geneset_file <- read.csv("geneset_go_biologicalprocess.csv", header = TRUE)
            colnames(geneset_file) <- c("Gene", "ID", "Description")
        }
        else if (name == "KEGG") {
            geneset_file <- read.csv("geneset_KEGG_Pathways.csv", header = TRUE)
            colnames(geneset_file) <- c("ID", "Gene", "Description")
        }
        else if (name == "DisGeNET") {
            geneset_file <- read.csv("geneset_DisGeNET.csv", header = TRUE)
            colnames(geneset_file) <- c("Gene", "ID", "Description")
        } 
        else {
            return()
        }
        
        if (isTRUE(unique)) {
            return(unique(geneset_file$Description))
        } 
        else {
            return(geneset_file)
        }
        
    }
    
    update_geneset_search <- function() {
        # Based on input sourge (e.g. KEGG), load descriptions and update input search list
        geneset_list_unique <- load_geneset(input$geneset_source, TRUE)
        geneset_list_unique <- c("", geneset_list_unique)
        updateSelectizeInput(inputId = "geneset_search",
                             choices = geneset_list_unique)
    }
    
    update_geneset_from_source <- function(source = F, descr = F) {
        # Check if variables are set
        if (source == F) return()
        if (descr == F | descr == "") return()
        
        # Load geneset from input source (e.g. KEGG, GO ontology...)
        geneset <- load_geneset(source) %>% filter(Description == descr)
        # Select genes, then parse the list to a character with linebreaks
        parsed_genes <- parse_genes(input = intersect(geneset$Gene, genes_in_merged_data()), nl2br = TRUE)
        # Update textarea..
        updateTextAreaInput(inputId = "geneset_genes",
                            value = parsed_genes)
    }
    
    observeEvent(input$geneset_source, {
        update_geneset_search()
    })
    
    observeEvent(input$geneset_search, {
        
        update_geneset_from_source(source = input$geneset_source, descr = input$geneset_search)
    })
    
    observeEvent(input$geneset_reset, {
        update_geneset_search()
        updateTextAreaInput(inputId = "geneset_genes",
                            value = geneset_sample_genes())
    })
    
    # Get gene annotations (names) from biomart - or local storage if the file exists
    get_annotations <- function() {
        # For now, assume we've downloaded and saved the annotations - otherwise, yse code below to download the annotations
        annot <- read.table("gene_annotations.csv", sep = ",")
        #if (is.null(annot)) {
        #    ## Build a biomart query
        #    # In the example below, I use the human gene annotation from Ensembl release 82 located on "sep2015.archive.ensembl.org", More about the ensembl_build can be found on "http://www.ensembl.org/info/website/archives/index.html"
        #    dataset = "hsapiens_gene_ensembl"
        #    mart = biomaRt::useMart(
        #        "ENSEMBL_MART_ENSEMBL",
        #        dataset = dataset,
        #        host = paste0("jan2020", ".archive.ensembl.org"),
        #        path = "/biomart/martservice",
        #        archive = FALSE
        #    )
        #    #listAttributes(mart)
        #
        #    # Check if annotations were loaded from biomart before - if so, don't load them.
        #    annot <-
        #        biomaRt::getBM(
        #            mart = mart,
        #            attributes = c("ensembl_gene_id", "external_gene_name"),
        #            useCache = FALSE
        #        ) #useCahche=FALSE -> need to update version of R / bioconductor
        #    annot <-
        #        dplyr::transmute(annot, ensembl_gene_id, external_gene_name)
        #    write.table(annot, "gene_annotations.csv", sep = ",")
        #}
        return(annot)
    }
    
    # Manual colour scale for all plots
    get_color_scale <-
        function(merged_colD = merge_colData(),
                 train_colD = train_colData_filtered()) {
            scale <- c(
                "CBD" = "#2CA02C",
                "Chol_HLC" = "#95D095",
                "Fetal_Hep" = "#FF7F00",
                "FHep_HLC" = "#FFBF7F",
                "Fib_HLC" = "#8FBBDA",
                "Fibroblast" = "#1F77B4",
                "Hep_HLC" = "#FBB4AE",
                "HepG2" = "#E7298A",
                "Liver" = "#B59516",
                "PHH" = "#D62728",
                "PSC" = "#9467BD",
                "PSC_HLC" = "#CAB3DE"
            )
            
            extra <- c(
                "#00BFC4",
                # Extra colors for user-uploaded data
                "#F394C5",
                "#FFD600",
                "#B8B8B8",
                "#271C6E",
                "#165E42",
                "#5E2019",
                "#4D4D4D",
                "#6E1C5E",
                "#6E601C",
                "white",
                # And more white colors in case the user uploads way too many types
                "white",
                "white",
                "white",
                "white",
                "white",
                "white",
                "white",
                "white",
                "white",
                "white",
                "white",
                "white",
                "white"
            )
            
            # Select unique types from merged and train lists
            m_colD <- unique(merged_colD$Type)
            t_colD <- unique(train_colD$Type)
            
            # Check for types not in train coldata
            unique_types <- setdiff(m_colD, t_colD)
            
            # Count types not in train coldata
            n = length(unique_types)
            
            if (n == 0) {
                # Return scale if it exactly fits the standard coldata
                return(scale)
            } else{
                # If n is >0 then return extra colors, with the correct 'names'
                extra <- extra[1:(n)]
                names(extra) <- unique_types
                return(c(scale, extra))
            }
            
        }
    
    ## Function user_file_errors()
    # Test for errors in the uploaded coldata and countdata files
    user_file_errors <-
        function(coldata = user_colData(),
                 countdata = user_Data(),
                 train_D = train_Data,
                 merged_colD = merge_colData()) {
            # Do not return error if only using training data
            if (input$use_training_data_only == TRUE)
                return(FALSE)
            
            ## Check upload status
            # Check if countdata and coldata files were uploaded. If not, throw error. However, user won't usually get to this status screen if not uploaded...
            if (is.null(coldata))
                return("Please upload valid coldata file. [e1]")
            if (is.null(countdata))
                return("Please upload valid countdata file. [e2]")
            
            ## Data checks of uploaded files
            # Check if coldata contains necessary columns
            if (ncol(coldata) != 4)
                return(
                    "The uploaded coldata does not contain 4 columns. \n\nPlease check and fix the file. Maybe you chose the incorrect separator? [e3]"
                )
            
            # Check if colnames of coldata are correct
            if ("id" %ni% colnames(coldata))
                return(
                    "The uploaded coldata file does not contain column: `id`.\n\nPlease check and fix the file. Maybe you chose the incorrect separator? [e4]"
                )
            if ("Author" %ni% colnames(coldata))
                return(
                    "The uploaded coldata file does not contain column: `Author`.\n\nPlease check and fix the file. Capital letters matter. Maybe you chose the incorrect separator? [e5]"
                )
            if ("Type" %ni% colnames(coldata))
                return(
                    "The uploaded coldata file does not contain column: `Type`.\n\nPlease check and fix the file. Capital letters matter. Maybe you chose the incorrect separator? [e6]"
                )
            if ("Type2" %ni% colnames(coldata))
                return(
                    "The uploaded coldata file does not contain column: `Type2`.\n\nPlease check and fix the file. Capital letters matter. Maybe you chose the incorrect separator? [e6]"
                )
            
            # Check if order of columns in coldata is correct
            
            
            # Check if all sample names in user colData are unique and not duplicating training coldata names
            if (anyDuplicated(merged_colD$id) != 0)
                return(
                    "The uploaded sample ids are not unique. Please check the sample ids. \n\nBeware they may not overlap with sample names in the trainig dataset. [e14]"
                )
            
            # Check if countdata is in correct format
            
            # Check if countdata first column contains 'ensembl_gene_name', aka 16-18 characters
            if (nchar(countdata[1, 1]) < 16 |
                nchar(countdata[1, 1]) > 18)
                return(
                    paste(
                        "It seems like the first column of countdata does not contain the ensembl gene id.\n\nPlease check and fix the file. [ nchar:",
                        nchar(countdata[1, 1]),
                        "] [e8]"
                    )
                )
            
            # Check if countdata first column is character, not numeric
            if (!is.character(countdata[, 1]))
                return(
                    "It seems like the first column of countdata does not contain the ensembl gene id.\n\nPlease check and fix the file. [e9]"
                )
            
            # Check if coldata row number matches countdata col number (minus 1 because of column ensembl_gene_name)
            if (nrow(coldata) != ncol(countdata) - 1)
                return(
                    paste(
                        "The uploaded coldata and countdata files do not contain the same number of samples [e10]. \nColdata:\t",
                        nrow(coldata),
                        "\nCountdata:\t",
                        ncol(countdata) - 1,
                        "\n\nPlease check and fix your files. \n\nMaybe you chose the incorrect separator?"
                    )
                )
            
            # Check if countdata has the same number of genes as the training count data
            if (nrow(countdata) < nrow(train_D) &
                isFALSE(override_count()))
                return(
                    paste(
                        "The uploaded countdata file contains fewer genes (",
                        nrow(countdata),
                        ") than our training countdata (",
                        nrow(train_D),
                        ").\n\nPlease check and fix your files. \n\nAlternatively, you may choose to continue, by activating the 'count override' toggle on the left. [e12]"
                    )
                )
            if (nrow(countdata) > nrow(train_D) &
                isFALSE(override_count()))
                return(
                    paste(
                        "The uploaded countdata file contains more genes (",
                        nrow(countdata),
                        ") than our training countdata (",
                        nrow(train_D),
                        ").\n\nPlease check and fix your files. \n\nAlternatively, you may choose to continue, by activating the 'count override' toggle on the left. [e13]"
                    )
                )
            
            # Check if all samplenames (column `id` in user_colData) exist as column names in countdata (user_Data)
            by(coldata, seq_len(nrow(coldata)), function(row) {
                if (row$id %ni% colnames(countdata)) {
                    return(
                        paste(
                            "The following sample in your coldata does not exist in your countdata:",
                            row$id,
                            "[e11.n]"
                        )
                    )
                }
            })
            
            # No errors
            return(FALSE)
        }
    
    user_file_incompatible <- function() {
        # If only using training data, return FALSE
        if (input$use_training_data_only == TRUE)
            return(FALSE)
        
        # Next, load user coldata
        user_colD <- user_colData()
        if ("PHH" %ni% user_colD$Type &
            "Liver" %ni% user_colD$Type)
            return("Uploaded data does not contain any of Type 'PHH' or 'Liver'. [e20]")
        
        # Next, test compatibility with our PHH/Liver
        
        load("WarningInput.RData")
        # Assign qn_vst_data to another variable/object and transpose it
        transposed_vst_data <- as.data.frame(t(qn_vst_data()))
        # Subset the predictor gene set
        transposed_vst_data <- transposed_vst_data[, predictor]
        # Run the Random Forest classifier
        class_prob <-
            predict(model, transposed_vst_data, type = "prob")
        class_prob <- as.data.frame(class_prob)
        
        not_matching <-
            class_prob %>% dplyr::filter(`PHH/Liver` < 0.5)
        
        if (nrow(not_matching) > 0) {
            # Get coldata, filter only PHH and liver samples
            ColD_PHHLiver <-
                merge_colData() %>% filter(Type %in% c("PHH", "Liver"))
            
            # Get rownames = sample names of samples with PHH/Liver < 0.5
            warning_samples <- rownames(not_matching)
            
            # Intersect, to only show PHH/Liver samples
            warning_list <-
                intersect(warning_samples, ColD_PHHLiver$id)
            #warning_list <- warning_samples
            #warning_list <- ColD_PHHLiver$id
            
            samples <- paste(shQuote(warning_list), collapse = ", ")
            message <-
                paste(
                    "Your PHH or Liver samples do not cluster with the training PHH or Liver samples. The following samples are strange (training samples may have been skewed by uploaded samples):\n ",
                    samples,
                    ". You may proceed at your own risk. [e21]"
                )
            return(message)
        }
        
        return(NULL)
        
    }
    
    ## Function merge_colData
    # To merge the training and user coldatasets
    merge_colData <-
        function(user_colD = user_colData(),
                 train_colD = train_colData_filtered()) {
            # Add 'source' = training
            train_colD$Source <- "Training"
            
            # If only training data is to be used, don't merge, but return training data
            if (input$use_training_data_only == TRUE)
                return(train_colD)
            
            # Otherwise: if user coldata is not set, then fail
            if (is.null(user_colD))
                return()
            
            # Add "filter" columns - otherwise cannot merge (because not same columns in both datasets)
            user_colD$Culture.Type <-
                user_colD$Personalized.Medicine <-
                user_colD$Expandability <-
                user_colD$Culture.Duration <-
                user_colD$Cell.Source <-
                user_colD$Viral.Transduction <-
                NA
            
            # Add 'source' = user
            user_colD$Source <- "User"
            
            # Merge the datasets
            merge_col <- rbind(train_colD, user_colD)
            #merge_col <- train_colD
            
            return(merge_col)
        }
    
    ## Function merge_Data
    # To merge the training and user count datasets
    merge_Data <-
        function(user_D = user_Data(),
                 train_D = train_Data_filtered()) {
            # If only training data is to be used, don't merge, but return training data
            if (is.null(user_D) &
                input$use_training_data_only == FALSE)
                return()
            
            # Merge the datasets (only if use_training_data_only is FALSE)
            if (input$use_training_data_only == FALSE) {
                merge_D <- merge(train_D, user_D, by = "ensembl_gene_name")
            } else{
                merge_D <- train_D
            }
            
            # Move 'ensembl_gene_name' to first column - otherwise, everything crashes
            col_idx <- grep("ensembl_gene_name", names(merge_D))
            merge_D <-
                merge_D[, c(col_idx, (1:ncol(merge_D))[-col_idx])]
            
            # (1) Removing the string after "." in the ensembl gene names
            removedot <-
                gsub("\\..*", "", merge_D[, 1]) #the gene name should be in the first column, otherwise adjust the value of "[,1]"
            merge_D[, 1] <- removedot
            
            # (2) Removing duplicate of ensembl genes
            merge_D <-
                merge_D[!duplicated(merge_D$ensembl_gene_name), ] # "ensembl_gene_name" should be the header for where the ensembl gene column, otherwise adjust it accordingly
            
            # (3) Removing genes with total row counts <10
            merge_D <- merge_D[rowSums(merge_D[, -1]) >= 10, ]
            
            return(merge_D)
        }
    
    vst_data <-
        function(data = merge_Data(),
                 colData = merge_colData()) {
            dds <- DESeqDataSetFromMatrix(
                countData = data,
                colData = colData,
                design =  ~ Type,
                #can be any random column name from the colData, does not matter
                tidy = T
            )
            
            #applying the VST transformation
            vst_data <- vst(dds, blind = FALSE)
            
            #change the vst data format from a "Large DESeqTransform" into a matrix
            vst_data <-
                as.matrix(assay(vst_data)) #this step can be combined with the previous step also if you like
            
            return(vst_data)
        }
    
    qn_vstdata_stored <- 0
    qn_vstdata_stored <- reactiveVal(NULL)
    
    qn_vst_data <- function(data = NULL) {
        # If vstdata exists, don't run deseq again... Waste of resources. Just load from variable.
        if (!is.null(qn_vstdata_stored())) {
            return(qn_vstdata_stored())
        }
        
        if (is.null(data))
            data <- vst_data()
        
        qn_vstdata <- normalize.quantiles(data)
        
        #after quantile normalization, the colnames and rownames are gone, so we need to retrieve them back
        rownames(qn_vstdata) <- rownames(data)
        colnames(qn_vstdata) <- colnames(data)
        
        qn_vstdata_stored(qn_vstdata)
        
        return(qn_vstdata)
    }
    
    qn_vstdata2_stored <- 0
    qn_vstdata2_stored <- reactiveVal(NULL)
    
    # Function to add gene names to qn_vst_data
    qn_vst_data_2 <-
        function(qn_vstdata = qn_vst_data(),
                 annot = get_annotations()) {
            # If vstdata2 exists, don't run this function again... Waste of resources. Just load from variable.
            if (!is.null(qn_vstdata2_stored())) {
                return(qn_vstdata2_stored())
            }
            
            # Filter and re-order gene.annotations to match the order in your input genes list
            gene_name <-
                annot %>% dplyr::filter(ensembl_gene_id %in% rownames(qn_vstdata))
            gene_name <-
                gene_name[order(match(gene_name$ensembl_gene_id, rownames(qn_vstdata))), ]
            
            # Merging annotation with dataset
            qn_vstdata2 <-
                merge(
                    as.data.frame(gene_name),
                    as.data.frame(qn_vstdata),
                    by.x = "ensembl_gene_id",
                    by.y = 0,
                    all = TRUE
                )
            
            # Tidy up the final dataframe
            qn_vstdata2 <-
                qn_vstdata2[!duplicated(qn_vstdata2$external_gene_name), ] #remove duplicate of gene name
            rownames(qn_vstdata2) <-
                qn_vstdata2[, 2] #make the HGNC as rownames
            qn_vstdata2 <-
                qn_vstdata2[, c(-1, -2)] #remove unnecessary column
            
            qn_vstdata2_stored(qn_vstdata2)
            
            return(qn_vstdata2)
        }
    
    # Function to check all genes in the merged dataset -  genes in this list may not correspond to all genes in biomart geneset
    genes_in_merged_data <-
        function(merged_D = merge_Data(),
                 annot = get_annotations()) {
            # Set rownames of merged_D (Ensembl ids)
            rownames(merged_D) <- merged_D[, 1]
            
            # Filter and re-order gene.annotations to match the order in your input genes list
            gene_name <-
                annot %>% dplyr::filter(ensembl_gene_id %in% rownames(merged_D))
            gene_name <-
                gene_name[order(match(gene_name$ensembl_gene_id, rownames(merged_D))), ]
            
            # Merging annotation with dataset
            gene_list <-
                merge(
                    as.data.frame(gene_name),
                    as.data.frame(merged_D),
                    by.x = "ensembl_gene_id",
                    by.y = 0,
                    all = TRUE
                )
            
            # Tidy up the final dataframe
            gene_list <-
                gene_list[!duplicated(gene_list$external_gene_name), ] #remove duplicate of gene name
            rownames(gene_list) <-
                gene_list[, 2] #make the HGNC as rownames
            
            gene_list <- rownames(gene_list)
            
            return(gene_list)
        }
    
    # Function to plot PCA plot
    plot_pca <-
        function(merged_colData = merge_colData(),
                 qn_vstdata = qn_vst_data()) {
            if (exists("var_plot_pca"))
                rm(var_plot_pca, pos = ".GlobalEnv")
            
            #specify the number genes to take into account in the pca plot
            # ntop
            
            #specify the labeling, we just take all the colname from merged_Coldata
            intgroup = colnames(merged_colData)
            
            #calculate row variance and ordered them decending to select the highest variance genes (ntop)
            rv <- rowVars(as.matrix(qn_vstdata[, ]))
            select <-
                order(rv, decreasing = TRUE)[seq_len(min(input$pca_ntop, length(rv)))]
            
            #create the PCA
            pca_res <-
                prcomp(t(qn_vstdata[select, ]), scale. = FALSE)
            percentVar <- pca_res$sdev ^ 2 / sum(pca_res$sdev ^ 2)
            
            #prepare the label for the PCA plot
            intgroup.df <-
                as.data.frame(merged_colData[, intgroup, drop = FALSE])
            
            #combine the PCA with the label
            d <- data.frame(
                PC1 = pca_res$x[, 1],
                PC2 = pca_res$x[, 2],
                intgroup.df,
                name = colnames(qn_vstdata[, ])
            )
            
            #plot the PCA
            plot <-
                ggplot(data = d, aes_string(
                    x = "PC1",
                    y = "PC2",
                    colour = "Type"
                )) +
                geom_point(aes(shape = Source), size = 3) +
                xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
                ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
                coord_fixed() +
                scale_colour_manual(values = get_color_scale())
            
            return(plot)
        }
    
    get_matrix <-
        function(plot = "heatmap",
                 qn_vstdata2 = qn_vst_data_2()) {
            #upload the gene to be selected
            # selectgene <- as.vector(read.csv("example_geneset.csv", header=TRUE)[,])
            
            if (plot == "heatmap") {
                selectgene <- geneset_genes_filtered()
                if (length(selectgene) < 2)
                    return()
            } else if (plot == "boxplot") {
                selectgene <- geneset_genes_filtered()
                if (length(selectgene) < 2)
                    return()
            } else{
                selectgene <- parse_genes(geneset_sample_genes())
            }
            
            #subset the dataset using the selected gene set
            matrix <- qn_vstdata2[selectgene, ]
            
            #remove NA containing rows (CAN WE GIVE WARNING FOR GENES THAT ARE NOT AVAILABLE?)
            matrix <- matrix[complete.cases(matrix), ]
            
            #apply mean centering for easier data interpretation
            matrix <- matrix - rowMeans(matrix)
            
            return(matrix)
        }
    
    # Function to plot heatmap plot
    plot_heatmap <-
        function(colD = merge_colData(),
                 matrix = get_matrix()) {
            if (exists("var_plot_heatmap"))
                rm(var_plot_heatmap, pos = ".GlobalEnv")
            
            #set rownames of coldata - needed for pheatmap (matrix colnames need to be coldata rownames)
            rownames(colD) <- colD[, 1]
            
            #define the sample information displayed for the heatmap legend
            coldata <- colD[, c(0, 2, 3, 11)]
            
            #set the upper and lower border of the heatmap intensity (CAN WE MAKE THIS ADJUSTABLE BY SLIDER?)
            breaksList = seq(
                #isolate(
                input$heatmap_intensity[1]
                #    )
                ,
                #isolate(
                input$heatmap_intensity[2]
                #    )
                ,
                by = 0.01
            )
            
            # Check if rownames are correct - but this should be, otherwise the heatmap won't show anyways...
            if (length(setdiff(colnames(matrix), rownames(coldata))) != 0)
                cat(
                    "Colnames matrix not rownames coldata [e15]. \nDifferences: ",
                    setdiff(colnames(matrix), rownames(coldata))
                )
            
            cluster_rows_TF <-
                ifelse(
                    !is.null(
                        #isolate(
                        input$heatmap_cluster_rows
                        #    )
                    ) &
                        #isolate(
                        input$heatmap_cluster_rows
                    #    ) 
                    == TRUE,
                    TRUE,
                    FALSE
                )
            cluster_cols_TF <-
                ifelse(
                    !is.null(
                        #isolate(
                        input$heatmap_cluster_cols
                        #    )
                    ) &
                        #isolate(
                        input$heatmap_cluster_cols
                    #    ) 
                    == TRUE,
                    TRUE,
                    FALSE
                )
            
            #create the heatmap (CAN WE MAKE THE CLUSTER OPTION TWEAKABLE BY USER?)
            plot <- pheatmap(
                matrix,
                color = colorRampPalette(c(
                    "dodgerblue3", "white", "firebrick1"
                ))(length(breaksList)),
                fontsize = 7,
                fontsize_row = 6,
                fontsize_col = 6,
                cellheight = 7,
                cellwidth = NA,
                breaks = breaksList,
                annotation_col = coldata,
                cluster_rows = cluster_rows_TF,
                cluster_cols = cluster_cols_TF,
                #cluster_rows = TRUE,
                #cluster_cols = TRUE,
                cutree_rows = 1,
                cutree_cols = 3,
                border_color = NA,
                annotation_colors = list(
                    Source = c("Training" = "snow3",
                               "User" = "red"),
                    Type = get_color_scale()
                ),
                drop_levels = TRUE
            )
            
            return(plot)
        }
    
    # Function to plot boxplot
    plot_boxplot <-
        function(matrix = get_matrix(plot = "boxplot", qn_vstdata2 = qn_vst_data_2()),
                 colD = merge_colData()) {
            if (exists("var_plot_boxplot"))
                rm(var_plot_boxplot, pos = ".GlobalEnv")
            
            #set rownames of coldata - needed
            rownames(colD) <- colD[, 1]
            
            #calculate euclidean distance between sample in a certain gene set
            euclid <-
                as.matrix(dist(t(matrix), method = "euclidean"))
            
            #subset only the distance to PHH and liver controls
            comdist <- euclid[grep("PHH|Liver", rownames(euclid)), ]
            
            #transform the euclidean distance into distance-based similarity
            comsim <- (max(comdist) - comdist) / max(comdist)
            
            #remove the rownames, it will mess up downstream table transformation
            rownames(comsim) <- NULL
            
            #we change the table format from "wide" into "long"
            box_dist <-
                melt(comsim, varnames = c("Number", "Names"))
            box_dist <-
                box_dist[, -1] #we don't need the first column
            
            #adding author and cell type information
            com_box_dist <-
                merge(box_dist,
                      colD,
                      by.x = "Names",
                      by.y = 0,
                      all = TRUE)
            
            #combining author and cell type so replicates will be plotted in the same box
            com_box_dist <- com_box_dist %>%
                unite(
                    com,
                    Author,
                    Type,
                    Type2,
                    sep = "_",
                    na.rm = T,
                    remove = F
                )
            
            plot <-
                ggplot(com_box_dist,
                       aes(
                           x = fct_reorder(com, value),
                           y = value,
                           fill = Type,
                           colour = factor(Source)
                       )) +
                geom_boxplot() +
                ylim(0, 1) +
                xlab("") +
                ylab(expression("Distance-based similarity (to PHH/Liver)")) +
                theme(axis.text.x = element_text(
                    angle = 270,
                    vjust = 0.5,
                    hjust = 0
                )) +
                scale_colour_manual(values = c("User" = "red",
                                               "Training" = "black"),
                                    name = "Data source") +
                scale_fill_manual(values = get_color_scale())
            return(plot)
        }
    
    # Status tab - check input
    output$status <- renderTable({
        errors <- user_file_errors()
        if (errors != F) {
            out <- as.data.frame(errors)
            colnames(out) <- c("Not OK!")
        } else{
            if (input$use_training_data_only == TRUE)
                out <-
                    as.data.frame(c("Using training data only.\n\nPlease proceed."))
            else
                out <-
                    as.data.frame(c("Your files passed all checks. \n\nPlease proceed."))
            colnames(out) <- c("OK!")
        }
        return(out)
    })
    
    output$table_training <- renderTable({
        if (is.null(train_colData_filtered()))
            return()
        train_colData_filtered()
    })
    
    output$table_user <- renderTable({
        if (user_file_errors() != FALSE)
            return("There are errors. Please check the 'status' tab. [e16]")
        if (input$use_training_data_only == TRUE)
            return("N/A")
        user_colData()
    }) 
    
    files_updated <- reactive({
        list(
            input$file_ext_coldata,
            input$file_ext_countdata,
            input$filter_culture_type,
            input$filter_personalized_medicine,
            input$filter_expandability,
            input$filter_culture_duration,
            input$filter_cell_source,
            input$filter_viral_transduction
        )
    })
    
    observeEvent(files_updated(), {
        # Stop if only using training data
        # (not any more now we can filter data)
        #if (input$use_training_data_only == TRUE)
        #    return()
        # New user coldata uploaded - if files were uploaded AND no errors with the files AND qn_vstdata exists - remove qn_vstdata
        #if (!is.null(input$file_ext_coldata) & !is.null(input$file_ext_countdata) & user_file_errors() == FALSE) {
        if (!is.null(qn_vstdata_stored()))
            qn_vstdata_stored(NULL)
        if (!is.null(qn_vstdata2_stored()))
            qn_vstdata2_stored(NULL)
        #}
        
        data_warnings <- function() {
            # Don't throw a warning if only trainig data is used
            if (input$use_training_data_only == TRUE)
                return()
            
            # Throw warning if something is wrong tho
            
            warnings <- user_file_incompatible()
            
            if (!is.null(warnings)) {
                warning <- as.data.frame(
                    paste(
                        "Your dataset may be incompatible with HLCompR. For more info, please check our publication. ",
                        warnings
                    )
                )
                colnames(warning) <- c("Warning")
                
                return(warning)
            }
            return()
        }
        
        output$data_warnings_pca <- renderTable({
            return(data_warnings())
        })
        output$data_warnings_heatmap <- renderTable({
            return(data_warnings())
        })
        output$data_warnings_boxplot <- renderTable({
            return(data_warnings())
        })
        
    })
    
    # Draw PCA plot
    output$PCA <- renderPlot({
        if (user_file_errors() != F) {
            return()
        }
        var_plot_pca <<-
            plot_pca()
        return(var_plot_pca)
    })
    
    # Draw heatmap
    output$heatmap <- NULL
    
    check_geneset_errors <- function() {
        # Display error if genes that were input don't exist in get_annotations()
        #annot <- get_annotations()
        genes <- genes_in_merged_data()
        submitted_genes <-
            parse_genes(input = input$geneset_genes)
        errors <- setdiff(submitted_genes, genes)
        out <- NULL
        if (length(errors) != 0)
            out <- as.data.frame(c(
                paste(
                    "The following genes do not exist:",
                    toString(errors)
                )
            ))
        
        if (length(submitted_genes) < 2) {
            out <- as.data.frame(c(paste(
                "Please input at least two genes."
            )))
        }
        
        if (!is.null(out))
            colnames(out) <- c("Error")
        else
            return()
        return(out)
    }
    
    output$geneset_errors <- renderTable({
        return(check_geneset_errors())
    }) 
    
    output$heatmap_errors <- renderTable({
        return(check_geneset_errors())
    }) 
    
    output$boxplot_errors <- renderTable({
        return(check_geneset_errors())
    }) 
    
    output$heatmap <- renderPlot({
        
        if (user_file_errors() != F)
            return()
        
        genes <- genes_in_merged_data()
        submitted_genes <-
            parse_genes(input = input$geneset_genes)
        errors <- setdiff(submitted_genes, genes)
        if (length(errors) != 0)
            return()
        
        var_plot_heatmap <<-
            plot_heatmap()
        return(var_plot_heatmap)
    })
    
    # Draw boxplot
    output$boxplot <- NULL
    
    output$boxplot <- renderPlot({
        if (user_file_errors() != F)
            return()
        
        genes <- genes_in_merged_data()
        submitted_genes <-
            parse_genes(input = input$geneset_genes)
        errors <- setdiff(submitted_genes, genes)
        if (length(errors) != 0)
            return()
        
        var_plot_boxplot <<-
            plot_boxplot()
        return(var_plot_boxplot)
    })
    
    
    output$download_PCA <- downloadHandler(
        filename = function() {
            paste0("pca-", Sys.Date(), ".pdf")
        },
        content = function(file) {
            if (!exists("var_plot_pca"))
                return()
            pdf(file)
            print(var_plot_pca)
            dev.off()
        }
    )
    
    output$download_heatmap <- downloadHandler(
        filename = function() {
            paste0("heatmap-", Sys.Date(), ".pdf")
        },
        content = function(file) {
            #showModal(modalDialog("Please wait a moment while your download is being prepared.", footer = NULL))
            #on.exit(removeModal())
            if (!exists("var_plot_heatmap"))
                return()
            pdf(file)
            print(var_plot_heatmap)
            #print(plot_heatmap())
            dev.off()
        }
    )
    
    output$download_boxplot <- downloadHandler(
        filename = function() {
            paste0("boxplot-", Sys.Date(), ".pdf")
        },
        content = function(file) {
            #showModal(modalDialog("Please wait a moment while your download is being prepared.", footer = NULL))
            #on.exit(removeModal())
            if (!exists("var_plot_boxplot"))
                return()
            pdf(file)
            print(var_plot_boxplot)
            #print(plot_boxplot())
            dev.off()
        }
    )
    
    output$user_file_input <- renderUI({
        if (input$use_training_data_only == FALSE) {
            fluidRow(
                column(
                    12,
                    fileInput(
                        "file_ext_coldata",
                        "Upload coldata (.csv)",
                        multiple = F,
                        accept = c(
                            "text/csv",
                            "text/comma-separated-values,text/plain",
                            ".csv"
                        )
                    ),
                    fileInput(
                        "file_ext_countdata",
                        "Upload countdata (.csv)",
                        multiple = F,
                        accept = c(
                            "text/csv",
                            "text/comma-separated-values,text/plain",
                            ".csv"
                        )
                    ),
                    
                    hr(),
                    fluidRow(
                        column(
                            6,
                            radioButtons(
                                inputId = 'col_sep',
                                label = 'Coldata separator',
                                choices = c(
                                    Comma = ',',
                                    Semicolon = ';',
                                    Tab = '\t',
                                    Space = ''
                                ),
                                selected = ','
                            )
                        ),
                        column(
                            6,
                            radioButtons(
                                inputId = 'count_sep',
                                label = 'Countdata separator',
                                choices = c(
                                    Comma = ',',
                                    Semicolon = ';',
                                    Tab = '\t',
                                    Space = ''
                                ),
                                selected = ','
                            )
                        )
                    ),
                    hr(),
                    helpText(
                        "The number of genes in your countdata file may differ from the number of genes in our training countdata file. If this occurs, an error will be shown upon upload. You may override this error using this toggle. However, we do not recommend this, as unexpected results may occur."
                    ),
                    prettySwitch(
                        inputId = "override",
                        label = "Override",
                        value = FALSE
                    )
                )
            )
        }
    })
    
    # Tabs
    output$tabs <- renderUI({
        if ((is.null(user_colData()) |
             is.null(user_Data())) &
            input$use_training_data_only == FALSE) {
        }
        else {
            tabsetPanel(
                tabPanel("1. Status", tableOutput("status")),
                tabPanel(
                    "2. Filters",
                    fluidRow(
                        #column(12,
                        #       h5("Filters")),
                        column(
                            4,
                            h6("Culture type"),
                            awesomeRadio(
                                inputId = "filter_culture_type",
                                choices = c("No preference",
                                            "2D",
                                            "3D"),
                                label = "",
                                selected = "No preference"
                            )
                        ),
                        column(
                            4,
                            h6("Personalized medicine"),
                            awesomeRadio(
                                inputId = "filter_personalized_medicine",
                                choices = c("No preference",
                                            "Yes",
                                            "No"),
                                label = "",
                                selected = "No preference"
                            )
                        ),
                        column(
                            4,
                            h6("Expandability"),
                            awesomeRadio(
                                inputId = "filter_expandability",
                                choices = c("No preference",
                                            "Yes",
                                            "No"),
                                label = "",
                                selected = "No preference"
                            )
                        ),
                        column(
                            4,
                            h6("Culture duration"),
                            awesomeRadio(
                                inputId = "filter_culture_duration",
                                choices = c("No preference",
                                            "Short-term",
                                            "Long-term"),
                                label = "",
                                selected = "No preference"
                            )
                        ),
                        column(
                            4,
                            h6("Cell source"),
                            awesomeRadio(
                                inputId = "filter_cell_source",
                                choices = c(
                                    "No preference",
                                    "Fibroblast",
                                    "Hepatocyte",
                                    "Cholangiocyte",
                                    "PSC",
                                    "Fetal hepatocyte"
                                ),
                                label = "",
                                selected = "No preference"
                            )
                        ),
                        column(
                            4,
                            h6("Viral transduction"),
                            awesomeRadio(
                                inputId = "filter_viral_transduction",
                                choices = c("No preference",
                                            "Yes",
                                            "No"),
                                label = "",
                                selected = "No preference"
                            )
                        )
                    )
                ),
                tabPanel("3. Samples", fluidRow(
                    column(
                        12,
                        h6("Uploaded coldata"),
                        tableOutput("table_user") %>% withSpinner(color =
                                                                      "lightgrey")
                    ),
                    column(
                        12,
                        h6("Training coldata"),
                        tableOutput("table_training") %>% withSpinner(color =
                                                                          "lightgrey")
                    )
                    
                )),
                tabPanel("4. PCA plot",
                         fluidRow(
                             column(
                                 3,
                                 sliderInput(
                                     "pca_ntop",
                                     "Number of genes",
                                     min = 100,
                                     max = 40000,
                                     value = 5000
                                 ),
                             ),
                             column(
                                 9,
                                 tableOutput("data_warnings_pca") %>% withSpinner(color = "lightgrey", proxy.height = "100px"),
                                 plotOutput("PCA", height = "600px") %>% withSpinner(color = "lightgrey")
                             ),
                             column(12,
                                    downloadButton("download_PCA", "Download"))
                         )),
                tabPanel("5. Geneset",
                         fluidRow(column(4,
                                         h4("1"),
                                         h6("Select source"),
                                         helpText("For the next figures, you may optionally choose from predetermined genesets. First select a source.")),
                                  column(4,
                                         h4("2"),
                                         h6("Search gene sets"),
                                         helpText("Please note: selecting a pahtway will overwrite all genes in the set on the right. Any manual edits will be lost.")),
                                  column(4,
                                         h4("3"),
                                         h6("Check genes"),
                                         helpText("You may add or remove genes to this list."),
                                         helpText("All changes are immediately saved."),
                                         helpText("Take care: spaces are not allowed."))),
                         fluidRow(column(12,
                                         hr())),
                         fluidRow(column(4,
                                         radioGroupButtons(
                                             inputId = "geneset_source",
                                             label = "",
                                             choices = c("GO Biological Process", 
                                                         "KEGG", 
                                                         "DisGeNET"),
                                             direction = "vertical",
                                             justified = TRUE
                                         )),
                                  column(4,
                                         selectizeInput(
                                             inputId = "geneset_search",
                                             label = "", 
                                             choices = c(c("")),
                                             options = list(
                                                 maxItems = 1,
                                                 placeholder = "Start typing to search...",
                                                 `live-search` = TRUE)
                                         )
                                  ),
                                  column(4,
                                         textAreaInput(
                                             "geneset_genes",
                                             "",
                                             value = geneset_sample_genes(),
                                             height = "300px"
                                         ),
                                         tableOutput("geneset_errors"),
                                         actionButton(
                                             "geneset_reset",
                                             " Reset to standard geneset",
                                             icon = icon("recycle")
                                         )
                                  )
                         )
                ),
                tabPanel("6. Heatmap",
                         fluidRow(
                             column(
                                 3,
                                 helpText(
                                     "This plot is based on the geneset chosen in tab 5."
                                 ),
                                 hr(),
                                 sliderInput(
                                     "heatmap_intensity",
                                     "Log2 color range (mean: 0)",
                                     min = -10,
                                     max = 10,
                                     value = c(-6, 6)
                                 ),
                                 hr(),
                                 awesomeCheckbox(
                                     inputId = "heatmap_cluster_rows",
                                     label = "Cluster rows",
                                     value = TRUE
                                 ),
                                 awesomeCheckbox(
                                     inputId = "heatmap_cluster_cols",
                                     label = "Cluster columns",
                                     value = TRUE
                                 ),
                                 #hr(),
                                 #actionButton(
                                 #    "heatmap_generate",
                                 #    "Generate heatmap",
                                 #    icon = icon("rocket")
                                 #)
                             ),
                             column(
                                 9,
                                 tableOutput("heatmap_errors"),
                                 tableOutput("data_warnings_heatmap") %>% withSpinner(color = "lightgrey", proxy.height = "100px"),
                                 plotOutput("heatmap", height =
                                                "600px") %>% withSpinner(color = "lightgrey")
                             ),
                             column(
                                 12,
                                 downloadButton("download_heatmap", "Download")
                             )
                         )),
                tabPanel("7. Boxplot",
                         fluidRow(
                             column(
                                 3,
                                 #textAreaInput(
                                 #    "boxplot_genes",
                                 #    "Geneset",
                                 #    value = heatmap_sample_genes(),
                                 #    height = "400px"
                                 #),
                                 helpText(
                                     "This plot is based on the geneset chosen in tab 5."
                                 ),
                                 hr(),
                                 #actionButton(
                                 #    "boxplot_generate",
                                 #    "Generate boxplot",
                                 #    icon = icon("rocket")
                                 #)
                             ),
                             column(
                                 9,
                                 tableOutput("boxplot_errors"),
                                 tableOutput("data_warnings_boxplot") %>% withSpinner(color = "lightgrey", proxy.height = "100px"),
                                 plotOutput("boxplot", height =
                                                "600px") %>% withSpinner(color = "lightgrey")
                             ),
                             column(
                                 12,
                                 downloadButton("download_boxplot", "Download")
                             )
                         ))
            )
        }
    })
}

## Start app
shinyApp(ui = ui, server = server)

## Deply app
# deployApp(appTitle="HLC_CompR", launch.browser=T, logLevel = "verbose")

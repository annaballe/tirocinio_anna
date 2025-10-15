# ###make reactive database given reference genome
# txdb= reactive({if (input$UCSC_Genome == "hg19"){
#   txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene}
#   else{
#     txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene}
# })

# #################make reactive dataframe given input file

# mydata <- reactiveVal()
# observeEvent(input$file1, {

#   # Read uploaded table (user says this is already like a data frame)
#   tmp <- read.table(input$file1$datapath,
#                     header = input$header, stringsAsFactors = FALSE)

#   # Ensure first three columns are named consistently
#   if (ncol(tmp) >= 3) {
#     colnames(tmp)[1:3] <- c("chromosome","start","end")
#     head(tmp)
#   }

#   # If the uploaded table already has sample_* columns, just use it as-is
#   if (any(grepl("^sample_", names(tmp)))) {
#     mydata(tmp)
#     return(invisible(NULL))
#   }

#   # Otherwise, if it has a 4th column as coverage, map it to one sample
#   sample_label <- tools::file_path_sans_ext(basename(input$file1$name))
#   cov <- if (ncol(tmp) >= 4) suppressWarnings(as.integer(tmp[[4]])) else rep(NA_integer_, nrow(tmp))
#   new_df <- data.frame(
#     chromosome = tmp$chromosome,
#     start = as.integer(tmp$start),
#     end = as.integer(tmp$end),
#     cov = cov,
#     stringsAsFactors = FALSE
#   )
#   colnames(new_df)[4] <- paste0("sample_", sample_label)

#   mydata(new_df)
# })

# observeEvent(input$pileup, {

#   tmp_pileup <- coverage_input()
#   colnames(tmp_pileup)[1:3] <- c("chromosome","start","end") #nucleotide

#   # Determine how many coverage columns we have and how many samples
#   data_cols <- setdiff(colnames(tmp_pileup), c("chromosome","start","end"))
#   ncols <- length(data_cols)
#   samples <- name_sample()
#   k <- length(samples)

#   if (ncols == k && k > 0) {
#     # One coverage column per sample: name them as sample_* only
#     colnames(tmp_pileup)[-1:-3] <- paste0("sample_", samples)
#   } else if (ncols == 2 * k && k > 0) {
#     # Two columns per sample (legacy BAM): interleave as sample_*, nucleotide_* (kept for backward compatibility)
#     interleaved <- as.vector(rbind(paste0("sample_", samples)))#, paste0("nucleotide_", samples)))
#     colnames(tmp_pileup)[-1:-3] <- interleaved
#   } else {
#     # Fallback: generic sample_1..sample_n
#     colnames(tmp_pileup)[-1:-3] <- paste0("sample_", seq_len(ncols))
#   }

#   mydata(tmp_pileup)
#   print(head(tmp_pileup))
# })

# mysample<-reactive({
#   dat <- mydata()
#   if (is.null(dat)) return(NULL)
#   num <- input$Sample
#   i <- paste0("sample_", num)
#   if (!(i %in% names(dat))) return(NULL)

#   out <- dat[, c("chromosome","start","end", i), drop = FALSE]
#   names(out)[4] <- "coverage"
#   out
# })

# filtered_low<- reactive ({

#   if (is.null(mysample()))
#     return(NULL)
#   df <- mysample()
#   thr <- input$coverage_co
#   if (identical(thr, "all")) {
#     dplyr::filter(df, chromosome == Chromosome())
#   } else {
#     dplyr::filter(df, chromosome == Chromosome(), coverage <= as.numeric(thr))
#   }
#   print(filtered_low)

# })

# filtered_high<- reactive ({
#   if (is.null(mysample()))
#     return(NULL)
#   df <- mysample()
#   thr <- input$coverage_co
#   if (identical(thr, "all")) {
#     dplyr::filter(df, chromosome == Chromosome())
#   } else {
#     dplyr::filter(df, chromosome == Chromosome(),
#                   coverage > as.numeric(thr))
#   }
# })
mydata <- reactiveVal()
observeEvent(input$file1, {
 
  tmp <- read.table(input$file1$datapath,
                    header = input$header, stringsAsFactors = FALSE)
  colnames(tmp)[1:3]=c("chromosome","start","end")
  #n=  length(colnames(tmp)[-1:-3])
  #m<-rep(name_sample(), each=2)
  #colnames(tmp)[-1:-3]=paste0(c("sample_","nucleotide_"), m[1:n])
  head(tmp)
  ## do whatever is needed to parse the data
  mydata(tmp)
})
 
observeEvent(input$pileup, {
 
  tmp_pileup <- coverage_input()
  colnames(tmp_pileup)[1:3]=c("chromosome","start","end")
  n=  length(colnames(tmp_pileup)[-1:-3])
  #a=rep(head(seq_along(tmp_pileup),-3), each=2)
  m<-rep(name_sample(), each=2)
  #colnames(tmp_pileup)[-1:-3]=paste0(c("sample_","nucleotide_"), a[1:n])
  colnames(tmp_pileup)[-1:-3]<-paste0(c("sample_"), m[1:n])
  ## do whatever is needed to parse the data
  mydata(tmp_pileup)
  print(head(tmp_pileup))
})
 
 
 
mysample<-reactive({
  if (is.null(mydata()))
    return(NULL)
  num= input$Sample
  #num<- name_sample()
  i=paste0("sample_",num)
  nucleodites= paste0("nucleotide_",num)
  mydata() %>%
    dplyr:: select(chromosome, start, end,i,nucleodites) %>%
    dplyr::rename(coverage=i) %>%
    dplyr::rename(counts= nucleodites)
})
 
 
filtered_low<- reactive ({
  if (is.null(mysample()))
    return(NULL)
  mysample() %>%
    dplyr::select(-c(counts)) %>%
    dplyr::filter(chromosome == Chromosome(),
                  coverage <= as.numeric(input$coverage_co))
 
})
 
filtered_high<- reactive ({
  browser()
  if (is.null(mysample()))
    return(NULL)
  mysample() %>%
    dplyr::select(-c(counts)) %>%
    dplyr::filter(chromosome == Chromosome(),
                  coverage > as.numeric(input$coverage_co))
})
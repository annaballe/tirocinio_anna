require(OrganismDbi)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(TxDb.Hsapiens.UCSC.hg38.knownGene)
require(org.Hs.eg.db)
require(rlist)
require(Rsamtools)


ff=reactiveValues()
gene_list <- reactive ({
   file_gene<-input$gene1
   if (is.null(input$gene1)) {return(NULL)}
   if(input$type_file=="target_bed"){
      tmp_gene1<-read.table(input$gene1$datapath, stringsAsFactors = F)
      tmp_gene<- tmp_gene1[1:4]
      colnames(tmp_gene)<- c('chr', 'start', 'end', 'SYMBOL')
      #print(head(tmp_gene))
      ff<- tmp_gene }else{
         tmp_gene<- scan(file_gene$datapath, character(), quote = "")
         ff<- tmp_gene}
   #print(tmp_gene)
})

# list_bam<- reactive({
#    if (is.null(input$bam_list)) {return(NULL)}
#    tmp_bam <- scan(input$bam_list$datapath, character(), quote = "")
#    #list_bam= (tmp_bam)
#    #print(tmp_bam)
# })
# Lista file di coverage (bam o bed)
list_coverage <- reactive({
    if (is.null(input$bam_list)) return(NULL)
    scan(input$bam_list$datapath, character(), quote = "")
})


all_gene<- reactive({ if (input$Genome == "hg19"){
   all_gene<- TxDb.Hsapiens.UCSC.hg19.knownGene}else{
      all_gene<-TxDb.Hsapiens.UCSC.hg38.knownGene}
})

no_entrID<- reactive({
   if (is.null(gene_list()))
      return(NULL)
   if (is.null(list_coverage()))
      return(NULL)
   if(input$type_file=="target_bed"){
      errore<- data.frame(matrix(ncol = 2, nrow = 0))
   }else{# STAGE 1: Convert gene names to ENTREZ IDs using org.Hs.eg.db
    cat("Stage 1: Converting gene names to ENTREZ IDs using org.Hs.eg.db...\n")

    # Handle 1:many mappings by taking the first ENTREZ ID for each gene
    my_gene_name <- OrganismDbi::select(org.Hs.eg.db, key = gene_list(), columns = "ENTREZID", keytype = "ALIAS") %>%
      tibble::as_tibble() %>%
      dplyr::filter(!is.na(ENTREZID)) %>%
      dplyr::group_by(ALIAS) %>%
      dplyr::slice(1) %>%  # Take first ENTREZ ID if multiple exist
      dplyr::ungroup()

    # Check for genes not found in org.Hs.eg.db
    found_genes_stage1 <- unique(my_gene_name$ALIAS)
    not_found_stage1 <- base::setdiff(gene_list(), found_genes_stage1)

    if (length(not_found_stage1) > 0) {
      cat("WARNING: These genes were not found in org.Hs.eg.db and will be excluded:\n")
      for(gene in not_found_stage1) {
      cat(paste("  - Gene '", gene, "' not found in org.Hs.eg.db\n", sep=""))
    }
  cat("\n") #
}

if (nrow(my_gene_name) == 0) {
  stop("ERROR: No valid genes found in org.Hs.eg.db. Please check your gene names.")
}

cat(paste("Stage 1 result:", nrow(my_gene_name), "gene-ID mappings found from", length(found_genes_stage1), "unique genes\n\n"))

# STAGE 2: Check if ENTREZ IDs exist in genome reference database
cat("Stage 2: Validating genes against reference genome database...\n")

all_entrez_ids <- unique(my_gene_name$ENTREZID)
genome_keys <- AnnotationDbi::keys(all_gene(), "GENEID")  # Fixed: added () to call the reactive

# Find genes present in both databases
valid_entrez_ids <- base::intersect(all_entrez_ids, genome_keys)
invalid_entrez_ids <- base::setdiff(all_entrez_ids, genome_keys)

if (length(invalid_entrez_ids) > 0) {
  # Get gene names for invalid ENTREZ IDs
  invalid_genes_stage2 <- my_gene_name %>% 
    dplyr::filter(ENTREZID %in% invalid_entrez_ids) %>%
    dplyr::pull(ALIAS) %>%
    unique()
  
  cat("WARNING: These genes were found in org.Hs.eg.db but not in the reference genome database and will be excluded:\n")
  for(gene in invalid_genes_stage2) {
    cat(paste("  - Gene '", gene, "' not found in reference genome database (", input$Genome, ")\n", sep=""))  # Fixed: changed 'genome' to 'input$Genome'
  }
  cat("\n")
} else {
  invalid_genes_stage2 <- character(0)
}

# Filter to keep only valid genes
my_gene_name <- my_gene_name %>%
  dplyr::filter(ENTREZID %in% valid_entrez_ids)

ID <- valid_entrez_ids

if (length(ID) == 0) {
  stop("ERROR: No genes passed both validation stages. Please check your gene names and genome version.")
}

cat(paste("Stage 2 result:", length(valid_entrez_ids), "genes validated against reference genome\n"))
cat(paste("FINAL RESULT: Proceeding with", length(ID), "validated genes for analysis\n"))

# Print summary of excluded genes
total_excluded <- length(not_found_stage1) + length(invalid_genes_stage2)
if (total_excluded > 0) {
  cat(paste("SUMMARY: ", total_excluded, " genes excluded (", 
            length(not_found_stage1), " from stage 1, ", 
            length(invalid_genes_stage2), " from stage 2)\n", sep=""))
}
cat("=== GENE VALIDATION COMPLETE ===\n\n")

# Return error data frame if there are issues, otherwise return empty data frame
if (length(ID) == 0) {
  errore <- data.frame(
    Gene = c("No valid genes found"),
    Issue = c("All genes failed validation"),
    stringsAsFactors = FALSE
  )
  return(errore)
} else {
  errore <- data.frame(matrix(ncol = 2, nrow = 0))
  return(errore)
}
  }
})



#       gene_list1<- gene_list()
#       my_gene_name<- OrganismDbi::select(org.Hs.eg.db, key= gene_list1, columns="ENTREZID",
#                                          keytype="ALIAS")
#       #print(my_gene_name)
#       ID<- my_gene_name$ENTREZID
#       b<- c()
#       for (i in ID){
#          if( ! is.element(i, keys(all_gene(), "GENEID")))
#             b[i]<- i
#       }
#       errore<- subset(my_gene_name, my_gene_name$ENTREZID %in% b)
#       colnames(errore)[1]= " the following HGNC official gene names are
#     unrecognizable, please choose a correct name and reload a file"
#       return(errore)
#       #print(head(errore,n=12))
#    }
# })


for_bed<- reactive({
   if (is.null(gene_list()))
      return(NULL)
   if (is.null(list_coverage()))
      return(NULL)
   if(nrow(no_entrID())!=0)
      return(no_entrID())
   if(input$type_file=="target_bed"){
      for_bed<- gene_list() }else{
         gene_list1<- gene_list()
         
         # Get valid gene mappings (same logic as in no_entrID but simplified)
         my_gene_name <- OrganismDbi::select(org.Hs.eg.db, key = gene_list1, columns = "ENTREZID", keytype = "ALIAS") %>%
           tibble::as_tibble() %>%
           dplyr::filter(!is.na(ENTREZID)) %>%
           dplyr::group_by(ALIAS) %>%
           dplyr::slice(1) %>%  # Take first ENTREZ ID if multiple exist
           dplyr::ungroup()
         
         # Get valid ENTREZ IDs that exist in the genome database
         all_entrez_ids <- unique(my_gene_name$ENTREZID)
         genome_keys <- AnnotationDbi::keys(all_gene(), "GENEID")
         valid_entrez_ids <- base::intersect(all_entrez_ids, genome_keys)
         
         # Filter to keep only valid genes
         my_gene_name <- my_gene_name %>%
           dplyr::filter(ENTREZID %in% valid_entrez_ids)
         
         ID <- valid_entrez_ids
         pre<- data.frame()
         
         for (i in ID){
            # Check if the ENTREZ ID is valid before using it
            if (i %in% genome_keys) {
               txid <- OrganismDbi::select(all_gene(), i, "TXNAME", "GENEID")[["TXNAME"]]
               cds <-OrganismDbi::exonsBy(all_gene(), by="tx", use.names=TRUE)
               coor<- as.data.frame(cds[names(cds) %in% txid])
               coordinate<- cbind(coor,i)
               colnames(coordinate)[11]<- "ENTREZID"
               cood<- dplyr::inner_join(coordinate, my_gene_name, by="ENTREZID")
               pre<- rbind(cood,pre)
               pre$start<- pre$start -10
               pre$end<- pre$end +10
            }
         }
         #print(head(pre))
         for_bed<- data.frame()
         for (row in 1:nrow(pre)){
            sel_cc<- data.frame(as.character(pre$seqnames[row]),
                                pre$start[row], pre$end[row],pre$ALIAS[row],
                                stringsAsFactors = FALSE)
            for_bed<- rbind(for_bed, sel_cc)

         }
         colnames(for_bed)<- c('chr', 'start', 'end', 'SYMBOL')
         for_bed<- unique(for_bed)
         if(input$notation == "number"){
            for_bed$chr<- gsub("^.{0,3}", "", for_bed$chr, perl =TRUE)}
            
         return(for_bed)
         print(for_bed)
      }
})


coverage_input <- reactive({
    if (is.null(gene_list())) return(NULL)
    if (!is.null(no_entrID()) && nrow(no_entrID()) != 0) return(no_entrID())
    if (is.null(list_coverage())) return(NULL)

    if (input$type_coverage == "bam") {
        # --- calcolo coverage dai BAM ---
        for_grange <- GenomicRanges::makeGRangesFromDataFrame(for_bed(), keep.extra.columns = TRUE)
        print("for_grange") 
        print(for_grange)
        param <- Rsamtools::ScanBamParam(which=for_grange)
        paste("param")
        print(param)

        p_param <- Rsamtools::PileupParam(
            distinguish_nucleotides=FALSE,
            distinguish_strands=FALSE,
            min_mapq=as.numeric(input$MAPQ),
            min_base_quality=as.numeric(input$base_qual),
            min_nucleotide_depth=1,
            max_depth=150000
        )

        df <- list()
        for (i in list_coverage()) {
            pileup_df <- Rsamtools::pileup(i, scanBamParam=param, pileupParam=p_param)
            df <- rlist::list.append(df, pileup_df)
        }

        lst1 <- lapply(df, function(x) transform(x[,-5]))
        lst2 <- lapply(lst1, function(x) transform(x[!duplicated(x),]))
        riarrange.df <- function(list_df) {
            list_df %>%
                dplyr::mutate(end=pos) %>%
                dplyr::group_by(seqnames,pos,end) %>%
                dplyr::summarise(count=sum(count)) %>%
                dplyr::arrange()
        }
        lst3 <- lapply(lst2, riarrange.df)
        pp <- Reduce(function(...) merge(..., by=c("seqnames","pos","end")), lst3)
        pp[is.na(pp)] <- 0
        colnames(pp)[1:3] <- c("seqnames","start","end")
        g<-pp
         print("g")
        print(g)
        return(pp)
    }
    if (input$type_coverage == "bed") {
      # 1. Get target gene regions (same as BAM)
       for_grange <- GenomicRanges::makeGRangesFromDataFrame(
          for_bed(), 
           keep.extra.columns = TRUE
      )
    
    # 2. Read and filter each BED file
       df_list <- list()
      for (file in list_coverage()) {
        
          cat("\n=== READING BED FILE ===\n")
         cat(paste("File:", file, "\n"))
        
        # Read BED - assuming 4 columns: chr, start, end, coverage
         # Don't assume header exists in BED files
         df <- read.table(file, 
                         header = FALSE,
                         stringsAsFactors = FALSE,
                         colClasses = c("character", "integer", "integer", "numeric"))
        
        # Assign column names
           colnames(df) <- c("seqnames", "start", "end", "count")

          cat("Before filtering - rows:", nrow(df), "\n")
           print(head(df))
        
        # Convert to GRanges for filtering
         df_gr <- GenomicRanges::makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)
         # Right after creating df_gr and before filtering
         cat("\n=== COORDINATE COMPARISON ===\n")
         cat("BED Coverage Ranges:\n")
         cat("  Chromosome:", as.character(unique(seqnames(df_gr))), "\n")
         cat("  Min position:", min(start(df_gr)), "\n")
         cat("  Max position:", max(end(df_gr)), "\n")
         cat("  Example regions:\n")
         print(head(df_gr, 3))

         cat("\nGene Target Ranges:\n")
         cat("  Chromosome:", as.character(unique(seqnames(for_grange))), "\n")
         cat("  Min position:", min(start(for_grange)), "\n")
         cat("  Max position:", max(end(for_grange)), "\n")
         cat("  Example regions:\n")
         print(head(for_grange, 3))

         cat("\n=== CHECKING CHROMOSOME NAMES ===\n")
         cat("BED seqlevels:", seqlevels(df_gr), "\n")
         cat("Gene seqlevels:", seqlevels(for_grange), "\n")
         cat("Do they match?", any(seqlevels(df_gr) %in% seqlevels(for_grange)), "\n")
         cat("============================\n\n")

         # Now do the filtering
         df_gr_filtered <- subsetByOverlaps(df_gr, for_grange, type="any")

        
           # Filter to only gene regions (like ScanBamParam does for BAM)
         df_gr_filtered <- IRanges::subsetByOverlaps(df_gr, for_grange, type="any")
        
          cat("After filtering - regions:", length(df_gr_filtered), "\n")
        
          # Convert back to data frame
           df_filtered <- as.data.frame(df_gr_filtered)
        
        # Keep only the columns we need: seqnames, start, end, count
        df_filtered <- df_filtered[, c("seqnames", "start", "end", "count")]
        
         cat("Filtered data frame structure:\n")
         print(str(df_filtered))
         print(head(df_filtered))
         cat("===========================\n\n")
        
           df_list <- append(df_list, list(df_filtered))
       }
    
    # 3. Merge all samples (exactly like BAM processing!)
       pp <- Reduce(function(...) merge(..., by=c("seqnames","start","end"), all=TRUE), df_list)
     pp[is.na(pp)] <- 0
    
    # 4. CRITICAL: Convert types to match BAM output exactly
    # BAM produces: seqnames=factor, start=int, end=int, count=int
       pp$seqnames <- factor(pp$seqnames)
       pp$start <- as.integer(pp$start)
       pp$end <- as.integer(pp$end)
    
    # Convert all count columns to integer
       count_cols <- grep("^count", colnames(pp), value = TRUE)
     pp[count_cols] <- lapply(pp[count_cols], as.integer)
    
    # Verify output structure matches BAM
       cat("\n=== FINAL OUTPUT STRUCTURE ===\n")
       cat("Column types:\n")
       print(sapply(pp, class))
       cat("First rows:\n")
       print(head(pp))
       cat("===============================\n\n")
    
      g <- pp
      print("g")
      print(g)
    
      return(pp)
   }



        
    
})

# pileup_input<- reactive({
#    if (is.null(gene_list()))
#       return(NULL)
#    if (is.null(list_bam()))
#       return(NULL)
#    if(nrow(no_entrID())!=0)
#       return(no_entrID())

#    for_grange<-GenomicRanges::makeGRangesFromDataFrame(for_bed(),
#                                                        keep.extra.columns = TRUE)
#    #print(head(for_grange))
#    param <- Rsamtools::ScanBamParam(which= for_grange)

#    p_param <- Rsamtools::PileupParam(distinguish_nucleotides=FALSE,
#                                      distinguish_strands=FALSE,
#                                      min_mapq=as.numeric(input$MAPQ),
#                                      min_base_quality=as.numeric(input$base_qual),
#                                      min_nucleotide_depth=1,
#                                      max_depth=150000)
#    df= list()
#    for (i in list_bam()){
#       pileup_df= Rsamtools::pileup(i, scanBamParam=param, pileupParam=p_param)

#       df=rlist::list.append(df, pileup_df)
#    }
#    lst1 <- lapply(df, function(x) transform(x[,-5]))
#    lst2<- lapply(lst1, function(x) transform(x[!duplicated(x),]))

#    riarrange.df <- function(list_df){
#       require(dplyr)
#       list_df %>%
#          dplyr::mutate(end= pos) %>%
#          dplyr::group_by(seqnames,pos,end) %>%
#          dplyr::summarise(count=sum(count)) %>%
#          dplyr::arrange() }

#          #dplyr::summarise(value= as.numeric(paste(sum(count))),
#                           #counts=paste(count, sep=':',collapse=';'
   
#    lst3<- lapply(lst2, riarrange.df)

#    pp<- Reduce(function(...) merge(...,by= c("seqnames", "pos", "end")), lst3)
#    pp[is.na(pp)] <- 0
#    pp<- as.data.frame(pp)
#    if(input$notation == "number"){
#       for (i in pp[1]){
#          Chromosome<-paste("chr", i, sep="")
#          pp<- cbind(Chromosome, pp)
#          pp[,2]<-NULL
#       }
#    }
#    colnames(pp)<-NULL
#    return(pp)
# })


name_sample<- reactive({
   c<-gsub(".*/","",list_coverage())
   return(c)
})

stat_summ<- reactive({
   if (is.null(gene_list()))
      return(NULL)
   if (is.null(list_coverage()))
      return(NULL)
   if(nrow(no_entrID())!=0)
      return(no_entrID())
   ppinp<- as.data.frame(coverage_input())
   print(ppinp)
   colnames(ppinp)[1:3]<-c("seqnames","start","end")
   n<- length(colnames(ppinp)[-1:-3])
   #c<-gsub(".*/","",list_bam())
   m<-rep(name_sample(), each=2)
   colnames(ppinp)[-1:-3]<-paste0(c("sample_","nucleotide_"), m[1:n])

   #colnames(ppinp)[-1:-2]<-paste0(c("sample_"), m[1:n])
   for_grange<-GenomicRanges::makeGRangesFromDataFrame(for_bed(),
                                                       keep.extra.columns = TRUE)

   for_range_pp=GenomicRanges::makeGRangesFromDataFrame(ppinp,
                                                        keep.extra.columns = TRUE)
   tp<- GenomicRanges::findOverlaps(query= for_range_pp,
                                    subject = for_grange, type="any",
                                    select = "all")
   sts_df <- data.frame(for_range_pp[queryHits(tp),], for_grange[subjectHits(tp),])
   print(sts_df)
   ### Code diagnostics 7/10/2025
   # Dopo: sts_df <- data.frame(for_range_pp[queryHits(tp),], for_grange[subjectHits(tp),])
   cat("\n=== OVERLAP DIAGNOSTICS ===\n")
   cat(paste("Number of overlaps found:", length(tp), "\n"))
   cat(paste("Unique genes with coverage:", length(unique(sts_df$SYMBOL)), "\n"))
   cat("Genes WITH coverage:\n")
   print(unique(sts_df$SYMBOL))
   cat("Genes WITHOUT coverage:\n")
   genes_no_coverage <- setdiff(unique(for_bed()$SYMBOL), unique(sts_df$SYMBOL))
   print(genes_no_coverage)
   cat("Coverage chromosomes:\n")
   print(unique(ppinp$seqnames))
   cat("Gene chromosomes:\n")
   print(table(for_bed()$chr))
   cat("===========================\n\n")
   ### end

   statistiche<- sts_df[!duplicated(sts_df$start),]
   statistiche<- subset(statistiche,
                        select = -c(width, strand, seqnames.1,
                                    start.1, end.1, width.1, strand.1))

   colnames(statistiche)[1:3]<- c("chromosome","start","end")
   merge_g<- dplyr::full_join(for_bed(),statistiche, by="SYMBOL")
   #merge_g <- dplyr::full_join(for_bed(), statistiche, by="SYMBOL")
   
   #col.sub<- colnames(merge_g[, grepl("sample_" , names(merge_g))])
   col_name= colnames(merge_g)

   col.sub= col_name[grepl("sample_", col_name)]
  merge_g[col.sub] <- sapply(merge_g[col.sub],as.numeric)

   x<- list()
   for (i in col.sub){
      x_df<- merge_g %>%
      dplyr::select("chromosome","start.x","end.x", "SYMBOL",i)
      #print(x_df)
      colnames(x_df)[5]<- "value"
      x_df$sample<- paste0(i)
      x<- list.append(x, x_df)
   }
   statistical_operation= function(df){
      df %>%
         dplyr::group_by(SYMBOL,sample) %>%
         dplyr::summarize(Total= sum(value, na.rm=TRUE),
                          Mean = mean(value, na.rm=TRUE),
                          Median= median(value, na.rm=TRUE),
                          number_of_position_under_20x = sum(value < 20),
                          percentage_under_20x= (sum(value < 20)/sum(value, na.rm=TRUE))*100) %>%
         as.data.frame(.) %>%
         dplyr::mutate_if(is.numeric, round, 3)
   }

   out_r<- do.call(rbind, lapply(x, function(x) statistical_operation(x)))
   return(out_r)


})

#
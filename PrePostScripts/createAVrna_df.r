suppressMessages(library(VariantAnnotation))
suppressMessages(library(stringr))
suppressMessages(library(vcfR))
suppressMessages(library(jsonlite))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript script.R path_to_hotspots.json")
json_path <- args[1]

# Use current working directory as base
base_dir <- getwd()
message("Using base directory: ", base_dir)


# ---------------------------
# Utilities / File Reading
# ---------------------------

get_sample_ids <- function(path) {
  stopifnot(dir.exists(path))
  files <- list.files(path, pattern = "_mpileup_annot.vcf$", full.names = FALSE)
  str_replace(files, "_mpileup_annot.vcf$", "")
}

load_hotspots_json <- function(json_path) {
  if (!file.exists(json_path)) stop("Hotspot JSON file not found: ", json_path)
  fromJSON(json_path, simplifyVector = FALSE)
}

hotspots_json <- load_hotspots_json(json_path)

hotspots_p =  unlist(lapply(hotspots_json, function(x) if(!is.null(x$protein_hotspots)) x$protein_hotspots else character(0)))
hotspots_c = unlist(lapply(hotspots_json, function(x) if(!is.null(x$nucleotide_hotspots)) x$nucleotide_hotspots else character(0)))


# ---------------------------
# VCF Extraction & Filtering
# ---------------------------

extract_hotspot_data <- function(vcf_path, protein_hotspots, nucleotide_hotspots, hotspots_json, sample_id) {
  print(sample_id)
  tryCatch({
    vcf <- readVcf(vcf_path)

     if (dim(info(vcf))[1] == 0) {
       na_df <- data.frame(matrix(NA, nrow = 1, ncol = 45, dimnames = list(NULL, 
        c("Row.names", "seqnames", "start",  "end", "width", "strand", "paramRangeID", "QUAL", "Depth", "Read_Alt", 
          "Chrom", "Mut", "Annotation", "Impact", "Eval", "Gene_ID", "Feature_Type", "Feature", "Biotype", 
           "Rank", "HGVS_c", "HGVS_p", "cDNA_pos", "CDS_pos", "Protein_pos", "Allele", "RefAlt", "Strand", 
            "Consequence", "Variant_Classification", "Other", "Other_Info", "mutationType", "NA1", 
          "NA2", "NA3", "NA4", "Tag", "refSeqID", "EnsemblProtID", "NA5", "ProteinEffect", "NA6", "NA7","Sample_ID"))))
      na_df$Sample_ID <- sample_id
      return(na_df)
    }

    annotglob = as.data.frame(vcf@rowRanges)
    annotglob$QUAL = rowRanges(vcf)$QUAL
    annotglob$Depth = info(vcf)$DP
    
    annotglob$Read_Alt = sapply(seq_along(info(vcf)$DP4), function(gene_idx) sum(info(vcf)$DP4[[gene_idx]][3:4]))
    annotglob$Chrom= sapply(rownames(annotglob), strsplit, split = ":") %>% sapply(function(x) x[1])
    annotglob$Mut = sapply(rownames(annotglob), strsplit, split = "_") %>% sapply(function(x) x[2])


    csq_info <- info(vcf)$CSQ
    
      
     csq_df <- do.call(rbind, lapply(csq_info, function(csq) {
      csq_records <- strsplit(as.character(csq), ",")
      tmp <- do.call(rbind, lapply(csq_records, function(record) {
        csq_split <- strsplit(as.character(record), "\\|")[[1]]
        if(length(csq_split)==31){
          csq_split[32] <- ""
        }
        data.frame(matrix(csq_split, ncol = length(csq_split), byrow = TRUE), stringsAsFactors = FALSE)
      }))
    }))
    colnames(csq_df) <- c("Annotation", "Impact", "Eval", "Gene_ID", "Feature_Type", "Feature", "Biotype",
                          "Rank", "HGVS_c", "HGVS_p", "cDNA_pos", "CDS_pos", "Protein_pos", "Allele", "RefAlt",
                          "Strand", "Consequence", "Variant_Classification", "Other", "Other_Info", 
                          "mutationType", "NA1", "NA2", "NA3", "NA4", "Tag", "refSeqID", "EnsemblProtID", "NA5", "ProteinEffect", "NA6", "NA7")
    csq_df$corresp <- rep(unlist(info(vcf)@rownames),  sapply(as.list(info(vcf)$CSQ),length))

    csq_df$HGVS_c <- sapply(csq_df$cDNA_pos, function(hgvs) if (!is.na(hgvs)) strsplit(hgvs, ":")[[1]][2] else NA)
    csq_df$HGVS_p <- sapply(csq_df$CDS_pos, function(hgvs) if (!is.na(hgvs)) strsplit(hgvs, ":")[[1]][2] else NA)

    mergedf = merge(annotglob, csq_df, by.x = "row.names", by.y = "corresp", all.x = TRUE)
        
    posstart = hotspots_json$TP53$range$start
    posend = hotspots_json$TP53$range$end

    filtered_df <- mergedf[!is.na(mergedf$HGVS_c) %in% nucleotide_hotspots | (mergedf$start > posstart & mergedf$start < posend),]

    if (nrow(filtered_df) == 0) {
      df <- mergedf[0, ]  # même colonnes
      df[1, ] <- NA
      df$Sample_ID <- sample_id
      return(df)
    }
    
    filtered_df$Sample_ID <- sample_id  
    return(filtered_df)
    
  }, error = function(e) {
    
    message("⚠️ Erreur sur ", sample_id, ": ", e$message)
    return(data.frame(Sample_ID = sample_id, stringsAsFactors = FALSE))
  })
}

# sample_idsIPMN = c("BIPMN_017_4114_L1_S40_L001", "BIPMN_130_20123_L2" ,  "BIPMN_00513_07_L1_S18", "BIPMN_00513_14_L1_S19", "BIPMN_106_20529_L1","BIPMN_137_21613_L1", "BIPMN_137_21613_L2")    
sample_ids <- get_sample_ids(file.path(base_dir, "VEP_output"))

vcf_df <- do.call(rbind, lapply(sample_ids, function(sample_id) {
  vcf_path <- file.path("./VEP_output", paste0(sample_id, "_mpileup_annot.vcf"))
 extract_hotspot_data(vcf_path, hotspots_p, hotspots_c, hotspots_json, sample_id)
}))
vcf_df = vcf_df[(vcf_df$Gene_ID %in% c("KRAS", "GNAS") & vcf_df$Tag == "MANE_Select") | (vcf_df$Gene_ID %in% c("TP53", "WRAP53") &  vcf_df$Eval == "HIGH" &  grepl("frameshift|stop_gained", vcf_df$Impact)),]


# ---------------------------
# Build result_list from Samtools depth
# ---------------------------
build_positions_from_json <- function(hotspots_json) {
  list(
    KRAS = if (!is.null(hotspots_json$KRAS$positions)) as.character(unlist(hotspots_json$KRAS$positions)) else character(0),
    GNAS = if (!is.null(hotspots_json$GNAS$positions)) as.character(unlist(hotspots_json$GNAS$positions)) else character(0),
    TP53 = if (!is.null(hotspots_json$TP53$range)) as.character(seq(hotspots_json$TP53$range$start, hotspots_json$TP53$range$end)) else character(0)
  )
}
posi <- build_positions_from_json(hotspots_json)

samtools_depth_files <- list.files("./Samtools_depth", pattern = "_depth.txt", full.names = TRUE)
non_empty_files <- samtools_depth_files[sapply(samtools_depth_files, function(f) file.info(f)$size > 0)]

result_list <- list()

for (file in samtools_depth_files) {
  file_name <- str_replace(basename(file), "_depth.txt", "")
  
  # Créer un dataframe avec toutes les positions depuis le JSON
  df_final <- data.frame(
    gene = rep(names(posi), times = sapply(posi, length)),
    position = unlist(posi),
    samtoolsDepth = NA,
    mut_pos = NA,
    mut_readAlt = NA,
    mut_Depth = NA,
    HGVS_c = NA
  )
  
  # Remplir la profondeur si le fichier n’est pas vide
  if (file %in% non_empty_files) {
    df <- read.delim(file, sep = "\t", header = FALSE, col.names = c("chr", "position", "samtoolsDepth"))
    df_final$samtoolsDepth <- df$samtoolsDepth[match(df_final$position, df$position)]
  }
  
  result_list[[file_name]] <- df_final
}

# Ajouter les informations de mutations depuis mut
for (sample_id in names(result_list)) {
  print("Add samtools depth results")
  if (sample_id %in% vcf_df$Sample_ID) {
    df_sample_mut <- vcf_df[vcf_df$Sample_ID == sample_id, ]
 
    idx <- match(as.numeric(result_list[[sample_id]]$position), df_sample_mut$start)
    
    print(head(idx))
    
    result_list[[sample_id]]$mut_pos <- df_sample_mut$start[idx]
    result_list[[sample_id]]$mut_readAlt <- df_sample_mut$Read_Alt[idx]
    result_list[[sample_id]]$mut_Depth <- df_sample_mut$Depth[idx]
    result_list[[sample_id]]$HGVS_c <- df_sample_mut$HGVS_c[idx]
  }
}

# ---------------------------
# Mutation statistics / summary
# ---------------------------
compute_gene_metrics <- function(df_gene, gene_name) {
  res <- list(
    is_mut = FALSE,
    nuclMut = NA_character_,
    nbReadMut = 0L,
    nbReadTot = 0L,
    VAF = NA_real_
  )
  
  if (nrow(df_gene) == 0) return(res)
  
  # Minimum depth filter
  mut_df <- df_gene[!is.na(df_gene$mut_pos) & df_gene$samtoolsDepth >= 10, ]

    if (nrow(mut_df) > 0) {
    res$nbReadTot <- sum(mut_df$mut_Depth, na.rm = TRUE)
  } else {
    res$nbReadTot <- sum(df_gene$samtoolsDepth, na.rm = TRUE)
  }
  res$nbReadMut <- sum(mut_df$mut_readAlt, na.rm = TRUE)
  res$nuclMut <- paste(unique(na.omit(mut_df$HGVS_c)), collapse = ";")

  # VAF calcul only when mutations are present
  if (res$nbReadMut > 0) {
    if (tolower(gene_name) == "TP53") {
      mut_df$VAF <- with(mut_df, ifelse(mut_Depth > 0, mut_readAlt / mut_Depth, NA))
      res$VAF <- max(mut_df$VAF, na.rm = TRUE)
    } else {
      res$VAF <- ifelse(res$nbReadTot > 0, res$nbReadMut / res$nbReadTot, NA_real_)
    }
  } else {
    res$VAF <- 0
  }
  # filtre is_mut
  if (res$nbReadMut > 0 && res$nbReadTot >= 10 && res$VAF >= 0.2) {
    res$is_mut <- TRUE
  }
  return(res)
}

genes <- c("KRAS", "GNAS", "TP53")
print("Add mutation conclusions")

summary_df <- data.frame(sample_id = names(result_list), stringsAsFactors = FALSE)
for (gene in genes) {
  for (col in c("is_mut", "nuclMut", "nbReadMut", "nbReadTot", "VAF")) {
    summary_df[[paste0(gene, "_", col)]] <- NA
  }
}

for (sample_id in names(result_list)) {
  df_sample <- result_list[[sample_id]]
  
  for (gene in genes) {
    df_gene <- df_sample[df_sample$gene == gene, ]
    metrics <- compute_gene_metrics(df_gene, gene)    
    for (col in names(metrics)) {
      summary_df[summary_df$sample_id == sample_id, paste0(gene, "_", col)] <- metrics[[col]]
    }
  }
}
write.table(summary_df, file = "mutation_summary.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

library(stringr)

get_ptm_position = function(sequence_in,protein_positions_in,peptide_position_in,probability_threshold) {
  sequence_start = as.numeric(str_split(protein_positions_in,";")[[1]][1]) - as.numeric(str_split(peptide_position_in,";")[[1]][1]) 
  sequence_location = 0
  ptm_site_location_peptide = numeric()
  sequence_in_split = unlist(str_split(sequence_in,""))
  for (i in 1:length(sequence_in_split)) {
    if (str_detect(sequence_in_split[i],"[A-Z]")) {
      sequence_location = sequence_location + 1
    }
    else if (str_detect(sequence_in_split[i],"\\(")) {
      ptm_site_location_peptide = c(ptm_site_location_peptide,sequence_location)
    }
  }
  aa_sequence = unlist(str_split(paste(unlist(str_split(sequence_in,"\\(([0-9]*\\.[0-9]*)\\)|\\(([0-9])\\)")),collapse=""),"")) #split by the probabilities in parenthesis and extract out the amino acids
  aa = aa_sequence[ptm_site_location_peptide] #get the modified amino acid
  modification_probs = as.numeric(unlist(str_extract_all(sequence_in,"([0-9]*\\.[0-9]*)|([0-9])")))
  ptm_site_location_protein = sequence_start + ptm_site_location_peptide
  modification_df = data.frame(aa=aa,ptm_site_location_peptide=ptm_site_location_peptide,ptm_site_location_protein=ptm_site_location_protein,modification_probs = modification_probs,stringsAsFactors=FALSE)
  modification_string = character()
  for (i in 1:nrow(modification_df)) {
    if (as.numeric(modification_df$modification_probs[i]) >= probability_threshold) {
      modification_string = c(modification_string,paste(c(modification_df[i,1],modification_df[i,3],"*"),collapse=""))
    } else {
      modification_string = c(modification_string,paste(c(modification_df[i,1],modification_df[i,3]),collapse=""))
    }
  }
  return(paste(modification_string,collapse=","))
}


average_replicates = function(df_in, metadata_in, missing_value) {
  sample_column = "sample"
  length_unique_sample_column = length(unique(metadata_in[, sample_column]))
  averaged_expression_matrix = matrix(data = NA, nrow = dim(df_in)[1], ncol = length_unique_sample_column)
  for (i in 1:length_unique_sample_column) {
    current_sample = unique(metadata_in[, sample_column])[i]
    current_sample_full_text = as.character(metadata_in[metadata_in[, sample_column] == current_sample, "expression_column"])
    current_sample_df = df_in[, current_sample_full_text]
    for (j in 1:nrow(current_sample_df)) {
      if (is.na(missing_value)) {
        current_peptide_missing = sum(is.na(current_sample_df[j, ]))
      } else {
        current_peptide_missing = sum(current_sample_df[j, ] == missing_value)
      }
      current_peptide_sum = sum(na.omit(as.numeric(current_sample_df[j, ])))
      current_sample_df_ncol = ncol(current_sample_df) 
      if (current_peptide_missing < current_sample_df_ncol) {
        averaged_expression_matrix[j, i] = current_peptide_sum/(current_sample_df_ncol - current_peptide_missing) # "average" the nonzero expression
      } else if (current_peptide_missing == current_sample_df_ncol) {
        averaged_expression_matrix[j, i] = NA
      }
    }
  }
  colnames(averaged_expression_matrix) = unique(metadata_in[, sample_column])
  return(log2(data.frame(averaged_expression_matrix, stringsAsFactors = FALSE)))
}


process_mq_data = function(df_in, protein_id_col, metadata_in, is_ptm_data = FALSE, ptm_threshold = 0.50, IRON = TRUE, find_median_col = "", missing_values = 0) {
  if (sum(c("expression_column", "sample", "replicate") %in% names(metadata_in)) != 3) {
    stop("Metadata is missing required columns (expression_column, sample, sample_fraction, or replicate)")
  }
  if (is_ptm_data == TRUE) {
    modification_prob_column = names(df_in)[str_detect(names(df_in), "Probabilities")]
    modification_count_column = names(df_in)[str_detect(names(df_in), "Number\\.of")]
  }
  
  #normalize data with IRON ----
  if (IRON == TRUE) {
    if (find_median_col == "") {
      ironed_df = iron_proteomics(df_in[, metadata_in$expression_column])
      metadata_in$pre_iron_expression = metadata_in$expression_column
      metadata_in$expression_column = names(ironed_df)
      df_in = cbind(df_in, ironed_df)
    }
  } else {
    ironed_df = iron_proteomics(df_in[, metadata_in$expression_column], median_sample = find_median_col)
    metadata_in$pre_iron_expression = metadata_in$expression_column
    metadata_in$expression_column = names(ironed_df)
    df_in = cbind(df_in, ironed_df)
  }
  
  #peptide_qc ----
  #flag non-human entries
  uniprot_entries = read.delim("data/uniprot_reviewed_human_refproteome_08-24-15.txt", header = T, stringsAsFactors = FALSE)
  df_in$qc_nothuman = sapply(df_in[, protein_id_col], function(x) {
    sum(unlist(str_extract_all(x, uniprot_regex(isoforms = "none"))) %in% uniprot_entries$Entry) == 0
  })
  
  #flag reverse only entries
  df_in$qc_reverse_only = sapply(df_in[, protein_id_col], function(x) {
    x_split = unlist(str_split(x, ";"))
    if (sum(str_detect(x_split, "REV__(CON__){0,1}|(REV__){0,1}CON__ENSEMBL|(REV__){0,1}CON__REFSEQ|(REV__){0,1}CON__H-INV")) == length(x_split)) {
      TRUE
    } else {
      FALSE
    }
  })
  
  #flag completely missing rows; na.omit() works for SILAC ratios because sum(na.omit(c(NA,NA,NA))) is 0.
  df_in$qc_missing_all_expression = apply(df_in, 1, function(x) {
    sum(na.omit(as.numeric(x[metadata_in$expression_column]))) == 0
  })
  
  #flag PEP score greater than threshold
  df_in$qc_pep = sapply(df_in$PEP, function(x) {
    if (x > 0.1) {
      TRUE
    } else {
      FALSE
    }
  })
  
  if (is_ptm_data == TRUE) {
    #flag missing modification count
    df_in$qc_mod_missing = sapply(df_in[, modification_count_column], function(x) {
      if (is.na(x)) {
        TRUE
      } else {
        FALSE
      }
    })
  }
  
  
  #get rid of rows that don't pass all qc
  if (is_ptm_data == TRUE) {
    qc_columns = c("qc_nothuman", "qc_reverse_only", "qc_missing_all_expression", "qc_pep", "qc_mod_missing")
    
  } else if (is_ptm_data == FALSE) {
    qc_columns = c("qc_nothuman", "qc_reverse_only", "qc_missing_all_expression", "qc_pep")
    
  }
  df_in$qc_pass = apply(df_in[, qc_columns], 1, function(x) { #weird things occurring with apply coersion; specify qc_columns in apply so that the booleans are read in as booleans and not characters
    sum(as.numeric(x[qc_columns]))
  })
  df_in_keep = subset(df_in, qc_pass == 0)
  
  #calculate number of PTMs for each site ----
  if (is_ptm_data == TRUE) {
    df_in_keep$PTM_position = ""
    for (i in 1:nrow(df_in_keep)) {
      df_in_keep$PTM_position[i] = get_ptm_position(df_in_keep[i, modification_prob_column], df_in_keep$Positions.within.proteins[i], df_in_keep$Position.in.peptide[i], ptm_threshold)
    }
  }
  
  #"average" expression by replicate but ignore missing values
  combined_expression = average_replicates(df_in_keep, metadata_in, missing_values)
  
  return(list(df_out_processed = cbind(df_in_keep, combined_expression), metadata_out_processed = metadata_in))
}

get_pep_from_peptides = function(peptide_ids_to_map, peptide_df) {
  pepscores = matrix(peptide_df[, "PEP"])
  row.names(pepscores) = peptide_df[, "id"]
  calculated_pep_score = vapply(peptide_ids_to_map, FUN.VALUE = numeric(1), function(x) {
    prod(pepscores[unlist(str_split(x, ";")),])
  })
  return(calculated_pep_score)
}


# sum_fractions = function(df_in, metadata_in) {
#     unique_samples = unique(metadata_in[, "sample"])
#     length_unique_sample_column = length(unique(metadata_in[, "sample"]))
#     summed_expression_matrix = matrix(data = NA, nrow = dim(df_in)[1], ncol = length_unique_sample_column)
#     for (i in 1:nrow(df_in)) {
#         for (j in 1:length_unique_sample_column) {
#             current_sample = unique_samples[j]
#             current_sample_fraction = unique(subset(metadata_in, sample == current_sample)$sample_fraction)
#             summed_expression_matrix[i,j] = sum(na.omit(as.numeric(df_in[i, current_sample_fraction])))
#         }
#     }
#     summed_expression_matrix = log2(summed_expression_matrix)
#     summed_expression_matrix[summed_expression_matrix == -Inf] = NA
#     colnames(summed_expression_matrix) = unique_samples
#     return(data.frame(summed_expression_matrix, stringsAsFactors = FALSE))
# }
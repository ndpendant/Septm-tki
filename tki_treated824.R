source("G:/Paul/SEPTM_tki/data_processing_functions.R")
cell_line_meta = read.delim("tki_treated_cell_meta.txt", header=TRUE, stringsAsFactors = F)

#data and IRON----

septm_data = c("phos", "acetyl", "ubi")
for (j in 1:length(septm_data)) {
  ptm_type = septm_data[j]
  
  if (ptm_type == "phos") {
    ptm_data = read.delim("Phospho (STY)Sites.txt", header = TRUE, stringsAsFactors = F)
  } else if (ptm_type == "acetyl") {
    ptm_data = read.delim("Acetyl (K)Sites.txt", header = TRUE, stringsAsFactors = F)
  } else if (ptm_type == "ubi") {
    ptm_data = read.delim("GlyGly (K)Sites.txt", header = TRUE, stringsAsFactors = F)
  } else {
    stop("No or bad ptm_type provided.")
  }
  
  intensity_column_names = names(ptm_data[, str_detect(names(ptm_data), "Intensity\\..*$")])
  intensity_meta = data.frame(
    stringsAsFactors = F, 
    column_name = intensity_column_names, 
    pass_missing = as.logical(apply(ptm_data[, intensity_column_names], 2, function(x) { sum(x != 0) > 0.01*nrow(ptm_data)})), 
    name = as.character(sapply(intensity_column_names, function(x) {str_match(x, paste(cell_line_meta$name, collapse="|"))[1, 1]})), 
    guolin_experiment = as.character(sapply(intensity_column_names, function(x) {str_match(x, "DMSO|Crizotinib|Afatinib|Erlotinib")[1, 1]})), 
    actual_experiment = ""
  )
  
  for (i in 1:nrow(intensity_meta)) {
    if (intensity_meta$guolin_experiment[i] != "DMSO") {
      intensity_meta$actual_experiment[i] = as.character(subset(cell_line_meta, name==intensity_meta$name[i])$treatment)
    } else {
      intensity_meta$actual_experiment[i] = "DMSO"
    }
  }
  
  if (ptm_type == "phos") {
    intensity_meta[intensity_meta$column_name == "Intensity.09152014_GZ_HCC827_Crizotinib.pY_run1", "actual_experiment"] = "DMSO"
    intensity_meta[intensity_meta$column_name == "Intensity.09152014_GZ_HCC827_Crizotinib.pY_run2", "actual_experiment"] = "DMSO"
    intensity_meta[intensity_meta$column_name == "Intensity.09152014_GZ_HCC827_DMSO.pY_run1", "actual_experiment"] = "Erlotinib"
    intensity_meta[intensity_meta$column_name == "Intensity.09152014_GZ_HCC827_DMSO.pY_run2", "actual_experiment"] = "Erlotinib"
  }
  #discard the failed/all 0 runs and normalize with IRON
  
  intensity_meta = subset(intensity_meta, intensity_meta$pass_missing == TRUE)
  intensity_column_names = intensity_meta$column_name
  #intensity_ironed = iron_proteomics(ptm_data[, intensity_meta$column_name])
  #intensity_meta$ironed_column_names = names(intensity_ironed)
  #ptm_data = cbind(ptm_data, intensity_ironed)
  
  #Sample QC----
  
  #correlations look reasonable, all > 0.5; majority > ~0.8
  for (i in 1:length(cell_line_meta$name)) {
    current_cell_line = cell_line_meta$name[i]
    cell_line_dmso = subset(intensity_meta, intensity_meta$name == current_cell_line & intensity_meta$actual_experiment == "DMSO")$column_name
    cell_line_treatment = subset(intensity_meta, intensity_meta$name == current_cell_line & intensity_meta$actual_experiment != "DMSO")$column_name
    dmso_cor_matrix = try(cor(ptm_data[, cell_line_dmso]))
    try( if (TRUE %in% (dmso_cor_matrix < 0.3)) {
      print("correlation in dmso_cor_matrix < 0.3")
      print(c("i=", i))
      print(current_cell_line)
      break
    })
    treatment_cor_matrix = try(cor(ptm_data[, cell_line_treatment]))
    try(if (TRUE %in% (treatment_cor_matrix < 0.3)) {
      print("correlation in treatment_cor_matrix < 0.3")
      print(c("i=", i))
      print(current_cell_line)
      break
    })
  }
  
  #Peptide QC----
  uniprot_entries = read.delim("uniprot_reviewed_human_refproteome_08-24-15.txt", header=T, stringsAsFactors = FALSE)
  ptm_data$qc_nothuman = sapply(ptm_data$Leading.Proteins, function(x) {
    sum(unlist(str_extract_all(x, uniprot_regex())) %in% uniprot_entries$Entry) == 0
  })
  
  #flag reverse-only hits
  ptm_data$qc_reverse_only = apply(ptm_data, 1, function(x) {
    str_detect(x["Leading.Proteins"], "^REV__[A-Z0-9]*$|^REV__[A-Z0-9]*;(REV__[A-Z0-9]*;)*REV__[A-Z0-9]*$")
  })
  
  #flag pep scores greater than 0.1
  ptm_data$qc_pep_gt_0.1 = apply(ptm_data, 1, function(x) {
    x["PEP"] > 0.1
  })
  
  #flag completely missing rows
  ptm_data$qc_missing_all_intensities = apply(ptm_data, 1, function(x) {
    sum(as.numeric(x[intensity_meta$column_name])) == 0
  })
  
  #flag if number of modifications is 0; need to use sapply to avoid apply coersion weirdness
  number_of_col_text = names(ptm_data)[str_detect(names(ptm_data), "Number\\.of\\..*")]
  ptm_data$qc_ptm_missing = sapply(ptm_data[, number_of_col_text], function(x) {
    as.numeric(x) == 0
  })
  
  ptm_probs_col_text = names(ptm_data)[str_detect(names(ptm_data), "Probabilities$")]
  ptm_data$qc_prob_missing = sapply(ptm_data[, ptm_probs_col_text], function(x) {
    x == ""
  })
  
  qc_columns = c("qc_nothuman", "qc_reverse_only", "qc_pep_gt_0.1", "qc_ptm_missing", "qc_missing_all_intensities", "qc_prob_missing")
  ptm_data$qc_pass = apply(ptm_data[, qc_columns], 1, function(x) { #weird things occurring with apply coersion; specify qc_columns in apply so that the booleans are read in as booleans and not characters
    sum(as.numeric(x[qc_columns]))
  })
  ptm_data_keep = subset(ptm_data, qc_pass == 0)
  
  if (ptm_type == "phos") {
    for (i in 1:nrow(ptm_data_keep)) {
      ptm_data_keep$PTM_position[i] = get_ptm_position(ptm_data_keep$Phospho..STY..Probabilities[i], ptm_data_keep$Positions[i], ptm_data_keep$Position.in.peptide[i], 0.5)
    }
  } else if (ptm_type == "acetyl") {
    for (i in 1:nrow(ptm_data_keep)) {
      ptm_data_keep$PTM_position[i] = get_ptm_position(ptm_data_keep$Acetyl..K..Probabilities[i], ptm_data_keep$Positions[i], ptm_data_keep$Position.in.peptide[i], 0.5)
    }
  } else if (ptm_type == "ubi") {
    for (i in 1:nrow(ptm_data_keep)) {
      ptm_data_keep$PTM_position[i] = get_ptm_position(ptm_data_keep$GlyGly..K..Probabilities[i], ptm_data_keep$Positions[i], ptm_data_keep$Position.in.peptide[i], 0.5)
    }
  }
  
  if (ptm_type == "phos") {
    ptm_data_keep[, "high_prob_Y_modifications"] = sapply(ptm_data_keep[, "PTM_position"], function(x) {
      sum(str_detect(string = unlist(str_split(x, ",")), pattern = "Y[0-9]*\\*"))
    })
    highprob_modification_col = "high_prob_Y_modifications"
  } else if (ptm_type == "acetyl") {
    ptm_data_keep[, "high_prob_K_modifications"] = sapply(ptm_data_keep[, "PTM_position"], function(x) {
      sum(str_detect(string = unlist(str_split(x, ",")), pattern = "K[0-9]*\\*"))
    })
    highprob_modification_col = "high_prob_K_modifications"
  } else if (ptm_type == "ubi") {
    ptm_data_keep[, "high_prob_K_modifications"] = sapply(ptm_data_keep[, "PTM_position"], function(x) {
      sum(str_detect(string = unlist(str_split(x, ",")), pattern = "K[0-9]*\\*"))
    })
    highprob_modification_col = "high_prob_K_modifications"
  }
  
  
  
  #Calcualte ratios between control and treatment pairs ----
  #going to use the "nonzero mean" for the calculation
  mean_ignore0 = function(vector_in) {
    vector_no0 = vector_in[vector_in != 0]
    return(mean(vector_no0))
  }
  
  #calculate ratio table by using mean_ignore0 function on all cell lines
  columns_to_include = c("Leading.Proteins", "Positions", "Gene.Names", "PTM_position", highprob_modification_col,"Modified.Sequence", "Protein.Names", "Localization.Prob", "Score.Diff", "PEP", "Score")
  ratio_table = ptm_data_keep[, c(columns_to_include)]
  for (i in 1:length(cell_line_meta$name)) {
    current_cell_line = cell_line_meta$name[i]
    dmso_colnames = subset(intensity_meta, name == current_cell_line & actual_experiment == "DMSO")$column_name
    treatment_colnames = subset(intensity_meta, name == current_cell_line & actual_experiment != "DMSO")$column_name
    dmso = ptm_data_keep[, dmso_colnames]
    treatment = ptm_data_keep[, treatment_colnames]
    #dmso or treatment may only have one dimension so use sapply, else apply
    if (is.vector(dmso) == TRUE) { 
      dmso_avg = sapply(dmso, mean_ignore0) #if one dimensional
    } else {
      dmso_avg = apply(dmso, 1, mean_ignore0)
    }
    if (is.vector(treatment) == TRUE) {
      treatment_avg = sapply(treatment, mean_ignore0) #if one dimensional
    } else {
      treatment_avg = apply(treatment, 1, mean_ignore0)
    }
    ratio_table[, paste(current_cell_line, "log2ratio", sep=".")] = log2(treatment_avg/dmso_avg)
    ratio_table[, paste(current_cell_line, "dmso", sep=".")] = log2(dmso_avg)
    ratio_table[, paste(current_cell_line, "treatment", sep=".")] = log2(treatment_avg)
  }
  dmso_names = names(ratio_table)[str_detect(names(ratio_table), "dmso")]
  treatment_names = names(ratio_table)[str_detect(names(ratio_table), "treatment")]
  ratio_names = names(ratio_table)[str_detect(names(ratio_table), "log2ratio")] #get the cell.line.log2ratio strings
  ratio_table$missing_ratios = apply(ratio_table[, ratio_names], 1, function(x) { sum(is.na(x))})
  ratio_table$overall_avg_log2_ratio = apply(ratio_table[, ratio_names], 1, function(x) { mean(na.omit(x))})
  
  cell_line_meta$ratio_names = ratio_names #add the cell.line.log2ratio strings to the cell line meta
  
  erlo_ratio_names = subset(cell_line_meta, treatment == "Erlotinib")$ratio_names
  crizo_ratio_names = subset(cell_line_meta, treatment == "Crizotinib")$ratio_names
  dasa_ratio_names = subset(cell_line_meta, treatment == "Dasatinib")$ratio_names
  afa_ratio_names = subset(cell_line_meta, treatment == "Afatinib")$ratio_names
  #summarize results for output
  #na breaks things so be sure to na.omit any applies
  ratio_table$crizotinib_avg_log2_ratio = apply(ratio_table[, crizotinib_ratio_names], 1, function(x) { mean(na.omit(x))})
  ratio_table$dasatinib_avg_log2_ratio = apply(ratio_table[, dasatinib_ratio_names], 1, function(x) { mean(na.omit(x))})
  ratio_table$erlotinib_avg_log2_ratio = apply(ratio_table[, erlotinib_ratio_names], 1, function(x) { mean(na.omit(x))})
  ratio_table = ratio_table[order(ratio_table$missing_ratios, ratio_table$overall_avg_log2_ratio), ]
  
      if (ptm_type == "phos") {
          write.table(ratio_table, "tki_treated_phospho_ratio_table_08-09-16.txt", sep="\t", row.names = FALSE)
      } else if (ptm_type == "acetyl") {
          write.table(ratio_table, "tki_treated_acetyl_ratio_table_08-09-16.txt", sep="\t", row.names = FALSE)
      } else if (ptm_type == "ubi") {
          write.table(ratio_table, "tki_treated_ubi_ratio_table_08-09-16.txt", sep="\t", row.names = FALSE)
      }
}
test1=cbind(ratio_table_phos$H4006.log2ratio,ratio_table_phos$HCC827.log2ratio,ratio_table_phos$PC9.log2ratio, checkmeout$incc2, ratio_table_phos$'#increased_2fold_erlo',checkmeout$decc2, ratio_table_phos$`#decreased_2fold_erlo` ,checkmeout$incc1.5, ratio_table_phos$`#increased_1.5fold_erlo`,
            checkmeout$decc1.5,ratio_table_phos$`#decreased_1.5fold_erlo` )
checkmeout<-genaddform(ratio_table_phos,erlo_ratio_names,c("incc1.5","decc1.5","incc2","decc2"),c('1.5i','1.5d','2i','2d'))


ratio_table_phos=genaddform(ratio_table_phos,erlo_ratio_names,c("#increased_1.5fold_erlo","#decreased_1.5fold_erlo","#increased_2fold_erlo","#decreased_2fold_erlo"),c('1.5i','1.5d','2i','2d'))
ratio_table_phos=genaddform(ratio_table_phos,crizo_ratio_names,c("#increased_1.5fold_crizo","#decreased_1.5fold_crizo","#increased_2fold_crizo","#decreased_2fold_crizo"),c('1.5i','1.5d','2i','2d'))
ratio_table_phos=genaddform(ratio_table_phos,dasa_ratio_names,c("#increased_1.5fold_dasa","#decreased_1.5fold_dasa","#increased_2fold_dasa","#decreased_2fold_dasa"),c('1.5i','1.5d','2i','2d'))
ratio_table_phos=genaddform(ratio_table_phos,afa_ratio_names,c("#increased_1.5fold_afa","#decreased_1.5fold_afa","#increased_2fold_afa","#decreased_2fold_afa"),c('1.5i','1.5d','2i','2d'))

ratio_table_ubi=genaddform(ratio_table_ubi,erlo_ratio_names,c("#increased_1.5fold_erlo","#decreased_1.5fold_erlo","#increased_2fold_erlo","#decreased_2fold_erlo"),c('1.5i','1.5d','2i','2d'))
ratio_table_ubi=genaddform(ratio_table_ubi,crizo_ratio_names,c("#increased_1.5fold_crizo","#decreased_1.5fold_crizo","#increased_2fold_crizo","#decreased_2fold_crizo"),c('1.5i','1.5d','2i','2d'))
ratio_table_ubi=genaddform(ratio_table_ubi,dasa_ratio_names,c("#increased_1.5fold_dasa","#decreased_1.5fold_dasa","#increased_2fold_dasa","#decreased_2fold_dasa"),c('1.5i','1.5d','2i','2d'))
ratio_table_ubi=genaddform(ratio_table_ubi,afa_ratio_names,c("#increased_1.5fold_afa","#decreased_1.5fold_afa","#increased_2fold_afa","#decreased_2fold_afa"),c('1.5i','1.5d','2i','2d'))

ratio_table_acetyl=genaddform(ratio_table_acetyl,erlo_ratio_names,c("#increased_1.5fold_erlo","#decreased_1.5fold_erlo","#increased_2fold_erlo","#decreased_2fold_erlo"),c('1.5i','1.5d','2i','2d'))
ratio_table_acetyl=genaddform(ratio_table_acetyl,crizo_ratio_names,c("#increased_1.5fold_crizo","#decreased_1.5fold_crizo","#increased_2fold_crizo","#decreased_2fold_crizo"),c('1.5i','1.5d','2i','2d'))
ratio_table_acetyl=genaddform(ratio_table_acetyl,dasa_ratio_names,c("#increased_1.5fold_dasa","#decreased_1.5fold_dasa","#increased_2fold_dasa","#decreased_2fold_dasa"),c('1.5i','1.5d','2i','2d'))
ratio_table_acetyl=genaddform(ratio_table_acetyl,afa_ratio_names,c("#increased_1.5fold_afa","#decreased_1.5fold_afa","#increased_2fold_afa","#decreased_2fold_afa"),c('1.5i','1.5d','2i','2d'))



# 
# 
# 
# #Ratio Table analysis----
ratio_table_phos = read.delim("tki_treated_phospho_ratio_table_08-09-16.txt", header = TRUE, stringsAsFactors = F)
ratio_table_acetyl = read.delim("tki_treated_acetyl_ratio_table_08-09-16.txt", header = TRUE, stringsAsFactors = F)
ratio_table_ubi = read.delim("tki_treated_ubi_ratio_table_08-09-16.txt", header = TRUE, stringsAsFactors = F)
uniprot_entries = read.delim("uniprot_reviewed_human_refproteome_08-24-15.txt", header=T, stringsAsFactors = FALSE)


#Testing for accuracy---------
checkmeout<-genaddform(ratio_table_phos,erlo_ratio_names,c("incc1.5","decc1.5","incc2","decc2","ratios"),c('1.5i','1.5d','2i','2d','!na.is'))

test1=cbind(ratio_table_phos$H4006.log2ratio,ratio_table_phos$HCC827.log2ratio,ratio_table_phos$PC9.log2ratio, checkmeout$incc2, ratio_table_phos$'#increased_2fold_erlo',checkmeout$decc2, ratio_table_phos$`#decreased_2fold_erlo` ,checkmeout$incc1.5, ratio_table_phos$`#increased_1.5fold_erlo`,
           checkmeout$decc1.5,ratio_table_phos$`#decreased_1.5fold_erlo`, checkmeout$ratios , ratio_table_phos$`#erlo_ratios`)

#--------------------------

ratio_table_phos=genaddform(ratio_table_phos,erlo_ratio_names,c("#increased_1.5fold_erlo","#decreased_1.5fold_erlo","#increased_2fold_erlo","#decreased_2fold_erlo", "#erlo_ratios"),c('1.5i','1.5d','2i','2d','!na.is'))
ratio_table_phos=genaddform(ratio_table_phos,crizo_ratio_names,c("#increased_1.5fold_crizo","#decreased_1.5fold_crizo","#increased_2fold_crizo","#decreased_2fold_crizo", "#crizo_ratios"),c('1.5i','1.5d','2i','2d','!na.is'))
ratio_table_phos=genaddform(ratio_table_phos,dasa_ratio_names,c("#increased_1.5fold_dasa","#decreased_1.5fold_dasa","#increased_2fold_dasa","#decreased_2fold_dasa", "#dasa_ratios"),c('1.5i','1.5d','2i','2d','!na.is'))
ratio_table_phos=genaddform(ratio_table_phos,afa_ratio_names,c("#increased_1.5fold_afa","#decreased_1.5fold_afa","#increased_2fold_afa","#decreased_2fold_afa", "#afa_ratios"),c('1.5i','1.5d','2i','2d','!na.is'))

ratio_table_ubi=genaddform(ratio_table_ubi,erlo_ratio_names,c("#increased_1.5fold_erlo","#decreased_1.5fold_erlo","#increased_2fold_erlo","#decreased_2fold_erlo", "#erlo_ratios"),c('1.5i','1.5d','2i','2d','!na.is'))
ratio_table_ubi=genaddform(ratio_table_ubi,crizo_ratio_names,c("#increased_1.5fold_crizo","#decreased_1.5fold_crizo","#increased_2fold_crizo","#decreased_2fold_crizo", "#crizo_ratios"),c('1.5i','1.5d','2i','2d','!na.is'))
ratio_table_ubi=genaddform(ratio_table_ubi,dasa_ratio_names,c("#increased_1.5fold_dasa","#decreased_1.5fold_dasa","#increased_2fold_dasa","#decreased_2fold_dasa", "#dasa_ratios"),c('1.5i','1.5d','2i','2d','!na.is'))
ratio_table_ubi=genaddform(ratio_table_ubi,afa_ratio_names,c("#increased_1.5fold_afa","#decreased_1.5fold_afa","#increased_2fold_afa","#decreased_2fold_afa", "#afa_ratios"),c('1.5i','1.5d','2i','2d','!na.is'))

ratio_table_acetyl=genaddform(ratio_table_acetyl,erlo_ratio_names,c("#increased_1.5fold_erlo","#decreased_1.5fold_erlo","#increased_2fold_erlo","#decreased_2fold_erlo", "#dasa_ratios"),c('1.5i','1.5d','2i','2d','!na.is'))
ratio_table_acetyl=genaddform(ratio_table_acetyl,crizo_ratio_names,c("#increased_1.5fold_crizo","#decreased_1.5fold_crizo","#increased_2fold_crizo","#decreased_2fold_crizo", "#crizo_ratios"),c('1.5i','1.5d','2i','2d','!na.is'))
ratio_table_acetyl=genaddform(ratio_table_acetyl,dasa_ratio_names,c("#increased_1.5fold_dasa","#decreased_1.5fold_dasa","#increased_2fold_dasa","#decreased_2fold_dasa", "#dasa_ratios"),c('1.5i','1.5d','2i','2d','!na.is'))
ratio_table_acetyl=genaddform(ratio_table_acetyl,afa_ratio_names,c("#increased_1.5fold_afa","#decreased_1.5fold_afa","#increased_2fold_afa","#decreased_2fold_afa", "#afa_ratios"),c('1.5i','1.5d','2i','2d','!na.is'))



# ratio_table_acetyl[, "#increased_1.5fold_erlo"] = apply(ratio_table_acetyl[, erlo_ratio_names], 1, function(x){sum(na.omit(x) > log2(1.5))})
# ratio_table_acetyl[, "#decreased_1.5fold_erlo"] = apply(ratio_table_acetyl[, erlo_ratio_names], 1, function(x){sum(na.omit(x) < log2(0.667))})
# ratio_table_phos[, "#increased_1.5fold_erlo"] = apply(ratio_table_phos[, erlo_ratio_names], 1, function(x){sum(na.omit(x) > log2(1.5))})
# ratio_table_phos[, "#decreased_1.5fold_erlo"] = apply(ratio_table_phos[, erlo_ratio_names], 1, function(x){sum(na.omit(x) < log2(0.667))})
# ratio_table_ubi[, "#increased_1.5fold_erlo"] = apply(ratio_table_ubi[, erlo_ratio_names], 1, function(x){sum(na.omit(x) > log2(1.5))})
# ratio_table_ubi[, "#decreased_1.5fold_erlo"] = apply(ratio_table_ubi[, erlo_ratio_names], 1, function(x){sum(na.omit(x) < log2(0.667))})
# 
# ratio_table_acetyl[, "#increased_1.5fold_crizo"] = apply(ratio_table_acetyl[, crizo_ratio_names], 1, function(x){sum(na.omit(x) > log2(1.5))})
# ratio_table_acetyl[, "#decreased_1.5fold_crizo"] = apply(ratio_table_acetyl[, crizo_ratio_names], 1, function(x){sum(na.omit(x) < log2(0.667))})
# ratio_table_phos[, "#increased_1.5fold_crizo"] = apply(ratio_table_phos[, crizo_ratio_names], 1, function(x){sum(na.omit(x) > log2(1.5))})
# ratio_table_phos[, "#decreased_1.5fold_crizo"] = apply(ratio_table_phos[, crizo_ratio_names], 1, function(x){sum(na.omit(x) < log2(0.667))})
# ratio_table_ubi[, "#increased_1.5fold_crizo"] = apply(ratio_table_ubi[, crizo_ratio_names], 1, function(x){sum(na.omit(x) > log2(1.5))})
# ratio_table_ubi[, "#decreased_1.5fold_crizo"] = apply(ratio_table_ubi[, crizo_ratio_names], 1, function(x){sum(na.omit(x) < log2(0.667))})
# 
# ratio_table_acetyl[, "#increased_1.5fold_dasa"] = apply(ratio_table_acetyl[, dasa_ratio_names], 1, function(x){sum(na.omit(x) > log2(1.5))})
# ratio_table_acetyl[, "#decreased_1.5fold_dasa"] = apply(ratio_table_acetyl[, dasa_ratio_names], 1, function(x){sum(na.omit(x) < log2(0.667))})
# ratio_table_phos[, "#increased_1.5fold_dasa"] = apply(ratio_table_phos[, dasa_ratio_names], 1, function(x){sum(na.omit(x) > log2(1.5))})
# ratio_table_phos[, "#decreased_1.5fold_dasa"] = apply(ratio_table_phos[, dasa_ratio_names], 1, function(x){sum(na.omit(x) < log2(0.667))})
# ratio_table_ubi[, "#increased_1.5fold_dasa"] = apply(ratio_table_ubi[, dasa_ratio_names], 1, function(x){sum(na.omit(x) > log2(1.5))})
# ratio_table_ubi[, "#decreased_1.5fold_dasa"] = apply(ratio_table_ubi[, dasa_ratio_names], 1, function(x){sum(na.omit(x) < log2(0.667))})
# 
# ratio_table_acetyl[, "#increased_1.5fold_afa"] = sapply(ratio_table_acetyl[, afa_ratio_names], function(x){sum(na.omit(x) > log2(1.5))})
# ratio_table_acetyl[, "#decreased_1.5fold_afa"] = sapply(ratio_table_acetyl[, afa_ratio_names], function(x){sum(na.omit(x) < log2(0.667))})
# ratio_table_phos[, "#increased_1.5fold_afa"] = sapply(ratio_table_phos[, afa_ratio_names], function(x){sum(na.omit(x) > log2(1.5))})
# ratio_table_phos[, "#decreased_1.5fold_afa"] = sapply(ratio_table_phos[, afa_ratio_names], function(x){sum(na.omit(x) < log2(0.667))})
# ratio_table_ubi[, "#increased_1.5fold_afa"] = sapply(ratio_table_ubi[, afa_ratio_names], function(x){sum(na.omit(x) > log2(1.5))})
# ratio_table_ubi[, "#decreased_1.5fold_afa"] = sapply(ratio_table_ubi[, afa_ratio_names], function(x){sum(na.omit(x) < log2(0.667))})
# 
# 
# ratio_table_acetyl[, "#increased_2fold_erlo"] = apply(ratio_table_acetyl[, erlo_ratio_names], 1, function(x){sum(na.omit(x) > 1)})
# ratio_table_acetyl[, "#decreased_2fold_erlo"] = apply(ratio_table_acetyl[, erlo_ratio_names], 1, function(x){sum(na.omit(x) < -1)})
# ratio_table_phos[, "#increased_2fold_erlo"] = apply(ratio_table_phos[, erlo_ratio_names], 1, function(x){sum(na.omit(x) > 1)})
# ratio_table_phos[, "#decreased_2fold_erlo"] = apply(ratio_table_phos[, erlo_ratio_names], 1, function(x){sum(na.omit(x) < -1)})
# ratio_table_ubi[, "#increased_2fold_erlo"] = apply(ratio_table_ubi[, erlo_ratio_names], 1, function(x){sum(na.omit(x) > 1)})
# ratio_table_ubi[, "#decreased_2fold_erlo"] = apply(ratio_table_ubi[, erlo_ratio_names], 1, function(x){sum(na.omit(x) < -1)})
# 
# ratio_table_acetyl[, "#increased_2fold_crizo"] = apply(ratio_table_acetyl[, crizo_ratio_names], 1, function(x){sum(na.omit(x) > 1)})
# ratio_table_acetyl[, "#decreased_2fold_crizo"] = apply(ratio_table_acetyl[, crizo_ratio_names], 1, function(x){sum(na.omit(x) < -1)})
# ratio_table_phos[, "#increased_2fold_crizo"] = apply(ratio_table_phos[, crizo_ratio_names], 1, function(x){sum(na.omit(x) > 1)})
# ratio_table_phos[, "#decreased_2fold_crizo"] = apply(ratio_table_phos[, crizo_ratio_names], 1, function(x){sum(na.omit(x) < -1)})
# ratio_table_ubi[, "#increased_2fold_crizo"] = apply(ratio_table_ubi[, crizo_ratio_names], 1, function(x){sum(na.omit(x) > 1)})
# ratio_table_ubi[, "#decreased_2fold_crizo"] = apply(ratio_table_ubi[, crizo_ratio_names], 1, function(x){sum(na.omit(x) < -1)})
# 
# ratio_table_acetyl[, "#increased_2fold_dasa"] = apply(ratio_table_acetyl[, dasa_ratio_names], 1, function(x){sum(na.omit(x) > 1)})
# ratio_table_acetyl[, "#decreased_2fold_dasa"] = apply(ratio_table_acetyl[, dasa_ratio_names], 1, function(x){sum(na.omit(x) < -1)})
# ratio_table_phos[, "#increased_2fold_dasa"] = apply(ratio_table_phos[, dasa_ratio_names], 1, function(x){sum(na.omit(x) > 1)})
# ratio_table_phos[, "#decreased_2fold_dasa"] = apply(ratio_table_phos[, dasa_ratio_names], 1, function(x){sum(na.omit(x) < -1)})
# ratio_table_ubi[, "#increased_2fold_dasa"] = apply(ratio_table_ubi[, dasa_ratio_names], 1, function(x){sum(na.omit(x) > 1)})
# ratio_table_ubi[, "#decreased_2fold_dasa"] = apply(ratio_table_ubi[, dasa_ratio_names], 1, function(x){sum(na.omit(x) < -1)})
# 
# ratio_table_acetyl[, "#increased_2fold_afa"] = sapply(ratio_table_acetyl[, afa_ratio_names], function(x){sum(na.omit(x) > 1)})
# ratio_table_acetyl[, "#decreased_2fold_afa"] = sapply(ratio_table_acetyl[, afa_ratio_names], function(x){sum(na.omit(x) < -1)})
# ratio_table_phos[, "#increased_2fold_afa"] = sapply(ratio_table_phos[, afa_ratio_names], function(x){sum(na.omit(x) > 1)})
# ratio_table_phos[, "#decreased_2fold_afa"] = sapply(ratio_table_phos[, afa_ratio_names], function(x){sum(na.omit(x) < -1)})
# ratio_table_ubi[, "#increased_2fold_afa"] = sapply(ratio_table_ubi[, afa_ratio_names], function(x){sum(na.omit(x) > 1)})
# ratio_table_ubi[, "#decreased_2fold_afa"] = sapply(ratio_table_ubi[, afa_ratio_names], function(x){sum(na.omit(x) < -1)})
# 



ratio_table_acetyl[, "#erlo_ratios"] = apply(ratio_table_acetyl[, erlo_ratio_names], 1, function(x){sum(!is.na(x))})
ratio_table_acetyl[, "#crizo_ratios"] = apply(ratio_table_acetyl[, crizo_ratio_names], 1, function(x){sum(!is.na(x))})
ratio_table_acetyl[, "#dasa_ratios"] = apply(ratio_table_acetyl[, dasa_ratio_names], 1, function(x){sum(!is.na(x))})
ratio_table_acetyl[, "#afa_ratios"] = sapply(ratio_table_acetyl[, afa_ratio_names], function(x){sum(!is.na(x))})

ratio_table_phos[, "#erlo_ratios"] = apply(ratio_table_phos[, erlo_ratio_names], 1, function(x){sum(!is.na(x))})
ratio_table_phos[, "#crizo_ratios"] = apply(ratio_table_phos[, crizo_ratio_names], 1, function(x){sum(!is.na(x))})
ratio_table_phos[, "#dasa_ratios"] = apply(ratio_table_phos[, dasa_ratio_names], 1, function(x){sum(!is.na(x))})
ratio_table_phos[, "#afa_ratios"] = sapply(ratio_table_phos[, afa_ratio_names], function(x){sum(!is.na(x))})

ratio_table_ubi[, "#erlo_ratios"] = apply(ratio_table_ubi[, erlo_ratio_names], 1, function(x){sum(!is.na(x))})
ratio_table_ubi[, "#crizo_ratios"] = apply(ratio_table_ubi[, crizo_ratio_names], 1, function(x){sum(!is.na(x))})
ratio_table_ubi[, "#dasa_ratios"] = apply(ratio_table_ubi[, dasa_ratio_names], 1, function(x){sum(!is.na(x))})
ratio_table_ubi[, "#afa_ratios"] = sapply(ratio_table_ubi[, afa_ratio_names], function(x){sum(!is.na(x))})


# write.table(ratio_table_acetyl, "output/ratio_table_acetyl_08-09-16.txt", sep="\t", row.names=F)
# write.table(ratio_table_phos, "output/ratio_table_phos_08-09-16.txt", sep="\t", row.names=F)
# write.table(ratio_table_ubi, "output/ratio_table_ubi_08-09-16.txt", sep="\t", row.names=F)

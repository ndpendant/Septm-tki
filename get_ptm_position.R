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

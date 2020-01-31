# LABEL PARSING
label_parser = function(data_dir) {
  
  labels = list()
  
  csv_file = 'stage_2_train.csv'
  labels[['train']] = parse.csv(data_dir, csv_file, data_split = 'train')
  
  csv_file = 'stage_2_sample_submission.csv'
  labels[['test']] = parse.csv(data_dir, csv_file, data_split = 'test')
  
  return(labels)
  
}

parse.csv = function(data_dir, csv_file, data_split = NA) {
  
  # The training data is provided as a set of image Ids and multiple labels, 
  # one for each of five sub-types of hemorrhage, plus an additional label for any, 
  # which should always be true if any of the sub-type labels is true.
  # 
  # There is also a target column, Label, indicating the probability of whether 
  # that type of hemorrhage exists in the indicated image.
  # 
  # There will be 6 rows per image Id. The label indicated by a particular row will look like 
  # [Image Id]_[Sub-type Name], as follows:
  #   
  # Id,Label
  # 1_epidural_hemorrhage,0
  # 1_intraparenchymal_hemorrhage,0
  # 1_intraventricular_hemorrhage,0
  # 1_subarachnoid_hemorrhage,0.6
  # 1_subdural_hemorrhage,0
  # 1_any,0.9
  
  filepath_to_read = file.path(data_dir, csv_file)
  cat(paste0('Reading in the ', data_split, ' split from "', csv_file, '"\n'))
  raw_csv = read.csv(filepath_to_read, stringsAsFactors = FALSE)
  cat(paste0(' ... number of rows = ', length(raw_csv[['ID']]), '\n'))
  cat(paste0(' ... number of files = ', length(raw_csv[['ID']])/6, '\n'))
  
  cat(paste0('       --> Make each label its own column\n'))
  labels_out = labels.rows.to.cols(raw_csv)
  labels_out[['counts']] = get.label.stats(labels_out[['labels']])
  
  return(labels_out)
  
}


labels.rows.to.cols = function(raw_csv) {
  
  no_of_rows =  length(raw_csv[['ID']])
  no_of_dcm_files = length(raw_csv[['ID']])/6
  
  col_names = get.colnames(one_dcm_file = raw_csv$ID[1:6])
  idxs = create.idxs(col_names, no_of_dcm_files)
  labels_as_cols = to.cols.from.idxs(idxs, raw_csv)
  
  return(labels_as_cols)
  
}

get.colnames = function(one_dcm_file) {
  
  ID = sapply(one_dcm_file, function(x) strsplit(x, "_")[[1]][2], USE.NAMES=FALSE)[1]
  col_names = sapply(one_dcm_file, function(x) tail(strsplit(x, "_")[[1]],1), USE.NAMES=FALSE)
  
  return(col_names)
  
}

create.idxs = function(col_names, no_of_dcm_files, no_of_labels = 6) {
  
  labels_to_rep = list()
  labels_empty = rep(FALSE, no_of_labels)
  labels_idxs = list()
  
  for (l in 1 : no_of_labels) {
    
    # set correct one as TRUE
    labels_to_rep[[col_names[l]]] = labels_empty
    labels_to_rep[[col_names[l]]][l] = TRUE
    
    # replicate
    labels_idxs[[col_names[l]]] = rep(labels_to_rep[[col_names[l]]], no_of_dcm_files)
    
  }
  
  return(labels_idxs)
  
}

to.cols.from.idxs = function(idxs, raw_csv) {
  
  col_names_in = names(idxs)
  labels_as_cols = list()
  
  for (n in 1 : length(col_names_in)) {
    colname = col_names_in[n]
    labels_as_cols[[colname]] = raw_csv[['Label']][idxs[[colname]]]
    
    if (n == 1) {
      raw_ids = raw_csv[['ID']][idxs[[colname]]]
      ids = sapply(raw_ids, function(x) strsplit(x, "_")[[1]][2], USE.NAMES=FALSE)
    }
  }
  
  # Add a label for control (healthy) slice from ANY, TODO!
  labels_as_cols[['healthy']] = 1 - labels_as_cols[['any']]
  
  labels_out = list()
  labels_out[['labels']] = labels_as_cols
  labels_out[['ID']] = ids
  
  return(labels_out)
  
}

get.label.stats = function(labels_as_cols, label_threshold = 0.5) {
  
  label_counts = list()
  label_means = list()
  col_names_in = names(labels_as_cols)
  
  for (n in 1 : length(col_names_in)) {
    colname = col_names_in[n]
    valid_rows = labels_as_cols[[colname]] > label_threshold
    non_zeros = labels_as_cols[[colname]][valid_rows]
    label_counts[[colname]] = length(non_zeros)
    # label_means[[colname]] = mean(non_zeros) # apparently everything is boolean, not as outlined in Kaggle
  }
  
  return(label_counts)
}

get.ids = function(labels_as_cols, raw_csv) {
  
  ids = list()
  col_names_in = names(labels_as_cols)
  
  for (n in 1 : length(col_names_in)) {
    colname = col_names_in[n]
    raw_ids = raw_csv[['ID']]
  }
  
}
call.hist.check() {
  path_air = file.path(out_dir_NII, '..', 'Kaggle_air_histograms')
  check.for.histogram.stats(path_air)
}



check.for.histogram.stats = function(path_air) {
  
  files_hist = list.files(path = path_air, pattern = '*.RData')
  airMode = as.integer(sapply(files_hist, function(x) strsplit(x, "_")[[1]][1], USE.NAMES=FALSE))
  non_zero = !airMode == 0
  
  files_to_check = files_hist[non_zero]
  codes_to_find = sapply(files_to_check, function(x) strsplit(x, "_")[[1]][4], USE.NAMES=FALSE)
  airModes = airMode[non_zero]
 
  df_out = data.frame(code = codes_to_find,
                      HU_mode = airModes) 
  
  write.csv(df_out, file = file.path(path_air, '..', 'nonZero_airs.csv'), row.names = FALSE)
  
}
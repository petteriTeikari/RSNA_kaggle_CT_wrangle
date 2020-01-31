call.hist.check() {
  path_headers = '/home/petteri/Dropbox/manuscriptDrafts/CThemorr/HOME_DATA/RSNA_Kaggle_Headers'
  check.for.header.quality(path_headers)
}



check.for.header.quality = function(path_headers) {
  
  files = list.files(path_headers, pattern = '*.RData', full.names = TRUE)
  fnames = list.files(path_headers, pattern = '*.RData', full.names = FALSE)
  file_idxs = seq(from = 1, to = length(files))
  
  imagePos_array = matrix(NA, nrow = length(files), ncol = 3)
  imageOrient_array = matrix(NA, nrow = length(files), ncol = 6)
  ID = matrix(NA, nrow = length(files), ncol = 1) # unique for each DCM
  study_ID = matrix(NA, nrow = length(files), ncol = 1) # unique for each VOLUME
  
  # Populate the preallocated matrices
  # TODO! this could be faster 
  for (f in 1 : length(files)) {
    
    if ((f %% 1000) == 0) { cat('f = ', f, '\n') }
    
    load(files[f]) 
    
    idx = grep('ImageOrientationPatient', hdr$name)
    vector = as.numeric(strsplit(hdr$value[idx], ' ')[[1]])
    imageOrient_array[f,] = vector
    
    idx = grep('ImagePositionPatient', hdr$name)
    vector = as.numeric(strsplit(hdr$value[idx], ' ')[[1]])
    imagePos_array[f,] = vector
    
    ID[f] = gsub('.RData', '', fnames[f])
    
    idx = grep('StudyInstanceUID', hdr$name)
    study_ID[f] = hdr$value[idx]
    
  }
  
  # Combine into a data frame to be written to disk
  data_out = cbind(study_ID, ID, imagePos_array, imageOrient_array)
  no_of_NAs = sum(is.na(data_out)) # = 0
  
  df = data.frame(data_out, stringsAsFactors = FALSE)
  colnames(df) = c('StudyID', 'ID', 'ImagePos_x', 'ImagePos_y', 'ImagePos_z', 
                   'ImageOrient_1', 'ImageOrient_2', 'ImageOrient_3',
                   'ImageOrient_4', 'ImageOrient_5', 'ImageOrient_6')
  
  # KAGGLE_imagePosAndOrientation_metadata.tar.xz
  write.csv(df, file = file.path(path_headers, '..', 'KAGGLE_imagePosAndOrientation_metadata.csv'))
  
}
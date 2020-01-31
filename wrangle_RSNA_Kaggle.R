call_RSNA_Kaggle_wrangling = function() {
  
  source('~/Dropbox/manuscriptDrafts/CThemorr/CT_hemorrhage/R_analysis/RSNA_Kaggle_label_parsing.R')
  if (!require("parallel")) install.packages("parallel"); library("parallel")
  if (!require("imager")) install.packages("imager"); library("imager")
  if (!require("dcm2niir")) install.packages("dcm2niir"); library("dcm2niir")
  # install_dcm2nii()
  if (!require("oro.dicom")) install.packages("oro.dicom"); library("oro.dicom")
  if (!require("RNifti")) install.packages("RNifti"); library("RNifti")
  if (!require("PET")) install.packages("PET"); library("PET")
  # sudo apt install libgsl-dev
  # if (!require("fmri")) install.packages("fmri"); library("fmri")
  # if (!require("dcemri")) install.packages("dcemri", repos="http://R-Forge.R-project.org"); library("dcemri")
  
  # https://www.kaggle.com/c/rsna-intracranial-hemorrhage-detection/data
  data_dir = '/mnt/wwn-0x5000c500b6b1dabd/rsna-intracranial-hemorrhage-detection'
  wrangle_RSNA_Kaggle(data_dir)  
  
}

wrangle_RSNA_Kaggle = function(data_dir) {

  labels = label_parser(data_dir)
  listing = get.dcm.listing(data_dir, labels, labels_to_use = 'train')
  create.volumes.from.dcm.slices(data_dir, listing)
  
  # Copy all the intraparencymal slices somewhere with Dicom -> Nifti conversion
  wildcard = 'intraparenchymal'
  copy.wildcard.scans.as.NII.to.another.folder(data_dir, labels,
                                               raw_csv, wildcard)
  
}  


get.dcm.listing = function(data_dir, labels, labels_to_use = 'train') {
  
  if (labels_to_use == 'train') {
    subdir = 'stage_2_train'
  } else {
    cat('Well we do not know the labels of the test set\n')
  }
  
  cat(paste0('Getting file listing from "', subdir, '"\n'))
  files = list.files(path = file.path(data_dir, subdir), pattern = '*.dcm')
  filepaths = list.files(path = file.path(data_dir, subdir), pattern = '*.dcm', full.names = TRUE)
  IDs_from_files = sapply(files, function(x) strsplit(gsub('.dcm', '', x), "_")[[1]][2], USE.NAMES=FALSE)
  cat(paste0(' ... number of files found = ', length(files), '\n'))
  
  cat(paste0('Getting metadata as well from the .dcms (takes a minute)\n'))
  hdr_path_out = '/home/petteri/Dropbox/manuscriptDrafts/CThemorr/HOME_DATA/RSNA_Kaggle_Headers'
  dir.create(hdr_path_out, showWarnings = FALSE)
  meta_out = get.metadata.from.listing(files, subdir, data_dir, IDs_from_files, hdr_path_out)
  
  reconstruct.volumes.from.slices(meta_out, files, filepaths, hdr_path_out, labels)
  
  listing = list()
  listing[['files']] = files
  listing[['IDs']] = IDs_from_files
  listing[['meta_out']] = meta_out
  
  return(listing)
}

reconstruct.volumes.from.slices = function(meta_out, files, filepaths, hdr_path_out, labels,
                                           no_cores_to_use = 8, check_done = FALSE,
                                           check_weird_airs = TRUE) {
  
  tmp_dcm_dir = file.path('/home/petteri/', 'temp_Kaggle')
  # dir.create(tmp_dcm_dir, showWarnings = FALSE) # useful if the dcm2niix would work
  
  out_dir_NII = file.path('/home/petteri/', 'Kaggle_NII')
  out_dir_NII = file.path('/mnt/wwn-0x5000c500b6b1dabd', 'Kaggle_NII')
  dir.create(out_dir_NII, showWarnings = FALSE)
  
  out_dir_NII_labels = file.path('/home/petteri/Dropbox/manuscriptDrafts/CThemorr/HOME_DATA', 'Kaggle_NII_labels')
  dir.create(out_dir_NII_labels, showWarnings = FALSE)
  
  out_dir_NII_mask = file.path('/mnt/wwn-0x5000c500b6b1dabd', 'Kaggle_NII_imageAreaLabels')
  dir.create(out_dir_NII_mask, showWarnings = FALSE)
  
  hdr_filepaths = list.files(hdr_path_out, pattern = '.*.RData', full.names = TRUE)
  
  # http://gdcm.sourceforge.net/wiki/index.php/Imager_Pixel_Spacing#Ordering_of_slice_to_reconstruct_a_volume
  # https://www.kaggle.com/anjum48/reconstructing-3d-volumes-from-metadata/data
  # Reconstructing 3D volumes from metadata
  # ---------------------------------------
  # "In this kernel we will use the StudyInstanceUID in the metadata to group 
  # together images from the same scan. We will then sort the images using 
  # ImagePositionPatient_2 and create 3D volumes. 
  study_IDs = meta_out[['metadata']][['StudyInstanceUID']]
  unique_studies = unique(study_IDs)
  
  if (check_done) {
    unique_studies = check.for.done.studies(unique_studies, done_dir = out_dir_NII)
  }
  
  if (check_weird_airs) {
    
    checked = recheck.for.weird.air.baselines(unique_studies,
                                                     done_dir = out_dir_NII,
                                                     mask_dir = out_dir_NII_mask,
                                                     LUT_file = '/mnt/wwn-0x5000c500b6b1dabd/nonZero_airs.csv',
                                                     correction_file = '/home/petteri/Dropbox/manuscriptDrafts/CThemorr/HOME_DATA/kaggle_CT_HUs.csv')
    unique_studies = checked$unique_studies
    air_values = checked$corrs$AIR
    
  } else {
    air_values = NA
  }
  
  
  if_pars = rep(NA, length(unique_studies))
  for (st in 1 :length(unique_studies)) {
    
    if_pars[st] = try({
      CT.study.process.wrapper(st, unique_studies, study_IDs, meta_out, 
                               files, filepaths, hdr_filepaths, labels, air_values, check_weird_airs,
                               tmp_dcm_dir, out_dir_NII, out_dir_NII_mask, out_dir_NII_labels)
    }, silent = TRUE)
    
    if (class(out) == "try-error") {
      cat('  error with reading = ', files[f], ', \n') # "ID_056d4b375c" / "error with reading =  ID_00032d440.dcm"
      cat('  ', out) 
      # Error in CT_array[, , s] <- DCM_slices[[s]] : 
      # number of items to replace is not a multiple of replacement length
    }
    
  }
  
  st = seq(from = 1, to = length(unique_studies))
  # st = seq(from = 1, to = 12)
  if_pars = mclapply(st, function(f){
                     CT.study.process.wrapper(st, unique_studies, study_IDs, meta_out, 
                                              files, filepaths, hdr_filepaths, labels, 
                                              tmp_dcm_dir, out_dir_NII, out_dir_NII_mask, out_dir_NII_labels)}, 
                     mc.cores = no_cores_to_use)
  
}

CT.study.process.wrapper = function(st, unique_studies, study_IDs, meta_out, 
                                    files, filepaths, hdr_filepaths, labels, air_values, check_weird_airs,
                                    tmp_dcm_dir, out_dir_NII, out_dir_NII_mask, out_dir_NII_labels,
                                    only_labels = FALSE) {
  
  if (st %% 1000 == 0) { cat(st, ' ') }
  
  study_in = unique_studies[st]
  if (is.na(air_values[1])) {
    air_value = NA
  } else {
    air_value = air_values[st]
  }
  
  no_of_studies = length(unique_studies)
  matching_idxs = study_IDs %in% study_in
  matching_idxs_linear = which(matching_idxs)
  no_of_slices = sum(matching_idxs)
  
  # Match to thte desired study_ID
  meta_matched = match.all.the.fields(matching_idxs, meta_out, filepaths, files, hdr_filepaths, labels)

  # Sort the slices according to the z-position
  meta_sorted = sort.all.the.fields(meta_matched)
  z_resolution_vector = diff(meta_sorted[['z_pos']], lag = 1) # there are some rounding error
  z_resolution = median(abs(z_resolution_vector))
  
  if (!only_labels) {
    combine.dcms(meta_sorted, z_resolution, no_of_slices, study_in, no_of_studies, 
                 air_value, check_weird_airs, 
                 tmp_dcm_dir, out_dir_NII, out_dir_NII_mask)  
  }
  
  # write the labels also to disk for volume
  # if_intraparench = export.labels.for.volume(meta_sorted, study_in, out_dir_NII_labels)
  
  return(if_intraparench)
  
}

recheck.for.weird.air.baselines = function(unique_studies, LUT_file, correction_file,
                                           done_dir, mask_dir, 
                                           if_inspect = FALSE) {
  
  weirdos = read.csv(LUT_file)
  # TODO! Continue!
  corrs = read.csv(correction_file, stringsAsFactors = FALSE)
  corrs$fname = gsub(' ', '', corrs$fname)
  corrs$codes = gsub('.nii.gz', '', corrs$fname)
  
  codes_to_rerun = unique_studies %in% corrs$codes
  unique_studies = unique_studies[codes_to_rerun]
  
  idx_not_kept = which(!corrs$codes %in% unique_studies)
  
  if (length(idx_not_kept) > 0) {
    cat('WARNING! For the following files we had manually defined thresholds but these were not found from "unique_studies" list?\n')
    for (kpt in 1 : length(idx_not_kept)) {
      cat(' ... ', corrs$codes[idx_not_kept[kpt]], '\n')     
    }
    corrs = corrs[codes_to_rerun,]
  } else {
    unique_studies = corrs$codes
  }
  
  list_out = list()
  list_out$unique_studies = unique_studies
  list_out$corrs = corrs
  
  if (if_inspect) {
    for (w in 1 : length(weirdos[['code']])) {
      fname = paste0('ID_', weirdos[['code']][w], '.nii.gz')  
      nii = readNifti(file.path(done_dir, fname))
      air_value_guess = nii[1,1,1]
      cat(paste0('#', w, ': ', fname, ' | corner_value (1,1,1) = ', air_value_guess, '\n'))
      if (air_value_guess > -900) {
        plot(as.cimg(nii[,,1]))
        readline(prompt="Press [enter] to continue")
      }
    }
  }
  
  return(list_out)
  
}

export.labels.for.volume = function(meta_sorted, study_in, out_dir_NII_labels) {
  
  slice_by_slice_labels = meta_sorted$labels_matched
  names(slice_by_slice_labels) <- gsub("healthy", "normal", names(slice_by_slice_labels))
  df_labels_out = data.frame(slice_by_slice_labels)
  
  # WRITE
  fname_out = paste0(study_in, '_sliceLabels.csv')
  path_out = file.path(out_dir_NII_labels, fname_out)
  write.csv(df_labels_out, file = path_out, row.names = FALSE)
  
  # if contains intraparenchymal (as we are interested in these)
  if_intraparench = sum(df_labels_out$intraparenchymal) > 0

  return(if_intraparench)
  
}

combine.dcms = function(meta_sorted, z_resolution, no_of_slices, study_in, no_of_studies, 
                        air_values, check_weird_airs, 
                        tmp_dcm_dir, out_dir_NII, out_dir_NII_mask, 
                        fromDCM = TRUE,
                        convert_with_DCM2NIIX = FALSE, clip_on = TRUE,
                        clip_range = c(-1024, 2047)) {
  
  # check that all the slices have the same dims and spacing
  x_res = rep(NA, no_of_slices)
  y_res = rep(NA, no_of_slices)
  x_spacing = rep(NA, no_of_slices)
  y_spacing = rep(NA, no_of_slices)
  
  DCM_slices = list()
  DCM_hdrs = list()
  
  if (fromDCM) {
    cat(paste0('Importing ', no_of_slices, ' slices, for study "', study_in, '" (#', st, ' / ', no_of_studies, ' studies)\n'))  
  } else {
    cat(paste0('Importing ', no_of_slices, ' headers (RData), for study "', study_in, '" (#', st, ' / ', no_of_studies, ' studies)\n'))
  }
    
  for (s in 1 : no_of_slices) {
      
    if (fromDCM) {
      readError = FALSE
      out = try({
        DCM <- readDICOMFile(meta_sorted[['filepaths']][s])
        hdr = dcmImage$hdr
        DCM_slices[[s]] = DCM$img
        DCM_hdrs[[s]] = hdr
      }, silent = TRUE)
      
      if (class(out) == "try-error") {
        cat('  error with reading = ', files[f], ', \n') # "ID_6431af929.dcm"
        cat('  ', out) 
      } else {
        readError = TRUE
      }
    } else {
      
      load(meta_sorted[['hdr_filepath']][s])
      
    }
    
    # make sure that the dimensions match, or more precisely have the biggest
    # matching x,y dimensions for all the slices of the volume
    y_res[s] = as.integer(hdr$value[grep('Rows', hdr$name)])
    x_res[s] = as.integer(hdr$value[grep('Columns', hdr$name)])
    y_spacing[s] = as.numeric(strsplit(hdr$value[grep('PixelSpacing', hdr$name)], split = ' ')[[1]][1])
    x_spacing[s] = as.numeric(strsplit(hdr$value[grep('PixelSpacing', hdr$name)], split = ' ')[[1]][2])
  
  }
  
  unique_x = unique(x_res)
  unique_y = unique(y_res)
  unique_x_spac = unique(x_spacing)
  unique_y_spac = unique(y_spacing)
  
  funky_volume = FALSE
  more_than_one = length(unique_x) > 1 | length(unique_y) > 1 | length(unique_x_spac) > 1 | length(unique_y_spac) > 1
  
  if (more_than_one) {
    cat('! NEED TO PAD SOME SLICES for = ', study_in, '\n')
  } else {
    
    if (convert_with_DCM2NIIX) {
      
      # THE DESIRED OPTION, however there are errors from conversion
      # https://www.kaggle.com/anjum48/reconstructing-3d-volumes-from-metadata/data#714729
      # Move the desired .dcms to temp folder and convert to NIFTI
      create.NII.from.tmp.DCMs(meta_sorted, 
                               filepaths_tmp = meta_sorted$filepaths,
                               filenames_tmp = meta_sorted$files,
                               tmp_dcm_dir, out_dir_NII)
      
    } else {
      
      CT_array = array(NA, dim = c(unique_x, unique_y, no_of_slices))
      for (s in 1 : no_of_slices) {
        CT_array[,,s] = DCM_slices[[s]]
      }
      
      hdr = DCM_hdrs[[s]]
      
      # Check range here
      min_HU = min(CT_array)
      max_HU = max(CT_array)
      img_mask = CT_array != -2000
      if (clip_on) {
        
        # NOTE! some of the input .dcms seem weirdly scaled so that outside circular "image area",
        # the values are -2000 whereas air has value 0!
        
        # Quick'n'dirty re-scaling
        hh = hist(CT_array)
        # plot(as.cimg(CT_array[,,1]))
        mode_bin = which.max(hh$counts)

        # we assume to find the most values to be air (this might fail?)
        HU_mode = hh$breaks[mode_bin]
        if (!is.na(air_value)) {
          
          if (air_value == 9999) {
            
            # something funky about these 5 files, clipped values?
            funky_volume = TRUE
              
          } else {
          
            # weird values, manually checked (268 volumes in total)
            CT_array = CT_array - (1024+air_value)
            cat(paste0(' ... Manually checked AIR THRESHOLD of ', air_value, ' HU! -> fixed (', study_in, ')\n'))
            
            # now we can do the "standard winsorization"
            CT_array[CT_array < clip_range[1]] = clip_range[1]
            CT_array[CT_array > clip_range[2]] = clip_range[2]  
            
          }
          
          
          
        } else {
          
          if (HU_mode > -10 & HU_mode < 10) {
            
            # air is at 0
            # so we only rescale if the mode is found around 0, this will fail if the 
            # head is very tightly cropped and the air is at wrong place
            # TODO! come back to finding these samples later maybe with something more powerful finding
            # the out-of-distribution volumes
            CT_array = CT_array - 1024
            cat(paste0(' ... Air was at 0 HU! -> fixed (', study_in, ')\n'))
            
          }
          
          path_air = file.path(out_dir_NII, '..', 'Kaggle_air_histograms')
          dir.create(path_air, showWarnings = FALSE)
          save(hh, file = file.path(path_air, paste0(HU_mode, '_HUMode_', study_in, '_histogram.RData')))
          
          # now we can do the "standard winsorization"
          CT_array[CT_array < clip_range[1]] = clip_range[1]
          CT_array[CT_array > clip_range[2]] = clip_range[2]
        }
      }
     
      if (!funky_volume) {
      
        # header info mostly lost 
        nii_out = updateNifti(rotate(CT_array,90))
        pixdims = c(unique_x_spac, unique_y_spac, z_resolution)
        pixdim(nii_out) = pixdims
        pixunits(nii_out) = c('mm', 's')
        min_HU = min(nii_out)
        max_HU = max(nii_out)
        
        filename_out = paste0(study_in, '.nii.gz')
        filepath = file.path(out_dir_NII, filename_out)
        writeNifti(nii_out, file = filepath, datatype = 'int16')
        
        # write the mask as well
        filename_out = paste0(study_in, '_imageAreaMask.nii.gz')
        filepath = file.path(out_dir_NII_mask, filename_out)
        writeNifti(img_mask, file = filepath, datatype = 'int16')
        
      } else {
        cat(paste0('Do not write the "funky volume": "', study_in, '" to disk\n'))
      }
    
    }
  }
}

create.NII.from.tmp.DCMs = function(meta_sorted, filepaths_tmp, filenames_tmp,
                                    tmp_dcm_dir, out_dir_NII) {
  
  source('~/Dropbox/manuscriptDrafts/CThemorr/CT_hemorrhage/R_analysis/wrangle_Qure500.R')
  
  cat(paste0('   ... copying .dcsm to temp directory (', tmp_dcm_dir, ') \n'))
  temp_to = rep(NA, length(filepaths_tmp))
  for (f in 1 : length(filepaths_tmp)) {
    from = filepaths_tmp[f]
    temp_to[f] = file.path(tmp_dcm_dir, filenames_tmp[f])
    o = file.copy(from, temp_to[f])
  }
  
  cat(paste0('    ... creating NifTI\n'))
  convert.dir.into.NII(tmp_dcm_dir, dcm_path = '/home/petteri/mricron/dcm2niix_nonequi')
  # TODO! See:
  # https://www.kaggle.com/anjum48/reconstructing-3d-volumes-from-metadata/data#714729
  
  
  cat(paste0('    .... moving the created NifTI and the associated JSON to "', out_dir_NII, '")\n'))
  nii_file = list.files(path = tmp_dcm_dir, pattern = '*.nii')
  nii_path = list.files(path = tmp_dcm_dir, pattern = '*.nii', full.names = TRUE)
  o1 = file.rename(nii_path, file.path(out_dir_NII, nii_file))
  
  json_dir = file.path(out_dir_NII, 'JSON')
  dir.create(json_dir, showWarnings = FALSE)
  json_file = list.files(path = tmp_dcm_dir, pattern = '*.json')
  json_path = list.files(path = tmp_dcm_dir, pattern = '*.json', full.names = TRUE)
  o2 = file.rename(json_path, file.path(json_dir, json_file))
  
  cat(paste0('     ... removing temporary .dcms\n'))
  o3 = file.remove(temp_to)
  
}

match.all.the.fields = function(matching_idxs, meta_out, filepaths, files, hdr_filepaths, labels) {
  
  metadata_fields = names(meta_out$metadata)
  meta_matched = list()
  for (m in 1 : length(metadata_fields)) {
    field_name = metadata_fields[m]
    all_values = meta_out[['metadata']][[field_name]]
    meta_matched[['metadata']][[field_name]] = all_values[matching_idxs]
  }
  
  meta_matched[['imagePositions']] = meta_out[['imagePositions']][matching_idxs,]
  meta_matched[['imageOrientations']] = meta_out[['imageOrientations']][matching_idxs,]
  
  meta_matched[['filepaths']] = filepaths[matching_idxs]
  meta_matched[['files']] = files[matching_idxs]
  meta_matched[['hdr_filepath']] = hdr_filepaths[matching_idxs]
  
  # LABELS
  labels_train = labels$train$labels
  label_names = names(labels_train)
  labels_matched = list()
  
  for (l in 1 : length(label_names)) {
    nam = label_names[l]
    labels_matched[[nam]] = labels_train[[nam]][matching_idxs]
  }
  meta_matched[['labels_matched']] = labels_matched
  
  # ImagePositionPatient_2 -> as in the last coordinate (0,1,2) in 
  # Use easier variable names for readability
  ImagePositionPatient_2 = meta_out$imagePositions[,3]
  meta_matched[['z_pos']] = as.numeric(ImagePositionPatient_2[matching_idxs])
  
  return(meta_matched)
  
}

sort.all.the.fields = function(meta_matched) {
  
  # the largest value is at the bottom of the stack (i.e. the 1st slice)  
  sorted = sort(meta_matched[['z_pos']], decreasing = FALSE, index.return = TRUE)
  meta_sorted = list()
  
  names_in = names(meta_matched)
  for (n in 1 : length(names_in)) {
    
    name = names_in[n]
    
    if (name == 'metadata'  | name == 'labels_matched') {
      subnames = names(meta_matched[[name]])
      for (sb in 1 : length(subnames) ) {
        subname = subnames[sb]
        meta_sorted[[name]][[subname]] = meta_matched[[name]][[subname]][sorted[['ix']]]
      }
      
    } else if (name == 'imagePositions' | name == 'imageOrientations') {
      meta_sorted[[name]] = meta_matched[[name]][sorted[['ix']],]
      
    } else {
      meta_sorted[[name]] = meta_matched[[name]][sorted[['ix']]]
    }
      
  }
  
  return(meta_sorted)
  
}

get.metadata.from.listing = function(files, subdir, data_dir, IDs_from_files, hdr_path_out,
                                     write_hdr = FALSE, load_hdr = TRUE, 
                                     load_meta_from_RData = FALSE,
                                     skip_done_files = FALSE) {
  
  if (skip_done_files) {
    
      # for debugging as there might be this happening (especially if you have corrupted .dcms)
      # Error in parsePixelData(fraw[(bstart + dcm$data.seek):fsize], hdr, endian,  : 
      # Number of bytes in PixelData does not match dimensions of image. 
      # ONE .DCM corrupt: "ID_6431af929.dcm"
      cat('DEBUG MODE where we skip processing done files, but this will result the "RSNA_kaggle_metadata.Rdata" not containing all the entries\n')
      cat('   ... thus run this loop again with "load_hdr = TRUE" when you get this working\n')
      files = remove.done.files(files, IDs_from_files, done_dir = hdr_path_out)
  }
  
  if (load_meta_from_RData) {
    load(file.path(hdr_path_out, '..', 'RSNA_kaggle_metadata.Rdata'))
    
  } else {
  
    # TAKES FOREVER (disk I/O intensive, so no need to reprocess every time)
    # https://www.kaggle.com/anjum48/reconstructing-3d-volumes-from-metadata
    # + windowing https://www.kaggle.com/wfwiggins203/eda-dicom-tags-windowing-head-cts
    no_of_files = length(files)
    metadata = list()
    
    imagePositions = matrix(NA, nrow = no_of_files, ncol = 3)
    imageOrientations = matrix(NA, nrow = no_of_files, ncol = 6)
    
    # https://mattstats.wordpress.com/2013/05/22/dicom/
    # https://cran.r-project.org/web/packages/oro.dicom/vignettes/dicom.pdf
    fields_to_get = c('PatientID', 'StudyInstanceUID', 'SeriesInstanceUID', 'SOPInstanceUID')
    
    # f = seq(from = 1, to = no_of_files)
    # file.dcm.wrapper(f, files, no_of_files, hdr_path_out, data_dir, subdir,
    #                  load_hdr, write_hdr)
    filename_skeleton = file.path(hdr_path_out, '..', gsub('.dcm', '_skeleton.RData', files[1]))
    
    for (f in 1 : no_of_files) {
      
      if (f %% 1000 == 0 | f == 1 | f == no_of_files)
        
      cat('Reading in file #', f, '/', no_of_files, '(', files[f], ')\n')
      filename_hdr = file.path(hdr_path_out, gsub('.dcm', '.RData', files[f]))
      file_fullpath = file.path(data_dir, subdir, files[f])
      
      readError = FALSE
      if (!load_hdr) {
        out = try({
                dcmImage <- readDICOMFile(file_fullpath)
                hdr = dcmImage$hdr
              }, silent = TRUE)
    
        if (f == 1) {
          save(hdr, file = filename_skeleton)
        }
        
        if (class(out) == "try-error") {
          
          readError = TRUE
          cat('  error with reading = ', files[f], ', \n') # "ID_6431af929.dcm"
          cat('  ', out)
          
          # manual fix with just the header imported from ImageJ?
          # TudorDICOM_Plugin download where, so the header could be exported from ImageJ
          # hdr_file = gsub('.dcm', '.hdr', files2[f])
          # hdr_file = file.path(data_dir, hdr_file)
          # nifri(hdr_file)
          cat('MANUAL HDR EXPORT with relevant fields from ImageJ, NOTE! Not all metadata fields are correct now!\n')
          load(filename_skeleton)
          img = matrix(as.integer(-1024), nrow=512, ncol = 512) # empty
          hdr = construct.manual.header(hdr, fields_to_get)
          readError = FALSE
          
        }
        
      } else {
      
        load(file = filename_hdr)
          
      }
      
      if (readError) {
        
        # can't really process the metadata / header info as we could not read in the .dcm file
        
      } else {
      
        if (write_hdr) {
          save(hdr, file = filename_hdr)
        }
        
        # https://www.kaggle.com/c/rsna-intracranial-hemorrhage-detection/discussion/109552
        for (fld in 1 : length(fields_to_get)) {
          fieldname = fields_to_get[fld]
          idx = which(hdr$name %in% fieldname)
          
          value_in = hdr$value[idx]
          
          # TODO! Preallocate these, this is computationally bad!
          if (f == 1) {
            metadata[[fieldname]] = value_in
          } else {
            metadata[[fieldname]] = c(metadata[[fieldname]], value_in)
          }
          
        }
        
        idx = which(hdr$name %in% 'ImagePositionPatient') #14
        value_in = hdr$value[idx]
        value_in = strsplit(value_in, ' ')[[1]]
        imagePositions[f,] =  value_in
        idx = which(hdr$name %in% 'ImageOrientationPatient')
        value_in = hdr$value[idx]
        value_in = strsplit(value_in, ' ')[[1]]
        imageOrientations[f,] = value_in
      }
      
    }
    

    unique_values = list()
    for (fld in 1 : length(fields_to_get)) {
      fieldname = fields_to_get[fld]
      unique_IDs = unique(metadata[[fieldname]])
      unique_values[[fieldname]] = length(unique_IDs)
    }
  
    meta_out = list()
    meta_out[['metadata']] = metadata
    meta_out[['imagePositions']] = imagePositions
    meta_out[['imageOrientations']] = imageOrientations
    meta_out[['unique_values']] = unique_values
    
    save(meta_out, file = file.path(hdr_path_out, '..', 'RSNA_kaggle_metadata.Rdata'))
    
  }
  
  return(meta_out)
}

construct.manual.header = function(hdr, fields_to_get) {
  
  # MANUAL FIX for"ID_6431af929.dcm"
  
  idx = which(hdr$name %in% 'ImagePositionPatient') #14
  values_orig = hdr$value[idx] # replace the original ones with desired values
  hdr$value[idx] = '-125.000000 -127.746 174.901'
  
  idx = which(hdr$name %in% 'ImageOrientationPatient') #15
  values_orig = hdr$value[idx] # replace the original ones with desired values
  hdr$value[idx] = '1.000 0.000 0.000 0.000 0.972370 -0.233445'
  
  # 'PatientID', 'StudyInstanceUID', 'SeriesInstanceUID', 'SOPInstanceUID')
  new_values = c('ID_bba2045f', 'ID_9180c688de', 'ID_863be16ddb', 'ID_6431af929')
  
  for (fld in 1 : length(fields_to_get)) {
    fieldname = fields_to_get[fld]
    idx = which(hdr$name %in% fieldname)
    values_orig = hdr$value[idx] # replace the original ones with desired values
    hdr$value[idx] = new_values[fld]
  }
  
  return(hdr)
  
}

  

copy.wildcard.scans.as.NII.to.another.folder = function(raw_csv, wildcard, labels_to_use) {
  
  
  label_list = labels[[labels_to_use]]
  no_of_labels_from_csv = length(label_list$labels[[1]])
  cat(paste0(' .... we had labels for ', no_of_labels_from_csv, ' files from .CSV\n'))
  
  label_idxs = label_list[['labels']][[wildcard]] == 1
  label_count = label_list[['counts']][[wildcard]]
  
  IDs_to_pick = label_list[['ID']][label_idxs]
  percentage = round(100*length(IDs_to_pick) / no_of_labels_from_csv,2)
  cat(paste0(' ..... ', length(IDs_to_pick), ' (', percentage, '% of all labels) of which corresponded to wildcard = "', wildcard, '"\n'))
  
  unique_IDs = unique(IDs_to_pick)
  no_of_unique_IDS = length(unique_IDs)
  cat(paste0(' ...... and ', no_of_unique_IDS, ' of those are unique IDs\n'))
  
  matching_file_idx = files %in% unique_IDs
  no_of_files_to_process = sum(matching_file_idx)  
  
  
}


remove.done.files = function(files, IDs_from_files, done_dir) {
  
  files_done = list.files(done_dir)
  IDs_done = sapply(files_done, function(x) strsplit(gsub('.RData', '', x), "_")[[1]][2], USE.NAMES=FALSE)
  
  undone_idxs = !(IDs_from_files %in% IDs_done)
  files = files[undone_idxs]
  
  return(files)
  
}

check.for.done.studies = function(unique_studies, done_dir) {
  
  files_done = list.files(path = done_dir, pattern = '*.nii.gz')
  studies_done = gsub('.nii.gz', '', files_done)
  
  idx_done = unique_studies %in% studies_done
  unique_studies = unique_studies[!idx_done]
  
  return(unique_studies)
  
}

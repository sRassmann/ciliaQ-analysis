# Functions for joining and analyzing output from Cilia Q

# Version v0.0.1

require('dplyr')

get_required_packages <- function(){
  return(c('stringi','dplyr', 'rlist', 'readr'))
}

get_file_list <- function(root_dir, file_pattern = '.*_CQ_.*s\\.txt$', filter_pattern = NA){
  # match the defined pattern in all files within all dirs in the defined parent dir and return the list of results
  # if filter_pattern is defined, files will be selected for containg any of the pattern contained in the vector
  all = c()
  for(dir in list.dirs(root_dir)){
    l <- list.files(dir)
    l <- l[stringi::stri_detect_regex(l, file_pattern)]
    if(!is.na(filter_pattern)){
      l <- l[stringi::stri_detect_regex(l, filter_pattern)]
    }
    l <- file.path(dir,l)
    all <- c(all,l)
  }
  return(all)
}

read_data_from_ciliaQ <- function(file_list){
  # extracts the tabular ciliaQ output data in all given files and joins all data into one big data frame. Columns are renamed for a smooth usage in R, 
  # however, note, that this is oriented at v0.0.10 of CiliaQ, if new columns are added, this might mess up the naming!
  # removes suffixes of CiliaQ Premerger and Channel Splitter
  # checks for duplicated entries (same cilia ID in images with the same name)
  df <- data.frame()
  for(file in file_list){
    tryCatch( {
      df <- rbind(df, suppressWarnings( suppressMessages(readr::read_delim(file, "\t", col_types = readr::cols(), escape_double = FALSE, col_names = FALSE, 
                                                                           trim_ws = TRUE, progress = T)[,c(1:51)])))
    }, error = function(e) {
      print(paste("Error in opening file:", file," - This file was not included. Probably the file is empty, check and rerun function."))
    }, finally = function(){
    })
  }
  colnames(df) <- c('filename', 'cilia_ID_unique','cilia_ID_in_image','center_x', 'center_y', 'center_z', 'volume_vox', 'volume_um', 'n_surface_vox', 'surface_um',
                    'shape_complex', 'sphere_radius_um', 'max_span_um', 'coloc_vol_A_um', 'coloc_vol_A_perc', 'coloc_vol_B_um',
                    'coloc_vol_B_perc', 'coloc_vol_A_BG_um',  'coloc_vol_A_BG_perc', 'coloc_vol_B_BG_um',  'coloc_vol_B_BG_perc',
                    'min_int_reconst', 'max_int_reconst', 'av_int_reconst_upper_ten_perc', 'av_int_reconst', 'sd_int_reconst',    
                    'min_int_A', 'max_int_A', 'av_int_A_upper_ten_perc', 'av_int_A', 'sd_int_A',
                    'min_int_B', 'max_int_B', 'av_int_B_upper_ten_perc', 'av_int_B', 'sd_int_B',
                    'n_found_skeletons', 'n_branches', 'tree_lenght_um', 'cilia_length_um', 'orient_vect_x',
                    'orient_vect_y', 'orient_vect_z', 'cilia_bend', 'int_thres_A', 'int_thres_B','empty2',                   
                    'integrated_int_A', 'av_int_A_center_line', 'integrated_int_B', 'av_int_B_center_line')
  set     <- colnames(df)[!(colnames(df) %in% c('filename', 'cilia_ID_unique','cilia_ID_in_exp'))]
  df[set] <- apply(df[set],2, as.double)
  df$cilia_ID_in_image <- sprintf('c%04d', df$cilia_ID_in_image)
  df$filename <- gsub(df$filename, pattern = '_CS', replacement = '') # delete suffix of ChannelSplitter
  df$filename <- gsub(df$filename, pattern = '_M.tif', replacement = '') # delte suffix of CiliaQ Premerger
  df$cilia_ID_unique <- paste(df$filename, df$cilia_ID_in_image, sep = ' - ')
  if (any(duplicated(paste(df$filename, df$cilia_ID_in_image)))){ # if a combination of filename and cilia ID appears duplicated in the list, stop the script
    msg <- df %>% 
      filter(duplicated(cilia_name)) %>% 
      pull(filename) %>% 
      unique() # print the name of the image containing duplicated cilia
    stop(paste('There are duplicated files in the list, check the following files and remove duplicates: ', msg))
  } else {
    rownames(df) <- df$cilia_ID_unique
  }
  df$empty2 <- NULL
  return(df)
}

convert_vector_to_regex <- function(regex){
  # converts a list or vector of regex patterns into a conjunction
  return( paste0('(',paste(regex, collapse = ")|("), ')') )
}

plot_pca_eigenvalues <- function(prcomp_obj){
    plot <- data.frame(sd = prcomp_obj$sdev^2, PC = as.factor(c(1:length(prcomp_obj$sdev)))) %>% 
            mutate(sd_perc = sd/sum(sd) * 100) %>% 
            ggplot(aes(y = sd_perc, x = PC)) + 
            geom_histogram(stat = 'identity')  +
            ylab('% of overall variance') + 
            ggtitle('Variance per PC')
    return(plot)
}

correlation_plot <- function(df, x, y, color, title, alpha = 0.5) {
  # Takes and dataFrame and generates a scatter plot of the two passed arguments, colored by the color argument 
  plot <- ggplot(df, aes_string(x = x, y = y)) +      # must use aes_string as the variables are passed as strings and not unquoted!
    geom_point(aes_string(color = color), alpha = alpha)
  if (!missing(title)) {
    plot <- plot +  ggtitle(title)    
  }          
  # TODO add correlation coefficient per group!
  return(plot)
}

calculate_tSNE <- function(df, cols, perpexlity = 15, max_iter = 100, check_duplicates = F){
  # handles pre-filtering, calculation of the tsne object, and re-joining of the inserted dataFrame 
  tsne_df <- tidyr::drop_na(df)
  tsne <- as.matrix(tsne_df[,cols]) %>% 
    Rtsne::Rtsne(perplexity = perpexlity, max_iter = max_iter, check_duplicates = check_duplicates)
  tsne_anno <- cbind(tsne_df, tsne$Y) %>% 
    rename('tSNE_1' = '1', 'tSNE_2' = '2')
  return(tsne_anno)
}


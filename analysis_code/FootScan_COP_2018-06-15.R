# COP data processing functions

# =============================================================================

## Import libraries
library(Morpho)
library(plyr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RNiftyReg)


# =============================================================================

#' Import Footscan data
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param pressure_filepath String. Filepath pointing to pressure file generated
#'   by FootScan system
import_footscan <- function(pressure_filepath) {
  # check filepath
  if (is.character(pressure_filepath) == FALSE)
    stop("filepath needs to be a character string")
  
  # import data to list of matrices
  x <- readLines(pressure_filepath)
  right_foot <- grep("Right foot data", x)
  left_foot <- grep("Left foot data", x)
  right_foot_dim <- as.numeric(c(x[right_foot + 8], x[right_foot + 10]))
  left_foot_dim <- as.numeric(c(x[left_foot + 8], x[left_foot + 10]))
  
  if (right_foot_dim[1] < left_foot_dim[1]) {
    x <- x[left_foot:(right_foot - 2)]
    sensor_matrix_dim <- left_foot_dim
  }
  if (left_foot_dim[1] < right_foot_dim[1]) {
    x <- x[right_foot:length(x)]
    sensor_matrix_dim <- right_foot_dim
  }
  
  start_point <- grep("Frame 0", x)
  no_frames <- (length(x) - (start_point - 2)) / (sensor_matrix_dim[2] + 2)
  frames <- list()
  start <- seq(from = start_point + 1, by = (sensor_matrix_dim[2] + 2), 
               length.out = no_frames + 1)
  for (i in 1:no_frames) {
    y <- x[start[i]:(start[i + 1] - 3)]
    tc_y <- textConnection(y)
    y <- read.table(tc_y, sep = "\t")
    y <- y[2:(ncol(y) - 1)]
    frames[[i]] <- y
    close(tc_y)
  }
  return(frames)
}


# =============================================================================

#' Generate center of pressure coordinates
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param frames list of pressure values
#' @param sens_x_size Numeric. The size of each pressure sensor (x direction)
#' @param sens_y_size Numeric. The size of each pressure sensor (y direction)
#' @param interpol Numeric. Include if the COP is to be interpolated to a set
#'   number of timepoints
gen_cop <- function(frames, sens_x_size = 0.00508, sens_y_size = 0.00762, 
                    interpol) {
  # Loading totals by column
  col_total <- list()
  for (i in seq_along(frames)) {col_total[[i]] <- colSums(frames[[i]])}
  
  # Loading totals by row
  row_total <- list()
  for (i in seq_along(frames)) {row_total[[i]] <- rowSums(frames[[i]])}
  
  # Sensor spacing in x direction
  sensor_spacing_x <- seq(from = sens_x_size / 2, by = sens_x_size,
                          length.out = length(col_total[[1]]))
  
  # Sensor spacing in y direction
  sensor_spacing_y <- seq(from = sens_y_size / 2, by = sens_y_size,
                          length.out = length(row_total[[1]]))
  
  # COP coordinates in x direction
  x_coord <- c()
  for (i in seq_along(frames)) {
    p_total <- sum(col_total[[i]])
    x_coord[i] <- (sum(sensor_spacing_x * col_total[[i]])) / p_total 
  }
  
  # COP coordinates in y direction
  y_coord <- c()
  for (i in seq_along(frames)) {
    p_total <- sum(row_total[[i]])
    y_coord[i] <- (sum(sensor_spacing_y * row_total[[i]])) / p_total 
  }
  
  # interpolate if required
  if (hasArg(interpol) == TRUE) {
    x_coord <- approx(x_coord, n = interpol)
    x_coord <- x_coord$y
    y_coord <- approx(y_coord, n = interpol)
    y_coord <- y_coord$y
  }
  
  # combine coordinates into dataframe
  COP_df <- data.frame(x_coord, y_coord)
  
  # return COP coordinates
  return(COP_df)
}


# =============================================================================

# Caluculate mean center of pressure line using NiftyReg algorithm
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param template_trial String. Filepath to trial that is to be used as the template to be matched to
#' @param list_of_trials List of strings. List of filepaths to trials that to be used to find the mean COP
#' @param plot_cops Logical. If TRUE, function will plot footprint outlines and center of pressure lines
mean_cop2 <- function(template_trial, list_of_trials, plot_cops = FALSE) {
  # import template file
  template_frames <- import_footscan(template_trial)
  
  # template file COP
  template_cop <- gen_cop(template_frames, interpol = 101)
  
  # template file footprint
  frames_mat_list <- list()
  rowcol_n <- dim(template_frames[[1]])
  for(i in seq_along(template_frames)) {
    frames_mat_list[[i]] <- matrix(as.numeric(unlist(template_frames[[i]])), 
                                   nrow = rowcol_n[1], ncol = rowcol_n[2])
  }
  temp_fprint <- apply(simplify2array(frames_mat_list), 1:2, max)
  temp_fprint[temp_fprint > 0.1] <- 1
  
  # generate aligned footprint COPs and chull data frame
  cop_list_og <- list()
  cop_list_trans <- list()
  #outline_df <- data.frame(outline_name = as.factor(NA), x = as.numeric(NA),
  #                         y = as.numeric(NA))
  for (i in seq_along(list_of_trials)) {
    # import trial
    trial_frame <- import_footscan(list_of_trials[i])
    
    # cop
    cop_list_og[[i]] <- gen_cop(trial_frame, interpol = 101)
    
    # footprint
    frames_mat_list <- list()
    rowcol_n <- dim(trial_frame[[1]])
    for(j in seq_along(trial_frame)) {
      frames_mat_list[[j]] <- matrix(as.numeric(unlist(trial_frame[[j]])), 
                                     nrow = rowcol_n[1], ncol = rowcol_n[2])
    }
    trial_fprint <- apply(simplify2array(frames_mat_list), 1:2, max)
    trial_fprint[trial_fprint > 0.1] <- 1
    
    # Align footprints
    reg <- niftyreg.linear(trial_fprint, temp_fprint, scope = "rigid")
    nifty_transform <- forward(reg)
    
    # generate transformed cop
    cop_df <- data.frame(x_coord = as.numeric(), y_coord = as.numeric())
    for (j in 1:length(trial_frame)) {
      fd <- RNiftyReg::applyTransform(nifty_transform, as.matrix(trial_frame[[j]]))
      fd <- list(fd[1:nrow(fd), 1:ncol(fd)])
      cop_df[j, ] <- gen_cop(fd)
    }
    
    x_coord <- approx(cop_df$x_coord, n = 101)
    x_coord <- x_coord$y
    y_coord <- approx(cop_df$y_coord, n = 101)
    y_coord <- y_coord$y
    
    cop_df <- data.frame(x_coord = x_coord, y_coord = y_coord)
    
    
    #alignment <- align_footprint(template_frames, trial_frame)
    #trans_matrix <- alignment[[1]]
    #outlines <- alignment[[2]]
    #outlines$outline_name <- revalue(outlines$outline_name, c("matched" = 
    #                                                            paste0("matched_", i)))
    #outlines$outline_name <- revalue(outlines$outline_name, c("match" = 
    #                                                            paste0("match_", i)))
    #outline_df <- bind_rows(outline_df, outlines)
    
    # transform COP
    #trial_cop_mat <- as.matrix(cop_list_og[[i]])
    #cop_trans <- applyTransform(trial_cop_mat, trans_matrix, inverse = TRUE)
    cop_list_trans[[i]] <- cop_df
      #data.frame(x_coord = cop_trans[, 1], 
      #                                y_coord = cop_trans[, 2])
  }
  
  # fix outlines
  #outline_df <- outline_df %>% filter(complete.cases(.))
  #outline_df$outline_name <- as.factor(outline_df$outline_name)
  
  # calculate mean center of pressure line
  mean_cop <- aaply(laply(cop_list_trans, as.matrix), c(2, 3), mean)
  
  ## plot center of pressure lines
  if (plot_cops == TRUE) {
    # make dataframes of all center of pressure lines
    mean_cop_df <- data.frame(x = mean_cop[, 1], y = mean_cop[, 2])
    cop_df <- bind_rows(cop_list_trans) 
    cop_og_df <- bind_rows(cop_list_og)
    Trial <- rep(letters[1: length(cop_list_trans)], each = 101) 
    cop_df <- cbind(Trial, cop_df)
    cop_og_df <- cbind(Trial, cop_og_df)
    
    # plot center of pressure lines
    g <- ggplot()
    g <- g + geom_point(data = cop_df, aes(x = x_coord, y = y_coord, colour = Trial), size = 1)
    g <- g + geom_point(data = mean_cop_df, x = mean_cop_df$x, y = mean_cop_df$y, colour = "black")
    #g <- g + geom_path(data = outline_df, aes(x = outline_df$x, y = outline_df$y, colour = outline_df$outline_name))
    g <- g + geom_path(data = cop_og_df, aes(x = cop_og_df$x, y = cop_og_df$y, colour = Trial), size = 0.5)
    g <- g + coord_fixed()
    print(g)
  }
  
  # return mean center of pressure coordinates
  return(as.data.frame(mean_cop))
}


# =============================================================================

# Plot footprint with COP
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param x formated pressure file
#' @param show_cop Logical. If TRUE, center of pressure line will also be 
#'   plotted
rsscan_plot <- function(x, show_cop = TRUE, show_outline = TRUE) {
  # find max footprint of template
  frames_mat_list <- list()
  rowcol_n <- dim(x[[1]])
  for(i in seq_along(x)) {
    frames_mat_list[[i]] <- matrix(as.numeric(unlist(x[[i]])), 
                                   nrow = rowcol_n[1], ncol = rowcol_n[2])
  }
  max_mat <- apply(simplify2array(frames_mat_list), 1:2, max)
  
  # sensor coordinates
  x_cor <- seq(from = 0.00508/2, by = 0.00508, length.out = ncol(max_mat))
  x_cor <- rep(x_cor, each = nrow(max_mat))
  y_cor <- seq(from = 0.00762/2, by = 0.00762, length.out = nrow(max_mat))
  y_cor <- rep(y_cor, times = ncol(max_mat))
  cor <- cbind(x_cor, y_cor)
  cor <- cbind(x_cor, y_cor, as.vector(max_mat))
  cor <- as.data.frame(cor)
  colnames(cor) <- c("x", "y", "z")
  cor$z <- (cor$z / (0.005 * 0.007)) / 1000
  colour <- c()
  
  # assign colours
  for (i in 1:(ncol(max_mat) * nrow(max_mat))) {
    if (cor$z[i] < 15) {colour = append(colour, 8)}
    if (cor$z[i] >= 15 & cor$z[i] < 40) {colour = append(colour, 1)}
    if (cor$z[i] >= 40 & cor$z[i] < 60) {colour = append(colour, 4)}
    if (cor$z[i] >= 60 & cor$z[i] < 100) {colour = append(colour, 5)}
    if (cor$z[i] >= 100 & cor$z[i] < 150) {colour = append(colour, 3)}
    if (cor$z[i] >= 150 & cor$z[i] < 220) {colour = append(colour, 7)}
    if (cor$z[i] >= 220 & cor$z[i] < 300) {colour = append(colour, 2)}
    if (cor$z[i] >= 300) {colour = append(colour, 6)}
  }
  cor <- cbind(cor, colour)
  # colours
  cols <- c("1" = "blue","2" = "orange", "3" = "green", "4" = "light blue",
            "5" = "light green", "6" = "red", "7" = "Yellow", "8" = "grey")
  
  # cop data
  cop_df <- gen_cop(x)
  
  # determine outline of max footprint of template
  outline <- footprint_outline(x)
  #P <- c(max_mat)
  #x_cor <- seq(from = 0.00508/2, by = 0.00508, length.out = ncol(max_mat))
  #x_cor <- rep(x_cor, each = nrow(max_mat))
  #y_cor <- seq(from = 0.00762/2, by = 0.00762, length.out = nrow(max_mat))
  #y_cor <- rep(y_cor, times = ncol(max_mat))
  #outline_df <- data.frame(x = x_cor, y = y_cor, P = P)
  #outline_df <- outline_df[which(P >= 0.1), ]
  #outline_df$P <- NULL
  
  # determine convex hull of template
  #chull_elements <- chull(x = outline_df$x, y = outline_df$y)
  #chull_temp_df <- outline_df[chull_elements, ]
  
  #plot data
  g <- ggplot()
  g <- g + geom_tile(data = cor, aes(x = x, y = y, fill = as.factor(colour)))
  g <- g + scale_fill_manual(values = cols)
  g <- g + geom_point(data = cop_df, aes(x = x_coord, y = y_coord), colour = "black")
  g <- g + geom_path(data = outline, aes(x = x, y = y), colour = "black")
  g <- g + coord_fixed()
  g <- g + theme_bw() + theme(legend.position = "none")
  print(g)
}


# =============================================================================

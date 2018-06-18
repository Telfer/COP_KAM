## Analysis of below shoe data from Dose Response study 
# Aims to determine if changes in the center of pressure measurement during the
# first half of stance are correlated with changes in the external knee
# adduction moment

# Notes 
# Primarily uses functions from FootScan_COP_2018-06-15.R file
# To use, rename "all_data_folder" and "RFP_KAM_Results3.csv" path to your 
# local directory

# =============================================================================


# Identify data files
subjects <- c("RFPC01", "RFPC02", "RFPC04", "RFPC05", "RFPC06", "RFPC07", 
              "RFPC08", "RFPC09", "RFPC10", "RFPC11", "RFPC13", "RFPP01",
              "RFPP02", "RFPP03", "RFPP04", "RFPP05", "RFPP08",
              "RFPP10", "RFPP11", "RFPP13")
conditions <- c("L06", "L04", "L02", "N00", "M02", "M04", "M06", "M08",
                "M10", "Shod")
sides <- c("left", "right", "right", "left", "left", "right", 
           "right", "right", "right", "right", "left", "right", 
           "right", "right", "right", "left", "left", 
           "right", "left", "right")
all_data_folder <- paste0("C:/Users/telfe/Dropbox/My_Grant_", 
                          "Applications/LWI_Fit/Preliminary data")


# For each subject, generate dataframe that is each trial aligned to baseline 
# footprint
list_of_cops <- list()
for (i in seq_along(subjects)) {
  subject_folder <- file.path(all_data_folder, subjects[i])
  colu <- rep(NA, times = 101)
  cond_df <- data.frame("L6_x" = colu, 
                        "L6_y" = colu, "L4_x" = colu, "L4_y" = colu, 
                        "L2_x" = colu, "L2_y" = colu, "N0_x" = colu, 
                        "N0_y" = colu, "M2_x" = colu, "M2_y" = colu, 
                        "M4_x" = colu, "M4_y" = colu, "M6_x" = colu, 
                        "M6_y" = colu, "M8_x" = colu, "M8_y" = colu, 
                        "M10_x" = colu, "M10_y" = colu,
                        "Shod_x" = colu, "Shod_y" = colu)
  
  # Get overall template
  folder <- paste0(subject_folder, "/N00")
  trials <- list.files(folder, full.names = TRUE)
  template <- trials[1]
  
  # trials by condition
  for (j in seq_along(conditions)) {
    folder <- paste0(subject_folder, "/", conditions[j])
    trials <- list.files(folder, full.names = TRUE)
    
    # generate mean center of pressure line
    cop_df <- mean_cop2(template, trials)
    cond_df[, ((j * 2) - 1)] <- cop_df$x_coord 
    cond_df[, (j * 2)] <- cop_df$y_coord
  }
    
  # save to list of all center of pressure lines
  list_of_cops[[i]] <- cond_df
}

## correlate with external knee adduction moment measurements
# import knee adduction moment data
KAM_df <- read.csv(paste0("C:/Users/telfe/Dropbox/My_Grant_Applications/", 
                          "LWI_Fit/Preliminary data/RFP_KAM_Results3.csv"))
KAM1_df <- KAM_df %>% filter(Variable == "KAM1")
KAM1_df <- KAM1_df %>% filter(Subject %in% subjects)

# cop position at time of max KAM1
maxKAM1 <- list()
for (i in seq_along(list_of_cops)) {
  KAM1_Pos_df <- KAM1_df %>% filter(Subject == subjects[i]) %>% select(Position)
  cop_df <- list_of_cops[[i]] %>%
    select(L6_x, L4_x, L2_x, N0_x, M2_x, M4_x, M6_x, M8_x, M10_x, Shod_x)
  pos_x <- c()
  for (j in seq_along(conditions)) {
    pos <- KAM1_Pos_df[j, ]
    pos_x <- c(pos_x, cop_df[pos, j])
  }
  
  pos_x <- pos_x - pos_x[4]
  if (sides[i] == "left") {maxKAM1[[i]] = pos_x}
  if (sides[i] == "right") {maxKAM1[[i]] = pos_x * -1}
}

# combine KAM data with COP data
combined_data <- list()
for (i in seq_along(subjects)) {
  KAM_vec <- KAM1_df %>% filter(Subject == subjects[i]) %>%
    select(Value)
  KAM_vec <- c(KAM_vec[1:10, 1])
  COP_vec <- as.vector(unlist(maxKAM1[[i]]))
  
  combined_data[[i]] <- data.frame(KAM = KAM_vec, COP = COP_vec)
}

for (i in seq_along(combined_data)) {
  # COP vs KAM scatterplot
  g <- ggplot(combined_data[[i]], aes(x = KAM, y = COP))
  g <- g + geom_point()
  g <- g + theme_bw()
  g <- g + labs(title = subjects[i])
  print(g)
  
  # stacked plot showing both
  df <- combined_data[[i]]
  df <- gather(df, "Measurement")
  df <- cbind(conditions[1:10], df)
  df$conditions <- factor(df$conditions, 
                          levels(df$conditions)[c(3, 2, 1, 9, 4, 5, 6, 7, 8, 10)])
  gg <- ggplot(df, aes(x = conditions, y = value))
  gg <- gg + geom_point()
  gg <- gg + facet_grid(Measurement ~ ., scales = "free")
  gg <- gg + ggtitle(subjects[i])
  print(gg)
}


## ============================================================================

## statistal modeling 
library(nlme)

## Genrate data frames with and without shod
for (i in seq_along(subjects)) {
  combined_data[[i]]$KAM <- combined_data[[i]]$KAM - combined_data[[i]]$KAM[4]
}
df <- bind_rows(combined_data)
group <- c(rep("control", times = 10 * 11), rep("patient", times = 10 * 9))
df <- cbind(df, group)
Subjects <- rep(subjects, each = 10)
Conditions <- rep(conditions[1:10], times = length(subjects))
df <- cbind(Subjects, Conditions, df)
df_shod <- df 
df_orthotic <- df %>% filter(Conditions != "Shod")

## Determine correlations
# individual
for (i in seq_along(subjects)) {
  sub <- subjects[i]
  df_ind <- df_orthotic %>% filter(Subjects == sub) 
  print(cor(df_ind$KAM, df_ind$COP))
}

# group
cor(df_orthotic$KAM, df_orthotic$COP)

# control group
df_con <- df_orthotic %>% filter(group == "control")
cor(df_con$KAM, df_con$COP)

# patient group
df_pat <- df_orthotic %>% filter(group == "patient")
cor(df_pat$KAM, df_pat$COP)

# check distribution of KAM1
hist(df_orthotic$KAM)

# run basic linear model
basic.lm <- lm(KAM ~ COP, data = df_orthotic)
summary(basic.lm)

# check residuals
plot(basic.lm, which = 1)

# check residuals
plot(basic.lm, which = 2)

# using nlme
mixed.lm_cop_null <- lme(KAM ~ 1 + Conditions + group, random = ~1 | Subjects/Conditions, 
                         data = df_orthotic, method = "ML")
mixed.lm_cop <- lme(KAM ~ COP + Conditions + group, random = ~1 | Subjects/Conditions, 
                    data = df_orthotic, method = "ML")
anova(mixed.lm_cop_null, mixed.lm_cop)

# For groups
group.lm <- lm(KAM ~ COP*group, data = df_orthotic)
anova(group.lm)

# shod vs L4
df_shod_s <- df_shod %>% filter(Conditions %in% c("Shod", "L04", "L06")) 


# =============================================================================

## Figures
# Figure 1: Footprint with COP overlay
example_directory <- paste0(all_data_folder, "/", subjects[11])
template <- list.files(paste0(example_directory, "/N00"), full.names = TRUE)[1]
x <- import_footscan(template)
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
cols <- c("1" = "grey90","2" = "grey40", "3" = "grey60", "4" = "grey80",
          "5" = "grey70", "6" = "grey30", "7" = "grey50", "8" = "white")

# cop data
x <- list_of_cops[[11]]

#plot data
g <- ggplot()
g <- g + geom_tile(data = cor, aes(x = x, y = y, fill = as.factor(colour)))
g <- g + scale_fill_manual(values = cols)
g <- g + geom_point(aes(x = x$Shod_x, y = x$Shod_y), colour = "black", size = 1)
g <- g + geom_point(aes(x = x$L6_x, y = x$L6_y), colour = "red", size = 1)
g <- g + geom_point(aes(x = x$L4_x, y = x$L4_y), colour = "pink", size = 1)
g <- g + geom_point(aes(x = x$L2_x, y = x$L2_y), colour = "orange", size = 1)
g <- g + geom_point(aes(x = x$N0_x, y = x$N0_y), colour = "yellow", size = 1)
g <- g + geom_point(aes(x = x$M2_x, y = x$M2_y), colour = "green", size = 1)
g <- g + geom_point(aes(x = x$M4_x, y = x$M4_y), colour = "cyan", size = 1)
g <- g + geom_point(aes(x = x$M6_x, y = x$M6_y), colour = "light blue", size = 1)
g <- g + geom_point(aes(x = x$M8_x, y = x$M8_y), colour = "blue", size = 1)
g <- g + geom_point(aes(x = x$M10_x, y = x$M10_y), colour = "purple", size = 1)
g <- g + theme_bw() + theme(axis.line=element_blank(),axis.text.x=element_blank(),
                            axis.text.y=element_blank(),axis.ticks=element_blank(),
                            axis.title.x=element_blank(),
                            axis.title.y=element_blank(),legend.position="none",
                            panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                            panel.grid.minor=element_blank(),plot.background=element_blank())
g <- g + coord_fixed()
g
ggsave("C:/Users/telfe/Dropbox/My_Projects/LWI_Fit/papers/COP_KAM/Figure 1.png", 
       g, width = 8, height = 8, dpi = 600)

# Figure 2: scatterplot KAM vs COP
g <- ggplot(df, aes(KAM, COP, colour = group))
g <- g + geom_point()
g <- g + geom_smooth(method = "lm")
g <- g + theme_bw()
g <- g + scale_colour_discrete(name = "Group",
                               breaks = c("control", "patient"),
                               labels = c("Control", "Patient"))
g <- g + labs(x = "KAM (%BW*H)", y = "COP (m)", legend.title = "Group")
g
ggsave("C:/Users/telfe/Dropbox/My_Projects/LWI_Fit/papers/COP_KAM/Figure 2.png", 
       g, width = 8, height = 8, dpi = 600)


###############################################################################
#------------------------------------------------------------------------------
##################################### END #####################################
#------------------------------------------------------------------------------
###############################################################################
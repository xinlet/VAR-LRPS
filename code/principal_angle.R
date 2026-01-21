library(ggplot2)

#lrps====
estimates_dir <- "~/Downloads/task-fmri/estimates"
estimates_dir <- "~/Downloads/task-fmri/estimates_cuedts"
estimates_dir <- "~/Downloads/movie_fmri/estimates"
estimates_dir <- "~/Downloads/resting_fmri/estimates"

# List all LRPS-related .csv files
lrps_files <- list.files(estimates_dir, pattern = "LRPS.*\\.csv$", full.names = TRUE)

# Step 1: Load reference matrix
reference_matrix <- as.matrix(read.csv(lrps_files[1], header = TRUE))
svd_ref <- svd(reference_matrix)
d <- 1
V_ref <- svd_ref$v[, 1:d]

# Step 2: Loop through remaining matrices and compute principal angles
principal_angles_deg <- numeric(length(lrps_files) - 1)
file_names <- basename(lrps_files)[-1]

for (i in 2:length(lrps_files)) {
  mat <- as.matrix(read.csv(lrps_files[i], header = TRUE))
  svd_mat <- svd(mat)
  V_mat <- svd_mat$v[, 1:d]
  
  # Compute inner product and principal angle
  inner_prod <- t(V_mat) %*% V_ref
  sigma <- svd(inner_prod)$d
  angle_rad <- acos(sigma)  
  principal_angles_deg[i - 1] <- angle_rad
}

# Step 3: Plot results
df <- data.frame(
  Subject = paste0("sub ", seq_along(principal_angles_deg)),
  PrincipalAngle = principal_angles_deg
)

ggplot(df, aes(x = Subject, y = PrincipalAngle)) +
  geom_segment(aes(x = Subject, xend = Subject, y = 0, yend = PrincipalAngle),
               color = "gray", size = 1) +
  geom_point(shape = 21, fill = "black") +
  xlab("Subject") +
  ylab("Principal Angle") +
  theme_minimal() +
  ylim(0,  1.7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.text.x = element_blank(),
      axis.ticks.x = element_blank())


mean(principal_angles_deg)


pdf("lrps_angle.pdf", width = 3.5, height = 1.5)
ggplot(df, aes(x = Subject, y = PrincipalAngle)) +
  geom_segment(aes(x = Subject, xend = Subject, y = 0, yend = PrincipalAngle),
               color = "gray", size = 1) +
  geom_point(shape = 21, fill = "black") +
  xlab("Subject") +
  ylab("Principal Angle") +
  theme_minimal() +
  ylim(0,  1.7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
dev.off()

#svar====

# List all related .csv files
lrps_files <- list.files(estimates_dir, pattern = "SVAR.*\\.csv$", full.names = TRUE)

# Step 1: Load reference matrix
reference_matrix <- as.matrix(read.csv(lrps_files[1], header = TRUE))
svd_ref <- svd(reference_matrix)
d <- 1
V_ref <- svd_ref$v[, 1:d]

# Step 2: Loop through remaining matrices and compute principal angles
principal_angles_deg <- numeric(length(lrps_files) - 1)
file_names <- basename(lrps_files)[-1]

for (i in 2:length(lrps_files)) {
  mat <- as.matrix(read.csv(lrps_files[i], header = TRUE))
  svd_mat <- svd(mat)
  V_mat <- svd_mat$v[, 1:d]
  
  # Compute inner product and principal angle
  inner_prod <- t(V_mat) %*% V_ref
  sigma <- svd(inner_prod)$d
  angle_rad <- acos(sigma)  
  principal_angles_deg[i - 1] <- angle_rad
}

# Step 3: Plot results
df <- data.frame(
  Subject = paste0("sub ", seq_along(principal_angles_deg)),
  PrincipalAngle = principal_angles_deg
)

ggplot(df, aes(x = Subject, y = PrincipalAngle)) +
  geom_segment(aes(x = Subject, xend = Subject, y = 0, yend = PrincipalAngle),
               color = "gray", size = 1) +
  geom_point(shape = 21, fill = "black") +
  xlab("Subject") +
  ylab("Principal Angle") +
  theme_minimal() +
  ylim(0,  1.7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

mean(principal_angles_deg)

pdf("svar_angle.pdf", width = 3.5, height = 1.5)
ggplot(df, aes(x = Subject, y = PrincipalAngle)) +
  geom_segment(aes(x = Subject, xend = Subject, y = 0, yend = PrincipalAngle),
               color = "gray", size = 1) +
  geom_point(shape = 21, fill = "black") +
  xlab("Subject") +
  ylab("Principal Angle") +
  theme_minimal() +
  ylim(0,  1.7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
dev.off()


#var====

# List all related .csv files
lrps_files <- list.files(estimates_dir, pattern = "VAR1.*\\.csv$", full.names = TRUE)

# Step 1: Load reference matrix
reference_matrix <- as.matrix(read.csv(lrps_files[1], header = TRUE))
svd_ref <- svd(reference_matrix)
d <- 1
V_ref <- svd_ref$v[, 1:d]

# Step 2: Loop through remaining matrices and compute principal angles
principal_angles_deg <- numeric(length(lrps_files) - 1)
file_names <- basename(lrps_files)[-1]

for (i in 2:length(lrps_files)) {
  mat <- as.matrix(read.csv(lrps_files[i], header = TRUE))
  svd_mat <- svd(mat)
  V_mat <- svd_mat$v[, 1:d]
  
  # Compute inner product and principal angle
  inner_prod <- t(V_mat) %*% V_ref
  sigma <- svd(inner_prod)$d
  angle_rad <- acos(sigma)  
  principal_angles_deg[i - 1] <- angle_rad
}

# Step 3: Plot results
df <- data.frame(
  Subject = paste0("sub ", seq_along(principal_angles_deg)),
  PrincipalAngle = principal_angles_deg
)

ggplot(df, aes(x = Subject, y = PrincipalAngle)) +
  geom_segment(aes(x = Subject, xend = Subject, y = 0, yend = PrincipalAngle),
               color = "gray", size = 1) +
  geom_point(shape = 21, fill = "black") +
  xlab("Subject") +
  ylab("Principal Angle") +
  theme_minimal() +
  ylim(0,  1.7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

mean(principal_angles_deg)

pdf("var_angle.pdf", width = 3.5, height = 1.5)
ggplot(df, aes(x = Subject, y = PrincipalAngle)) +
  geom_segment(aes(x = Subject, xend = Subject, y = 0, yend = PrincipalAngle),
               color = "gray", size = 1) +
  geom_point(shape = 21, fill = "black") +
  xlab("Subject") +
  ylab("Principal Angle") +
  theme_minimal() +
  ylim(0,  1.7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
dev.off()


#Hybrid====

# List all related .csv files
lrps_files <- list.files(estimates_dir, pattern = "Hybrid.*\\.csv$", full.names = TRUE)

# Step 1: Load reference matrix
reference_matrix <- as.matrix(read.csv(lrps_files[1], header = TRUE))
svd_ref <- svd(reference_matrix)
d <- 1
V_ref <- svd_ref$v[, 1:d]

# Step 2: Loop through remaining matrices and compute principal angles
principal_angles_deg <- numeric(length(lrps_files) - 1)
file_names <- basename(lrps_files)[-1]

for (i in 2:length(lrps_files)) {
  mat <- as.matrix(read.csv(lrps_files[i], header = TRUE))
  svd_mat <- svd(mat)
  V_mat <- svd_mat$v[, 1:d]
  
  # Compute inner product and principal angle
  inner_prod <- t(V_mat) %*% V_ref
  sigma <- svd(inner_prod)$d
  angle_rad <- acos(sigma)  
  principal_angles_deg[i - 1] <- angle_rad
}

# Step 3: Plot results
df <- data.frame(
  Subject = paste0("sub ", seq_along(principal_angles_deg)),
  PrincipalAngle = principal_angles_deg
)

ggplot(df, aes(x = Subject, y = PrincipalAngle)) +
  geom_segment(aes(x = Subject, xend = Subject, y = 0, yend = PrincipalAngle),
               color = "gray", size = 1) +
  geom_point(shape = 21, fill = "black") +
  xlab("Subject") +
  ylab("Principal Angle") +
  theme_minimal() +
  ylim(0,  1.7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

mean(principal_angles_deg)





pdf("hybrid_angle.pdf", width = 3.5, height = 1.5)
ggplot(df, aes(x = Subject, y = PrincipalAngle)) +
  geom_segment(aes(x = Subject, xend = Subject, y = 0, yend = PrincipalAngle),
               color = "gray", size = 1) +
  geom_point(shape = 21, fill = "black") +
  xlab("Subject") +
  ylab("Principal Angle") +
  theme_minimal() +
  ylim(0,  1.7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
dev.off()

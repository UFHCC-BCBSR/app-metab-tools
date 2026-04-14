set.seed(42)

n_samples <- 60
n_features <- 1000

batch <- rep(c("Batch_A", "Batch_B", "Batch_C"), each = 20)
bio   <- rep(rep(c("Control", "Treatment"), each = 10), times = 3)

# Tight base matrix on log scale
log_mat <- matrix(
  rnorm(n_features * n_samples, mean = 10, sd = 0.3),
  nrow = n_features,
  ncol = n_samples
)
rownames(log_mat) <- paste0("feature_", 1:n_features)
colnames(log_mat) <- paste0("sample_", 1:n_samples)

# Batch effect on log scale — large, affects all features
# Different offset per batch, plus per-feature random scaling
batch_offset <- ifelse(batch == "Batch_A", 0,
                       ifelse(batch == "Batch_B", 3.0, 6.0))

# Per-feature batch sensitivity — some features more affected than others
batch_sensitivity <- rnorm(n_features, mean = 1, sd = 0.4)
batch_effect <- outer(batch_sensitivity, batch_offset)

log_mat <- log_mat + batch_effect

# Biological effect on log scale — real but smaller than batch
# Affects all features but with smaller offset
bio_offset <- ifelse(bio == "Treatment", 4.0, 0)
bio_responders <- sample(1:n_features, size = 300)
bio_sensitivity <- rep(0, n_features)
bio_sensitivity[bio_responders] <- abs(rnorm(300, mean = 2, sd = 1))
bio_effect <- outer(bio_sensitivity, bio_offset)
log_mat <- log_mat + bio_effect

# Back-transform to raw intensities for the app input
raw_mat <- 2^log_mat

# Add ~5% missing values
missing_idx <- sample(length(raw_mat), size = round(length(raw_mat) * 0.05))
raw_mat[missing_idx] <- NA

count_df <- tibble::rownames_to_column(as.data.frame(raw_mat), var = "feature_id")
metadata_df <- data.frame(
  sample_name    = paste0("sample_", 1:n_samples),
  batch          = batch,
  biological_var = bio,
  stringsAsFactors = FALSE
)

write.csv(count_df,    "test-data/counts-data-metabo.csv",  row.names = FALSE)
write.csv(metadata_df, "test-data/sample-data-metabo.csv",  row.names = FALSE)
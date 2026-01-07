# Target sum: n*(n+1)*(2n+1)/6
sum_sq_1n <- function(n) {
  n * (n + 1) * (2 * n + 1) / 6
}

# 1. The Colleague's Matrix Approach
# Uses O(M*N) matrix product and O(M*N) rbind copy
benchmark_matrix <- function(perm_list, n) {
  nc <- sum_sq_1n(n)

  # Initialize with first permutation
  closed <- matrix(perm_list[[1]], nrow = 1)
  unique_count <- 1

  for (i in 2:length(perm_list)) {
    pxc <- perm_list[[i]]

    # Dot product with every seen permutation
    # Complexity: O(unique_perms * n)
    dots <- closed %*% pxc

    if (max(dots) < nc) {
      # Append new permutation
      # Complexity: O(unique_perms * n) due to copy-on-write
      closed <- rbind(closed, pxc)
      unique_count <- unique_count + 1
    }
  }
  return(unique_count)
}

# 2. The Hashed Environment Approach
# Uses O(N) string key and O(1) hash lookup
benchmark_hash <- function(perm_list, n) {
  closed <- new.env(hash = TRUE, parent = emptyenv(), size = 2000L)
  unique_count <- 0

  for (i in seq_along(perm_list)) {
    pxc <- perm_list[[i]]

    # Create key: O(n)
    key <- paste(pxc, collapse = ",")

    # Hash lookup: O(1)
    if (!exists(key, envir = closed, inherits = FALSE)) {
      closed[[key]] <- TRUE
      unique_count <- unique_count + 1
    }
  }
  return(unique_count)
}

# --- Execution & Benchmarking ---
# Setup dummy data matching your description: n=105, 2000 unique, with repeats
n <- 105
unique_perms <- replicate(2000, sample(n), simplify = FALSE)
# Create a sequence with 50% repeats
test_data <- sample(unique_perms, 4000, replace = TRUE)

cat("Timing Colleague (Matrix):\n")
print(system.time(benchmark_matrix(test_data, n)))

cat("\nTiming Proposal (Hash):\n")
print(system.time(benchmark_hash(test_data, n)))
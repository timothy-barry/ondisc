########
# Test 1
########
test_that("extract single, contiguous submatrix", {
  for (i in 1:n_datasets) {
    Mat <- r_mats[[i]]
    on_disc_mat <- cov_odms[[i]]@ondisc_matrix
    for (j in 1:n_reps) {
      # Generate consecutive indices
      col_idxs_range <- sample(x = 1:ncol(Mat), size = 2, replace = FALSE) %>% sort()
      col_idxs <- col_idxs_range[1]:col_idxs_range[2]
      row_idxs_range <- sample(x = 1:nrow(Mat), size = 2, replace = FALSE) %>% sort()
      row_idxs <- row_idxs_range[1]:row_idxs_range[2]
      # Compare on ordered indexes
      compare_Mat_on_disc_extract(Mat = Mat, on_disc_mat = on_disc_mat, col_idxs = col_idxs, row_idxs = row_idxs)
      # Compare on randomly shuffled indexes
      col_idxs <- sample(col_idxs)
      row_idxs <- sample(row_idxs)
      compare_Mat_on_disc_extract(Mat = Mat, on_disc_mat = on_disc_mat, col_idxs = col_idxs, row_idxs = row_idxs)
    }
  }
})

########
# Test 2
########
test_that("extract multiple, contiguous submatrices", {
  for (i in 1:n_datasets) {
    Mat <- r_mats[[i]]
    on_disc_mat <- cov_odms[[i]]@ondisc_matrix
    if (nrow(Mat) >= 6 && ncol(Mat) >= 6) {
      for (j in 1:n_reps) {
        # Generate consecutive indices
        col_idxs_range <- sample(x = 1:ncol(Mat), size = 6, replace = FALSE) %>% sort()
        col_idxs <- c(col_idxs_range[1]:col_idxs_range[2], col_idxs_range[3]:col_idxs_range[4], col_idxs_range[5]:col_idxs_range[6])
        row_idxs_range <- sample(x = 1:nrow(Mat), size = 6, replace = FALSE) %>% sort()
        row_idxs <- c(row_idxs_range[1]:row_idxs_range[2], row_idxs_range[3]:row_idxs_range[4], row_idxs_range[5]:row_idxs_range[6])
        # extract sub-matrix by column
        compare_Mat_on_disc_extract(Mat = Mat, on_disc_mat = on_disc_mat, col_idxs = col_idxs, row_idxs = row_idxs)
        # Compare on randomly shuffled indexes
        col_idxs <- sample(col_idxs)
        row_idxs <- sample(row_idxs)
        compare_Mat_on_disc_extract(Mat = Mat, on_disc_mat = on_disc_mat, col_idxs = col_idxs, row_idxs = row_idxs)
      }
    }
  }
})

########
# Test 3
########
test_that("Extract arbitrary submatrices", {
  for (i in 1:n_datasets) {
    Mat <- r_mats[[i]]
    on_disc_mat <- cov_odms[[i]]@ondisc_matrix
    for (j in 1:n_reps) {
      subset_size_col <- sample(1:(ceiling(ncol(Mat)/30)), 1)
      subset_size_row <- sample(1:(ceiling(nrow(Mat)/30)), 1)
      col_idxs <- sample(x = 1:ncol(Mat), size = subset_size_col)
      row_idxs <- sample(x = 1:nrow(Mat), size = subset_size_row)
      compare_Mat_on_disc_extract(Mat = Mat, on_disc_mat = on_disc_mat, col_idxs = col_idxs, row_idxs = row_idxs)
    }
  }
})

########
# Test 4
########
test_that("Illegal subsets and extracts", {
  for (i in 1:n_datasets) {
    on_disc_mat <- cov_odms[[i]]@ondisc_matrix
    # index OOB
    expect_error(on_disc_mat[,ncol(on_disc_mat) + 10])
    expect_error(on_disc_mat[nrow(on_disc_mat) + 10,])
    # duplicate indexes
    expect_error(on_disc_mat[c(1,1),])
    expect_error(on_disc_mat[,c(1,1)])
    # extracting nothing
    expect_error(on_disc_mat[[,]])
  }
})

########
# Test 5
########
test_that("Test correct dimensions after subset", {
  for (i in 1:n_datasets) {
    Mat <- r_mats[[i]]
    on_disc_mat <- cov_odms[[i]]@ondisc_matrix
    subset_size_col <- sample(1:(ceiling(ncol(Mat)/30)), 1)
    subset_size_row <- sample(1:(ceiling(nrow(Mat)/30)), 1)
    col_idxs <- sample(x = 1:ncol(Mat), size = subset_size_col)
    row_idxs <- sample(x = 1:nrow(Mat), size = subset_size_row)

    t1 <- on_disc_mat[,col_idxs]
    t2 <- Mat[,col_idxs,drop=FALSE]
    expect_equal(dim(t1), dim(t2))

    t1 <- on_disc_mat[,-col_idxs]
    t2 <- Mat[,-col_idxs,drop=FALSE]
    expect_equal(dim(t1), dim(t2))

    t1 <- on_disc_mat[row_idxs,]
    t2 <- Mat[row_idxs,,drop=FALSE]
    expect_equal(dim(t1), dim(t2))

    t1 <- on_disc_mat[-row_idxs,]
    t2 <- Mat[-row_idxs,,drop=FALSE]
    expect_equal(dim(t1), dim(t2))

    t1 <- Mat[row_idxs, col_idxs,drop=FALSE]
    t2 <- on_disc_mat[row_idxs, col_idxs]
    expect_equal(dim(t1), dim(t2))

    t1 <- Mat[-row_idxs, -col_idxs,drop=FALSE]
    t2 <- on_disc_mat[-row_idxs, -col_idxs]
    expect_equal(dim(t1), dim(t2))
  }
})

########
# Test 6
########
test_that("Extract arbitrary submatrices after subset", {
  for (i in 1:n_datasets) {
    Mat <- r_mats[[i]]
    on_disc_mat <- cov_odms[[i]]@ondisc_matrix
    for (j in 1:n_reps) {
      subset_size_col <- sample(1:(ceiling(ncol(Mat)/10)), 1)
      subset_size_row <- sample(1:(ceiling(nrow(Mat)/10)), 1)
      col_idxs <- sample(x = 1:ncol(Mat), size = subset_size_col)
      row_idxs <- sample(x = 1:nrow(Mat), size = subset_size_row)
      # Perform subset
      # Matrix first
      Mat_row_sub <- Mat[row_idxs,,drop=FALSE]
      Mat_col_sub <- Mat[,col_idxs,drop=FALSE]
      Mat_sub <- Mat[row_idxs, col_idxs,drop=FALSE]
      # on_disc_matrix second
      on_disc_mat_row_sub <- on_disc_mat[row_idxs,]
      on_disc_mat_col_sub <- on_disc_mat[,col_idxs]
      on_disc_mat_sub <- on_disc_mat[row_idxs, col_idxs]
      # Next, generate row and column indexes to extract on for the subset matrices
      col_idxs_sub <- sample(x = 1:ncol(Mat_sub), size = sample(1:ncol(Mat_sub), 1))
      row_idxs_sub <- sample(x = 1:nrow(Mat_sub), size = sample(1:nrow(Mat_sub), 1))
      # Test extracts
      compare_Mat_on_disc_extract(Mat = Mat_row_sub, on_disc_mat = on_disc_mat_row_sub, col_idxs = col_idxs, row_idxs = row_idxs_sub) # Row subset
      compare_Mat_on_disc_extract(Mat = Mat_col_sub, on_disc_mat = on_disc_mat_col_sub, col_idxs = col_idxs_sub, row_idxs = row_idxs)
      compare_Mat_on_disc_extract(Mat = Mat_sub, on_disc_mat = on_disc_mat_sub, col_idxs = col_idxs_sub, row_idxs = row_idxs_sub)
    }
  }
})

########
# Test 7
########
test_that("Subset/extract corner cases", {
  for (i in 1:n_datasets) {
    Mat <- r_mats[[i]]
    on_disc_mat <- cov_odms[[i]]@ondisc_matrix
    # subset nothing
    on_dist_mat_sub <- on_disc_mat[]
    expect_identical(on_disc_mat, on_dist_mat_sub)
    # find a row of all zeros
    zero_rows <- which(Matrix::rowSums(m) == 0)
    if (length(zero_rows) >= 1) {
      idx <- zero_rows[1]
      zero_extract <- as.numeric(on_disc_mat[[idx,]])
      expect_equal(zero_extract, rep(0, ncol(m)))
    }
  }
})

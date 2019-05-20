args = commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
	print("ERROR: Provide two command-line arguments specifying which datasets to compare, and a third specifying size of header measured as number of lines.")
	quit('no')
}

options(digits=17)

# C=25
# C=20
# C=10
# C=5
# C=4
C=3
# C=2
# C=1
# eps=0.000000000000002
# largest_diff = 0.0

header_line_length = as.numeric(args[3])

use_custom_eps=FALSE
if (length(args) == 4) {
	use_custom_eps=TRUE
	eps <- as.numeric(args[4])
	print(paste("Comparing points with acceptable error of", eps))
}

diff_found <- FALSE

master <- read.table(paste0(args[1]), sep="", header=FALSE, skip=header_line_length)
test <- read.table(paste0(args[2]), sep="", header=FALSE, skip=header_line_length)
for (c in 1:ncol(master)) {
	master_col <-master[, c]
	test_col <- test[, c]

	if (length(master_col) != length(test_col)) {
		print(paste("Master has", length(master_col), "rows but test has", length(test_col), "rows."))
		diff_found = TRUE
		break
	}

	master_nans <- is.nan(master_col)
	if (sum(master_nans) > 0) {
		print(paste("Master contains",sum(master_nans),"/",length(master_col),"NaN's at:"))
		L <- min(C, sum(master_nans))
		indices <- (1:length(master_nans))[master_nans][1:L]
		print(indices[1:L]-1)
		diff_found = TRUE
		# break
	}

	test_nans <- is.nan(test_col)
	if (sum(test_nans) > 0) {
		print(paste("Test contains",sum(test_nans),"/",length(test_col),"NaN's at:"))
		L <- min(C, sum(test_nans))
		indices <- (1:length(test_nans))[test_nans][1:L]
		print(indices[1:L]-1)
		diff_found = TRUE
		# break
	}

	diff <- master_col-test_col

	# print(diff)
	# quit('no')

	# largest_diff = max(c(largest_diff, max(abs(diff))))
	# nan_indices = ((1:length(test_col))[is.na(diff)])
	# if (max(abs(diff)) > eps) {
	# 	print(paste("Column",c,"has a discrepancy:"))
	# 	print(paste("Largest", C, args[1], "differences:"))
	# 	print(diff[order(-abs(diff))][1:C])

	# 	print(paste("First", C, args[1], "differences:"))
	# 	print(diff[abs(diff)>eps][1:C])
	# 	print("Indices:")
	# 	print((1:nrow(master))[abs(diff)>eps][1:C] - 1)
	# 	print("")

	# 	first_idx = (1:nrow(master))[abs(diff)>eps][1]
	# 	print("First differing value is:")
	# 	print(test_col[first_idx])
	# 	print("..., it should be:")
	# 	print(master_col[first_idx])
	# 	break
	# }

	if (use_custom_eps) {
		max_allowed_error <- rep(eps, length(master_col))
	} else {
		# If floating-point operations have been reordered, then a difference
		# is expected due to rounding-errors, but the difference should
		# be smaller than the following:
		#   1 x 10 ^ ( E -17 + N )
		# Where E = exponent of master value
		#       N = largest difference in exponents of any floating-point
		#           arithmetic operation performed

		# N represents how many of the least-significant base-10 digits
		# of the floating-point mantissa are allowed to differ due to 
		# FP arithmetic reordering. Its value is guessed, as to set it
		# accurately would require a trace of all floating-point operation
		# outputs during the runs.
		N = 8

		max_allowed_error <- abs(master_col * (10^(-17+N)))

		# # Ignore any errors smaller than 3e-19:
		# max_allowed_error <- pmax(max_allowed_error, rep(3e-19, length(master_col)))

		# Ignore any errors smaller than 3e-17:
		max_allowed_error <- pmax(max_allowed_error, rep(3e-17, length(master_col)))
	}
	unexpected_errors_mask <- (abs(diff) > max_allowed_error) | (is.nan(test_col) & !is.nan(master_col))

	# false_errors <- unexpected_errors_mask & (master_col == test_col)

	# When there is no difference, R will be comparing 0.0 with a very small negative
	# value which equates to TRUE, so need to mask out rows which match exactly:
	# unexpected_errors_mask <- unexpected_errors_mask & !(master_col == test_col)
	matching_mask <- (master_col == test_col) & !is.nan(test_col) & !is.nan(master_col)
	unexpected_errors_mask <- unexpected_errors_mask & !matching_mask

	# print(unexpected_errors_mask[2])
	# print(abs(diff)[2] > max_allowed_error[2])
	# quit('no')

	unexpected_errors <- diff[unexpected_errors_mask]

	if (length(unexpected_errors) > 0) {
		print(paste(length(unexpected_errors), "DIFFS FOUND IN COLUMN", paste0(c, "/", ncol(master))))
		diff_found <- TRUE

		L <- min(C, sum(unexpected_errors_mask))
		indices <- (1:nrow(master))[unexpected_errors_mask]

		if (L < C) {
			print(paste(L, args[1], "differences:"))
		} else {
			print(paste("First", L, args[1], "differences:"))
		}
		print(unexpected_errors[1:L])
		print("0-indexed indices:")
		print(indices[1:L]-1)

		# range_start=indices[1]
		# range_end_idx=1
		# range_end=indices[1]
		# while (range_end_idx != (L-1)) {
		# 	if (indices[range_end_idx+1] == range_end+1) {
		# 		range_end = range_end+1
		# 		range_end_idx = range_end_idx+1
		# 	}
		# 	else {
		# 		if (range_end == range_start) {
		# 			print(range_start)
		# 		} else {
		# 			print(paste(range_start, "to", range_end))
		# 		}
		# 		range_end_idx = range_end_idx+1
		# 		range_start = indices[range_end_idx]
		# 		range_end = indices[range_end_idx]
		# 	}
		# }
		# print(paste(range_start, "to", range_end))

		print("Master values are:")
		print(master_col[indices[1:L]])

		print("Incorrect values are:")
		print(test_col[indices[1:L]])

		# print("All indices:")
		# print(indices)

		# print("Permitted error is:")
		# print(max_allowed_error[indices[1:L]])

		# break
	}
	# else {
	# 	print(paste("NO DIFFS FOUND IN COLUMN", c))
	# }

	if (diff_found) {
		break;
	}
}


# if (largest_diff > eps) {
if (diff_found == TRUE) {
	# Difference is too large to be attributed to reordering of 
	# floating-point arithmetic, so must be caused by incorrect code.
	file.create("mismatch.out")
}

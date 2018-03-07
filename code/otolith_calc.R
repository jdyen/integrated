# load otolith data and convert it to probabilities for a
#   prior in the inverse population model
otolith_prior_fun <- function(oti_data, nstage) {
  
  nbreaks <- nstage + 1
  
  # calculate size_now and size_next
  size_now <- NULL
  size_next <- NULL
  id <- NULL
  sizemax <- max(tapply(oti_data$growth, oti_data$id, sum, na.rm = TRUE),
                 na.rm = TRUE)
  for (i in seq_along(unique(oti_data$id))) {
    oti_sub <- oti_data[which(oti_data$id == unique(oti_data$id)[i]), ]
    size_tmp <- cumsum(oti_sub$growth[order(oti_sub$year)])
    for (j in seq_len((length(size_tmp) - 1))) {
      size_now <- c(size_now, size_tmp[j])
      size_next <- c(size_next, size_tmp[j + 1])
      id <- c(id, unique(oti_sub$id)[1])
    }
  }
  oti_clean <- data.frame(id = id,
                          size_now = (size_now / sizemax),
                          size_next = (size_next / sizemax))
  
  # calculate breaks
  break_set <- c(0, quantile(oti_clean$size_now,
                             p = seq(0.1, 0.9, length = (nbreaks - 2))), 1)

  label_set <- seq_len(length(break_set) - 1)
  oti_clean$bin_now <- as.numeric(cut(oti_clean$size_now,
                                      breaks = break_set,
                                      labels = label_set))
  oti_clean$bin_next <- as.numeric(cut(oti_clean$size_next,
                                       breaks = break_set,
                                       labels = label_set))

  # calculate transition probabilities  
  transition_mat <- matrix(0, nrow = (nbreaks - 1), ncol = (nbreaks - 1))
  for (i in seq_len(nrow(oti_clean))) {
    xind <- oti_clean$bin_next[i]
    yind <- oti_clean$bin_now[i]
    transition_mat[xind, yind] <- transition_mat[xind, yind] + 1
  }
  
  transition_mat
  
}

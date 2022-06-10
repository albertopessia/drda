# log-logistic test data
lltd <- list(
  D = data.frame(
    x = rep(
      c(0, 2, 4, 6, 8, 10, 100, 1000), times = c(3, 3, 2, 4, 3, 1, 2, 1)
    ),
    y = c(
      0.8527, 0.7571, 0.937, 0.8275, 0.7778, 0.8846, 0.5561, 0.7048, 0.5453,
      0.486, 0.5351, 0.4399, 0.3553, 0.3132, 0.3499, 0.1421, 0.0168, 0.1386,
      0.1227
    ),
    w = c(
      0.718, 0, 0.5804, 1.4958, 0.8458, 1.1844, 0, 1.4749, 0, 1.332, 1.039,
      0.6865, 0.6791, 1.1611, 0.7153, 1.3647, 0.6719, 0.4121, 1.1564
    )
  ),
  theta_6 = c(0.8, -1, 2, 2, 4, 3),
  theta_5 = c(1.0, -0.85, 2, 5, 2),
  theta_4 = c(1.0, -0.85, 2, 5),
  theta_2_d = c(1.0, -1.0, 2, 5),
  theta_2_i = c(0.0, 1.0, 2, 5),
  theta_g = c(1.0, -0.85, 2, 5),
  sigma = 0.05
)

lltd$stats_1 <- matrix(
  c(
    sort(unique(lltd$D$x)),
    as.numeric(table(lltd$D$x)),
    by(lltd$D$y, lltd$D$x, mean),
    by(lltd$D$y, lltd$D$x, function(z) sum((z - mean(z))^2) / length(z))
  ),
  nrow = length(unique(lltd$D$x)),
  ncol = 4,
  dimnames = list(NULL, c("x", "n", "m", "v"))
)

lltd$stats_1_i <- matrix(
  c(
    sort(unique(lltd$D$x)),
    as.numeric(table(lltd$D$x)),
    by(rev(lltd$D$y), lltd$D$x, mean),
    by(rev(lltd$D$y), lltd$D$x, function(z) sum((z - mean(z))^2) / length(z))
  ),
  nrow = length(unique(lltd$D$x)),
  ncol = 4,
  dimnames = list(NULL, c("x", "n", "m", "v"))
)

lltd$stats_2 <- matrix(
  c(
    sort(unique(lltd$D$x)),
    by(lltd$D, lltd$D$x, function(z) sum(z$w)),
    by(lltd$D, lltd$D$x, function(z) sum(z$w * z$y) / sum(z$w)),
    by(lltd$D, lltd$D$x, function(z) {
      sum(z$w * (z$y - sum(z$w * z$y) / sum(z$w))^2) / sum(z$w)
    })
  ),
  nrow = length(unique(lltd$D$x)),
  ncol = 4,
  dimnames = list(NULL, c("x", "n", "m", "v"))
)

lltd$stats_2_i <- matrix(
  unlist(
    by(
      cbind(x = lltd$D$x, y = rev(lltd$D$y), w = lltd$D$w),
      lltd$D$x,
      function(z) {
        c(
          z$x[1], sum(z$w), sum(z$w * z$y) / sum(z$w),
          sum(z$w * (z$y - sum(z$w * z$y) / sum(z$w))^2) / sum(z$w)
        )
      }
    )
  ),
  nrow = length(unique(lltd$D$x)),
  ncol = 4,
  byrow = TRUE,
  dimnames = list(NULL, c("x", "n", "m", "v"))
)

# logistic test data
ltd <- list(
  D = data.frame(
    x = rep(
      c(-100, -30, -6, -2, 2, 6, 50, 100),
      times = c(3, 3, 2, 4, 3, 1, 2, 1)
    ),
    y = c(
      0.8527, 0.7571, 0.937, 0.8275, 0.7778, 0.8846, 0.5561, 0.7048, 0.5453,
      0.486, 0.5351, 0.4399, 0.3553, 0.3132, 0.3499, 0.1421, 0.0168, 0.1386,
      0.1227
    ),
    w = c(
      0.718, 0, 0.5804, 1.4958, 0.8458, 1.1844, 0, 1.4749, 0, 1.332, 1.039,
      0.6865, 0.6791, 1.1611, 0.7153, 1.3647, 0.6719, 0.4121, 1.1564
    )
  ),
  theta_6 = c(0.8, -1, 2, log(2), 4, 3),
  theta_5 = c(0.9, -0.7, 0.1, -4, 2),
  theta_4 = c(0.9, -0.7, 0.1, -4),
  theta_2_d = c(1.0, -1.0, 0.1, -4),
  theta_2_i = c(0.0, 1.0, 0.1, -4),
  theta_g = c(0.9, -0.7, 0.1, -4),
  sigma = 0.05
)

ltd$stats_1 <- matrix(
  c(
    sort(unique(ltd$D$x)),
    as.numeric(table(ltd$D$x)),
    by(ltd$D$y, ltd$D$x, mean),
    by(ltd$D$y, ltd$D$x, function(z) sum((z - mean(z))^2) / length(z))
  ),
  nrow = length(unique(ltd$D$x)),
  ncol = 4,
  dimnames = list(NULL, c("x", "n", "m", "v"))
)

ltd$stats_1_i <- matrix(
  c(
    sort(unique(ltd$D$x)),
    as.numeric(table(ltd$D$x)),
    by(rev(ltd$D$y), ltd$D$x, mean),
    by(rev(ltd$D$y), ltd$D$x, function(z) sum((z - mean(z))^2) / length(z))
  ),
  nrow = length(unique(ltd$D$x)),
  ncol = 4,
  dimnames = list(NULL, c("x", "n", "m", "v"))
)

ltd$stats_2 <- matrix(
  c(
    sort(unique(ltd$D$x)),
    by(ltd$D, ltd$D$x, function(z) sum(z$w)),
    by(ltd$D, ltd$D$x, function(z) sum(z$w * z$y) / sum(z$w)),
    by(ltd$D, ltd$D$x, function(z) {
      sum(z$w * (z$y - sum(z$w * z$y) / sum(z$w))^2) / sum(z$w)
    })
  ),
  nrow = length(unique(ltd$D$x)),
  ncol = 4,
  dimnames = list(NULL, c("x", "n", "m", "v"))
)

ltd$stats_2_i <- matrix(
  unlist(
    by(
      cbind(x = ltd$D$x, y = rev(ltd$D$y), w = ltd$D$w),
      ltd$D$x,
      function(z) {
        c(
          z$x[1], sum(z$w), sum(z$w * z$y) / sum(z$w),
          sum(z$w * (z$y - sum(z$w * z$y) / sum(z$w))^2) / sum(z$w)
        )
      }
    )
  ),
  nrow = length(unique(ltd$D$x)),
  ncol = 4,
  byrow = TRUE,
  dimnames = list(NULL, c("x", "n", "m", "v"))
)

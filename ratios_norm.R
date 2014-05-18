standardize_ratios <- function(ratios) {
  t(scale(t(ratios),center = apply(ratios, 1, median, na.rm = T), scale = apply(ratios, 1, sd, na.rm = T)))
}


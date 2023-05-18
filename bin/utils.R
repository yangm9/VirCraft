#
generate_named_vector <- function(keys, values) {
  n <- length(keys)
  named_vector <- setNames(values[1:n], keys)
  return(named_vector)
}

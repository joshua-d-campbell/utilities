interleave = function(a, b) {
  c(a,b)[ order( c(seq_along(a), seq_along(b)))]
}  
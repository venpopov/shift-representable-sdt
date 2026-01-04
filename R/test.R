library(profvis)

prof_out = profvis({
  fitSR(list(y1$total,y1$total))
})
# Original script by John Dunn, retrieved from https://osf.io/n62y4/files/osfstorage
# Minor edits for my workflow by Ven Popov

readexp <- function(expno = NULL, dir = NULL) {
  file <- sprintf("Experiment%d.csv", expno)
  filepath <- paste0(dir, file)
  exp <- read.csv(filepath)
  y <- list()
  n <- length(unique(exp$Participant))
  for (i in 1:n) {
    s <- exp[exp$Participant == i, ]
    y[[i]] <- as.matrix(s[, 3:8])
  }
  total <- Reduce("+", y)
  output <- NULL
  output$participants <- y
  output$total <- total
  output$mean <- total / n
  return(output)
}

# shift-representable-sdt

This repository contains optimized code for the analyses presented in Dunn, J. C., Anderson, L. M., & Stephens, R. G. (2025). A shift representable signal detection model of recognition memory. The original code was written by John Dunn and was obtained from https://osf.io/n62y4/. The current version provides a 10x speed-up, but is otherwise identical.

To use, first clone the repository:

```bash
git clone https://github.com/venpopov/shift-representable-sdt.git
```

Before running the code, you might need to install some packages. I use `renv` for package management. `renv` should install itself due to code in the .Rprofile file when you open the project in Rstudio or open an R console in the main directory. Afterwards run this line in the R console to install any missing dependencies:

```R
renv::restore()
```

Then you can follow John's tutorial from the `fitSR tutorial.pdf` file
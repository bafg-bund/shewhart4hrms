

# Create new directory and install all dependencies and shewhart4hrms
tempDir <- withr::local_tempdir()

oldWd <- getwd()

renv::init(project = tempDir, bare = T, restart = F)
renv::restore(project = tempDir, lockfile = file.path(oldWd, "renv.lock"), exclude = "ntsworkflow", prompt = F)

# install ntsworkflow
renv::install(file.path(oldWd, "tests", "testthat", "fixtures", "ntsworkflow_0.2.5.zip"), project = tempDir, prompt = F)

# install shewhart4hrms
renv::install(file.path(oldWd, "tests", "testthat", "fixtures", "shewhart4hrms_0.2.zip"), project = tempDir, prompt = F)

# test package loading
if (suppressWarnings(suppressMessages(require(shewhart4hrms))))
  message("Test passed")

detach("package:shewhart4hrms", unload = T)
setwd(oldWd)
renv::activate(oldWd)

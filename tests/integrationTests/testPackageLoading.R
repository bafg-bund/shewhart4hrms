
if (suppressMessages(suppressWarnings(require(xcms))))
  message("test passed") else stop("test failed")

detach("package:xcms", unload = T)

if (suppressMessages(suppressWarnings(require(ntsworkflow))))
  message("test passed") else stop("test failed") 

detach("package:ntsworkflow", unload = T)


con <- file("test.log")
sink(con, append=TRUE, type="message")

setwd(dirnam(rstudioapi::getActiveDocumentContext()$path))

# This will echo all input and not truncate 150+ character lines...
source("package_builder.R", echo=TRUE, max.deparse.length=10000)

# Restore output to console
sink(type="message")


# And look at the log...
cat(readLines("test.log"), sep="\n")


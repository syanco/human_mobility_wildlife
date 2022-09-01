'
Calculate space use before and during COVID-19 lockdowns using previously 
estimated dBBMMs

Usage:
calc-space-use.r [-c]
calc-space-use.r (-h | --help)

Parameters:

  
Options:
-h --help     Show this screen.
-v --version     Show version.
-c --continue   Indicates if script should first check for previous entries before calculating metrics
' -> doc
library(docopt)
ag <- docopt(doc, version = '0.1\n')

.continue <- as.logical(ag$continue)

message(.continue)

if(.continue){
  message("yeah, it worked")
}

if(.continue == T){
  message("yeah, that worked too")
}
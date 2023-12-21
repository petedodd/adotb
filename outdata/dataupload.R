library(here)
library(data.table)
library(glue)
library(googlesheets4)


## NOTE authors only
## create an ID to access the googlesheets results sheet
yourl <- "https://docs.google.com/spreadsheets/d/1IyqshINvFHWy1Gqk5NV__W2KKNX-TJPnYpbq6KJSu0o/edit#gid=0"
shid <- as.character(as_sheets_id(yourl))


## utility function
upload.to.sheets <- function(basename,filename,sheetid
                             ){
  filename <- gsub("\\.csv$","",filename) #safety in case csv included at and
  fn <- glue(basename) + filename + ".csv"
  tmp <- fread(file=fn)
  write_sheet(tmp,sheetid,sheet=filename)
}

## upload relevant table data
upload.to.sheets(here('outdata/'),'cftab',shid) #first will need to re-authenticate

## rest can be run as block
flz1 <- c(
  'CDR.csv',
  'cdro.csv',
  'cdry.csv',
  ## 'cftab.csv',
  'cfWS.csv',
  'ltbi.crange.csv',
  'ltbi.PC.age.csv',
  'ltbi.PC.noage.csv',
  'ltbi.tot.age.csv',
  'ltbi.tot.noage.csv',
  'nhiv.csv',
  'out.inc.csv',
  'out.ltbi.csv',
  'PAF.csv',
  'PC_ranked_1014.csv',
  'PC_ranked_1519.csv',
  'RRstats.csv',
  'smy.csv',
  'thinness.RRs.csv'
)

for( fn in flz1)
  upload.to.sheets(here('outdata/'),fn,shid)

Sys.sleep(120) #wait a bit so as not to annoy google

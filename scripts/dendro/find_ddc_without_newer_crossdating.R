## Find DDC cores that were 
#    (a) not reviewed by Paige
#    (b) not crossdated by someone else

setwd("~/UC Davis/Research Projects/Sierra dendro/sierra-dendro")

library(tidyverse)
library(lubridate)

crossdating = read.csv("data/dendro/crossdating-records/Crossdating progress records (Responses).csv")

crossdating$done_date = mdy_hms(crossdating$Timestamp)
crossdating$Core.number = toupper(crossdating$Core.number)

cd_ddc = crossdating %>%
  filter(Derek.double.check.needed. != "")

cd_no_ddc = crossdating %>%
  filter(Derek.double.check.needed. == "")


### get the most recent date of the DDC cores
cd_ddc = cd_ddc %>%
  group_by(Core.number) %>%
  summarize(done_date = max(date))



# for each DDC core, was there a more recent no_ddc record?

ddc_w_more_recent = NULL

for(i in 1:nrow(cd_ddc)) {
  
  # get its name and date
  core = cd_ddc[i,]$Core.number
  focal_done_date = cd_ddc[i,]$done_date
  
  # filter non-ddc cores to be after that date, and have the same core number
  cd_no_ddc_after = cd_no_ddc %>%
    filter(done_date > focal_done_date) %>%
    filter(Core.number == core)
  
  if(nrow(cd_no_ddc_after) > 0) {
    ddc_w_more_recent = c(ddc_w_more_recent,core)
  }
  
  
}
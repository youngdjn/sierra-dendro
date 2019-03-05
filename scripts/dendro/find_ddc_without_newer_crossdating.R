## Find DDC cores that were 
#    (a) not reviewed by Paige
#    (b) not crossdated by someone else

setwd("~/UC Davis/Research Projects/Sierra dendro/sierra-dendro")

library(tidyverse)
library(lubridate)
library(readxl)

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
  dplyr::summarize(done_date = max(done_date))

ddc_cores = cd_ddc$Core.number %>% toupper()



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

ddc_w_more_recent = ddc_w_more_recent %>% toupper()

## now ddc_w_more_recent has list of cores that have been re-crossdated without a DDC request

## now we need to look up all the cores that Cameron and Paige looked over
already_done = read_excel("data/dendro/crossdating-records/core_correlations_with_ref_chron.xlsx")

already_done = already_done %>%
  filter(!is.na(comment))

files_done = already_done$filename

files_done = str_replace(files_done,"t.pos","")
cores_done = str_replace(files_done,".pos","")


ddc_reduced = setdiff(ddc_cores,ddc_w_more_recent)
ddc_reduced2 = setdiff(ddc_reduced,cores_done)


### get the DDC comments for the remaining DDC cores

ddc_cores

remaining_ddc = crossdating %>%
  filter(Derek.double.check.needed. != "") %>%
  filter(done_date %in% cd_ddc$done_date) %>%
  filter(Core.number %in% ddc_reduced2)






## now we need to see if there were any very recent crossdating records requesting a DDC (in case this was recorded when checking the poor-correlation cores)

crossdating_recent = crossdating %>%
  filter(done_date > ymd("2015-12-30"))

remaining_ddc = remaining_ddc %>%
  filter(!(Core.number %in% crossdating_recent$Core.number))



write.csv(remaining_ddc,"data/dendro/crossdating-records/remaining_derek_double_check_cores.csv")

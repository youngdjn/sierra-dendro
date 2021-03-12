d = read.csv("data/compiled-for-analysis/trees.csv")
y = read.csv("data/compiled-for-analysis/years.csv")

## how many have an na between 2000 and 2013?

y_summ = y %>%
  group_by(year) %>%
  summarize(depth = sum(!is.na(rwi))) %>%
  ungroup() %>%
  filter(between(year,2000,2013))

write.csv(y_summ,"samp_depth.csv")

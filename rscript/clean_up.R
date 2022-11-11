# clean script 50 cores
library(here)
library(purrr)
library(furrr)
library(glue)
library(stringr)
library(parallel)



mapsets1 <- list.files('/home/tcrnbgh/Scratch/tmp/grassdb',
           pattern='mapset*',
           full.names = T)

numPerCore <- ceiling(length(mapsets1) / 20)

mapsetList <- 
  seq(1,length(mapsets1),by=numPerCore) %>% 
  map(~.x:((.x-1)+numPerCore))

mapsetList[[length(mapsetList)]] <- 
  mapsetList[[length(mapsetList)]][which(mapsetList[[length(mapsetList)]] <= length(mapsets1))]

mapsetList.path <- mapsetList %>% 
  map(~mapsets1[unlist(.x)])

cl <- parallel::makePSOCKcluster(20)
datOut <- parallel::clusterApply(cl, mapsetList.path, function(x) {
  library(glue)
  for (p in x) {
    system(glue('rsync -a --delete /home/tcrnbgh/Scratch/empty_dir/ {p}/'))
  }
})

future::plan('multisession',workers=20)
future.apply::future_lapply(mapsetList.path, function(x) {
  library(glue)
  for (p in x) {
    system(glue('rsync -a --delete /home/tcrnbgh/Scratch/empty_dir/ {p}/'))
  }
})



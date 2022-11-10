


timeStamps <- read.delim('logs/viz_out_o3.txt') %>% 
  as.data.frame() %>% 
  .[10:nrow(.),] %>% 
  .[seq(1,length(.),by=2)] %>%
  lubridate::as_datetime()

1:length(timeStamps) %>% map(~difftime(timeStamps[.x+1],timeStamps[.x])) %>% 
  as.numeric() %>% 
  mean(.,na.rm=T)

120*60  / 1000

ln <- stringr::str_split("2022-11-04 16:24:25.555726
point number 622393...
2022-11-04 16:24:37.999578
point number 701252...
2022-11-04 16:24:39.732596
point number 361259...
2022-11-04 16:24:50.248361
point number 392621...
2022-11-04 16:24:55.012520
point number 756587...
2022-11-04 16:25:02.636655
point number 171079...
2022-11-04 16:25:10.301242
point number 153713...
2022-11-04 16:25:15.138475
point number 749477...
2022-11-04 16:25:25.815576
point number 429901...
2022-11-04 16:25:27.726034
point number 539612...
2022-11-04 16:25:39.973502
point number 290622...
2022-11-04 16:25:41.203988
point number 522959...
2022-11-04 16:25:52.380951
point number 483287...
2022-11-04 16:25:56.679913
point number 684077...
2022-11-04 16:26:04.926240
point number 26451...
2022-11-04 16:26:11.697741",
                         '\n')
timeStamps <- unlist(ln)[seq(1,31,2)] %>% 
  lubridate::as_datetime()

diffs <- 1:length(timeStamps) %>% 
  map(~difftime(timeStamps[.x+1],timeStamps[.x]))

mean(as.numeric(diffs),na.rm=T)

zoo::rollapply(timeStamps,FUN=difftime,width=2)
timeStamps[1]-timeStamps[2]

cvaTimes <- list.files('logs',pattern='analysis_cva',full.names = T) %>% 
  map_chr(~tail(readLines(.x),n=1))

vsTimes <- list.files('logs',pattern='analysis_sinkout*',full.names = T) %>% 
  map_chr(~tail(readLines(.x),n=1))

data.frame(cvaTimes,vsTimes)

((3 * 40000) / 60) / 60

under20 <- r.tiles.length < 20000
sum(under20)
923-233


secsPerPoint <- (6.8 / 1000) * 60

outRas <- list.files('outputs',full.names = T) %>% 
  map(rast)

rMerged <- do.call('merge', outRas)
plot(rMerged)
writeRaster(rMerged,'merged.tif')






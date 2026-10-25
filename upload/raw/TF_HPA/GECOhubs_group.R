rm(list=ls())

cts = c('alpha','CD4NC','CD4ET','CD8NC','CD8ET','NK','BIN','BMem')
hubs_set = list()
for (i in 2:8) {
  hubs_set[[i]] = readLines(paste0('data/onek1k/hubs_', cts[i], '.txt'))
}
hubs_set[[1]] = readLines('data/onek1k/GECO_Shared.txt')
names(hubs_set) = c('Shared','CD4NC','CD4ET','CD8NC','CD8ET','NK','BIN','BMem')

hubs_group = list('Shared' = hubs_set[['Shared']],
                  'Tcell'  = unique(unlist(hubs_set[c('CD4NC','CD4ET','CD8NC','CD8ET')])),
                  'NKcell' = hubs_set[['NK']],
                  'Bcell'  = unique(unlist(hubs_set[c('BIN','BMem')])))

for (n in names(hubs_group)) {
  write.table(hubs_group[[n]], file.path('data/4a-TF_HPA/GECOhubs_group', paste0(n, '.txt')),
              row.names = FALSE, col.names = FALSE, sep = ' ', quote = FALSE)
}

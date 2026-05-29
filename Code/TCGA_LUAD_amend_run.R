library(igraph)
library(AMEND)

path.results <- "Results/"
checkpoint.file <- paste0(path.results, "LUAD_net_checkpoint.rds")
chk <- readRDS(checkpoint.file)
input_graphs <- chk$input_graphs
data.list <- chk$data.list
brw.list <- chk$brw.list

N <- 50
jp <- c(mirna=1, mrna=0.5, methyl=1)
sp <- list(mrna= c(string=0.5, path=0.5,aracne=0.5,ohmnet=0.5))
nw <- c(mirna=0.4, mrna=0.15, methyl=0.45)
lw <- list(mrna=c(string=0.25, path=0.15, aracne=0.3, ohmnet=0.3))
seed.function <- 'exp'
fun.params <- list(DOI=0)

#run 1 with BRW
message(paste("Run 1 with BRW : started at", Sys.time()))
subnet.brw <- run_AMEND(graph=input_graphs, n=N, data=data.list, brw.attr=brw.list,
                        aggregate.multiplex=NULL, degree.bias=NULL,
                        multiplex=TRUE, heterogeneous=TRUE, normalize="degree",
                        FUN=seed.function, FUN.params=fun.params,
                        jump.prob=jp, net.weight=nw, switch.layer.prob=sp,
                        layer.weight=lw, verbose=TRUE, in.parallel=FALSE)
saveRDS(subnet.brw, file=paste0(path.results, "LUAD_n", N, "_db.none_brw.rds"))

#run 2 without BRW/IN/aggr
message(paste("Run 2 without BRW/IN/aggr : started at", Sys.time()))
subnet.nobrw <- run_AMEND(graph = input_graphs, n=N, data=data.list, brw.attr=NULL,
                           aggregate.multiplex=NULL, degree.bias=NULL,
                           multiplex=TRUE, heterogeneous=TRUE, normalize="degree",
                           FUN=seed.function, FUN.params=fun.params,
                           jump.prob=jp, net.weight=nw, switch.layer.prob=sp,
                           layer.weight=lw, verbose=TRUE, in.parallel=FALSE)
saveRDS(subnet.nobrw, file=paste0(path.results, "LUAD_n", N, "_db.none.rds"))

#run 3 with BRW + IN degree bias + aggregate.multiplex (original analysis)
message(paste("Run 3: BRW + IN + aggregate.multiplex : started at", Sys.time()))
subnet.main <- run_AMEND(graph=input_graphs, n=N, data=data.list, brw.attr=brw.list,
                          aggregate.multiplex= list(primary="mrna_string", agg.method="mean"),
                          degree.bias= list(method="IN", component="mrna_string"),
                          multiplex=TRUE, heterogeneous=TRUE, normalize="degree",
                          FUN=seed.function, FUN.params=fun.params,
                          jump.prob=jp, net.weight=nw, switch.layer.prob=sp,
                          layer.weight=lw, verbose=TRUE, in.parallel=FALSE)
saveRDS(subnet.main, file=paste0(path.results, "LUAD_n", N, "_IN_agg_brw.rds"))
message(paste("All done at", Sys.time()))

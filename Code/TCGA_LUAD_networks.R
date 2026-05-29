#molecular interaction networks for TCGA-LUAD

rm(list= ls())

library(igraph)
library(readxl)
library(dplyr)
library(graphite)
library(aracne.networks)
library(org.Hs.eg.db)
library(AnnotationDbi)

extract_string <- function(x, k, pos) unlist(lapply(strsplit(x, k), function(y) y[pos]))
largest_connected_component <- function(g){
  if(!igraph::is_connected(g)){
    message("Taking largest connected component.")
    comps <- igraph::components(g)
    largest_comp_id <- which.max(comps$csize)
    g <- igraph::induced_subgraph(g, which(comps$membership == largest_comp_id))
  }
  g
}

set.seed(2707083)

ppi_edge_threshold <- 0.7

path.networks <- "Data/Interaction Networks/"
path.results <- "Results/"
dir.create(path.results, showWarnings= FALSE)




#load data
omics.dat <- readRDS(paste0(path.results, "LUAD_multi_omics.rds"))
de.res <- readRDS(paste0(path.results, "LUAD_DE_results.rds"))
hr <- readRDS(paste0(path.results, "LUAD_hazard_ratios.rds"))

rownames(de.res$mirna) <- extract_string(rownames(de.res$mirna), "\\|", 1)
names(hr$mirna) <- extract_string(names(hr$mirna), "\\|", 1)

# fix MIR AS naming in methylation
id <- which(grepl('AS', rownames(de.res$methyl)) & grepl('MIR', rownames(de.res$methyl)))
for(i in id){
  tips <- extract_string(rownames(de.res$methyl)[i], '-', 1)
  ends <- extract_string(rownames(de.res$methyl)[i], '-', 2)
  tips <- paste0(tips, "B")
  ends <- gsub('AS', '', ends)
  rownames(de.res$methyl)[i] <- paste0(tips, ends)
}
id <- which(grepl('AS', names(hr$methyl)) & grepl('MIR', names(hr$methyl)))
for(i in id){
  tips <- extract_string(names(hr$methyl)[i], '-', 1)
  ends <- extract_string(names(hr$methyl)[i], '-', 2)
  tips <- paste0(tips, "B")
  ends <- gsub('AS', '', ends)
  names(hr$methyl)[i] <- paste0(tips, ends)
}




#mirna-mrna bipartite (miRTarBase)
dat <- read_xlsx(path= paste0(path.networks, "miRTarBase_SE_WR.xlsx")) %>%
  dplyr::filter(substr(miRNA, 1, 3)== "hsa" &
                  `Species (miRNA)`== "Homo sapiens" &
                  `Species (Target Gene)`== "Homo sapiens") %>%
  dplyr::rename(ID=`miRTarBase ID`,
                Symbol = `Target Gene`,
                EntrezID=`Target Gene (Entrez ID)`) %>%
  dplyr::select(ID, miRNA, Symbol, EntrezID) %>%
  dplyr::distinct()

mirna.el <- as.matrix(dat[, c('miRNA', 'EntrezID')])
mirna.el[,2] <- gsub(pattern= " ", replacement= "", x= mirna.el[,2])
bp.mirna.mrna <- igraph::graph_from_edgelist(el= mirna.el, directed= FALSE)
V(bp.mirna.mrna)$node_type <- ifelse(V(bp.mirna.mrna)$name %in% mirna.el[,1], "mirna", "mrna")
E(bp.mirna.mrna)$weight <- ppi_edge_threshold

vertex_attr(bp.mirna.mrna, "native_nom") <- ifelse(V(bp.mirna.mrna)$node_type == "mirna", "miRID", "EntrezID")
id <- which(V(bp.mirna.mrna)$node_type=="mirna")
vertex_attr(bp.mirna.mrna, "miRID", id) <- V(bp.mirna.mrna)$name[id]
id <- which(V(bp.mirna.mrna)$node_type == "mrna")
vertex_attr(bp.mirna.mrna, "EntrezID", id) <- V(bp.mirna.mrna)$name[id]




#ohmnet lung ppi
k.ppi <- data.table::fread(file= paste0(path.networks, "PPT-Ohmnet_tissues-combined.txt")) %>%
  dplyr::filter(tissue=="lung") %>%
  dplyr::select(-tissue) %>%
  as.matrix()
k.ppi <- matrix(as.character(k.ppi), ncol= 2, dimnames= dimnames(k.ppi))
k.ppi <- cbind(k.ppi, as.character(ppi_edge_threshold))
g.ohmnet <- igraph::graph_from_edgelist(el= k.ppi[,1:2], directed= FALSE)
E(g.ohmnet)$weight <- as.numeric(k.ppi[,3])
vertex_attr(g.ohmnet,"native_nom") <- "EntrezID"
vertex_attr(g.ohmnet, "EntrezID") <- V(g.ohmnet)$name



#reactome (graphite)
hsR <- graphite::pathways("hsapiens", "reactome")
print(length(hsR))
for(i in seq_along(hsR)){
  if(i %% 100 == 0)
    message(i)
  g.tmp <- graphite::pathwayGraph(pathway= hsR[[i]], which= "proteins")
  g.tmp <- igraph::graph_from_graphnel(graphNEL= g.tmp)
  g.tmp <- igraph::as.undirected(graph= g.tmp, mode= "collapse")
  if(i == 1) g.path <- g.tmp
  else g.path <- igraph::union(g.path, g.tmp)
}
edge.weight.attrs <- igraph::edge_attr_names(g.path)[grepl("weight", igraph::edge_attr_names(g.path))]
for(i in seq_along(edge.weight.attrs)){
  g.path <- igraph::delete_edge_attr(graph= g.path, name= edge.weight.attrs[i])
}
E(g.path)$weight <- ppi_edge_threshold
g.path <- igraph::delete_vertices(graph= g.path, v= which(!grepl('UNIPROT:', V(g.path)$name)))
V(g.path)$name <- gsub("UNIPROT:", "", V(g.path)$name)
vertex_attr(g.path, "native_nom") <- "UniProt"
vertex_attr(g.path, "UniProt") <- V(g.path)$name




#aracne grn (TCGA-LUAD)
if(1){
  data(regulonluad)
  el <- do.call(rbind, lapply(names(regulonluad), function(tf){
    reg <- regulonluad[[tf]]
    data.frame(Regulator=tf,
               Target = names(reg$tfmode),
               MoA=unname(reg$tfmode),
               likelihood = unname(reg$likelihood),
               stringsAsFactors= FALSE)
  }))
  write.table(el, file= paste0(path.networks, "aracne_el_luad.txt"),
              sep= "\t", row.names= FALSE, quote= FALSE)
}
grn <- data.table::fread(file= paste0(path.networks, "aracne_el_luad.txt"), sep= "\t", header= TRUE)
grn <- grn[grn$likelihood >= ppi_edge_threshold, -3]
grn <- as.matrix(grn)
grn <- matrix(as.character(grn), ncol= 3, dimnames= dimnames(grn))
g.aracne <- igraph::graph_from_edgelist(el= grn[,1:2], directed= FALSE)
E(g.aracne)$weight <- as.numeric(grn[,3])
vertex_attr(g.aracne, "native_nom") <- "EntrezID"
vertex_attr(g.aracne, "EntrezID") <- V(g.aracne)$name



#string ppi v12
ppi <- data.table::fread(file= paste0(path.networks, "9606.protein.physical.links.v12.0.txt.gz"), header= TRUE) %>%
  dplyr::mutate(combined_score= combined_score / 1000) %>%
  dplyr::filter(combined_score >= ppi_edge_threshold) %>%
  dplyr::mutate(protein1= extract_string(protein1, "\\.", 2),
                protein2= extract_string(protein2, "\\.", 2)) %>%
  dplyr::select(protein1, protein2, combined_score) %>%
  dplyr::distinct() %>%
  as.matrix()

tmp <- character(nrow(ppi))
for(i in seq_along(tmp)){
  tmp[i] <- paste(sort(c(ppi[i,1], ppi[i,2])), collapse= "|")
}
scores <- as.numeric(ppi[,3])
dat <- stats::aggregate(scores, by= list(tmp), FUN= min)
string.ppi <- matrix(c(extract_string(dat[,1], "\\|", 1),
                        extract_string(dat[,1], "\\|", 2),
                        dat[,2]), ncol= 3)
g.string <- igraph::graph_from_edgelist(el= string.ppi[,1:2], directed= FALSE)
E(g.string)$weight <- as.numeric(string.ppi[,3])
vertex_attr(g.string, "native_nom") <- "ENSPID"
vertex_attr(g.string, "ENSPID") <- V(g.string)$name

#id mapping (org.Hs.eg.db, offline)
uniq.ensp <- unique(c(string.ppi[,1], string.ppi[,2]))
ensp.vec <- AnnotationDbi::mapIds(org.Hs.eg.db, keys=uniq.ensp,
                                  column="ENTREZID", keytype="ENSEMBLPROT",
                                  multiVals="first")
ensp.vec <- ensp.vec[!is.na(ensp.vec)]
ensp2entrez <- data.frame(ensembl_peptide_id=names(ensp.vec),
                          entrezgene_id = as.character(ensp.vec),
                          stringsAsFactors = FALSE)

uniq.uniprot <- unique(V(g.path)$name)
uniprot.vec <- AnnotationDbi::mapIds(org.Hs.eg.db, keys=uniq.uniprot,
                                     column = "ENTREZID", keytype="UNIPROT",
                                     multiVals = "first")
uniprot.vec <- uniprot.vec[!is.na(uniprot.vec)]
uniprot2entrez <- data.frame(uniprotswissprot=names(uniprot.vec),
                             entrezgene_id = as.character(uniprot.vec),
                             stringsAsFactors = FALSE)




#apply mappings
id <- match(V(g.string)$name, ensp2entrez[,1])
V(g.string)$name[!is.na(id)] <- make.unique(ensp2entrez[id[!is.na(id)], 2], sep=".")
V(g.string)$EntrezID[!is.na(id)] <- ensp2entrez[id[!is.na(id)], 2]

id <- match(V(g.path)$name, uniprot2entrez[,1])
V(g.path)$name[!is.na(id)] <- make.unique(uniprot2entrez[id[!is.na(id)], 2], sep=".")
V(g.path)$EntrezID[!is.na(id)] <- uniprot2entrez[id[!is.na(id)], 2]




#methyl-mrna bipartite
uniq.symbols <- unique(c(rownames(de.res$mrna), rownames(de.res$methyl)))
symbol.vec <- AnnotationDbi::mapIds(org.Hs.eg.db, keys=uniq.symbols,
                                    column = "ENTREZID", keytype="SYMBOL",
                                    multiVals = "first")
symbol.vec <- symbol.vec[!is.na(symbol.vec)]
symbol2entrez <- data.frame(hgnc_symbol=names(symbol.vec),
                            entrezgene_id = as.character(symbol.vec),
                            stringsAsFactors = FALSE)

id <- match(rownames(de.res$mrna), symbol2entrez[,1])
de.res$mrna$EntrezID[!is.na(id)] <- symbol2entrez[id[!is.na(id)], 2]

uniq.ppi <- unique(c(V(g.string)$name, V(g.aracne)$name, V(g.ohmnet)$name, V(g.path)$name))

id <- match(rownames(de.res$methyl), symbol2entrez[,1])
de.res$methyl$EntrezID[!is.na(id)] <- symbol2entrez[id[!is.na(id)], 2]

methyl.entrez.dat <- de.res$methyl[!is.na(de.res$methyl$EntrezID),]
methyl.entrez.id <- methyl.entrez.dat$EntrezID
common <- intersect(methyl.entrez.id, uniq.ppi)

bp.el <- matrix(c(paste(rownames(methyl.entrez.dat)[match(common, methyl.entrez.dat$EntrezID)], "methyl", sep= "|"),
                   paste(common, "mrna", sep= "|")), ncol= 2)
bp.methyl.mrna <- igraph::graph_from_edgelist(el= bp.el, directed= FALSE)
V(bp.methyl.mrna)$node_type <- extract_string(V(bp.methyl.mrna)$name, "\\|", 2)
V(bp.methyl.mrna)$name <- extract_string(V(bp.methyl.mrna)$name, "\\|", 1)

vertex_attr(bp.methyl.mrna, "native_nom") <- ifelse(V(bp.methyl.mrna)$node_type == "mrna", "EntrezID", "Symbol")
id <- which(V(bp.methyl.mrna)$node_type == "mrna")
vertex_attr(bp.methyl.mrna, "EntrezID", id) <- V(bp.methyl.mrna)$name[id]
id <- which(V(bp.methyl.mrna)$node_type == "methyl")
vertex_attr(bp.methyl.mrna, "Symbol", id) <- V(bp.methyl.mrna)$name[id]
id <- which(V(bp.methyl.mrna)$node_type == "methyl" & V(bp.methyl.mrna)$name %in% rownames(methyl.entrez.dat))
vertex_attr(bp.methyl.mrna, 'EntrezID', id) <- methyl.entrez.dat$EntrezID[match(V(bp.methyl.mrna)$name[id], rownames(methyl.entrez.dat))]




#mirna-methyl bipartite
mirna.meth.nm <- rownames(de.res$methyl)[grep("MIR", rownames(de.res$methyl))]
mirna.de.nm <- rownames(de.res$mirna)
mirna.g.nm <- V(bp.mirna.mrna)$name[V(bp.mirna.mrna)$node_type == "mirna"]

new.names <- gsub(pattern= '-[35]p', replacement= '', x= mirna.g.nm)
new.names <- gsub(pattern= 'hsa-miR', replacement= 'MIR', x= new.names)
new.names <- gsub(pattern= 'hsa-let', replacement= 'MIRLET', x= new.names)
new.names <- toupper(new.names)
tips <- extract_string(new.names, "-", 1)
ends <- unlist(lapply(strsplit(new.names, "-"), function(x) paste(x[2:length(x)], collapse= "-")))
hyphen <- grepl('-', ends)
ABC.b4.hyphen <- grepl('[ABCDEFGHIJKLMNOPQRSTUVWXYZ]', extract_string(ends, '-', 1))
ends[hyphen & ABC.b4.hyphen] <- gsub(pattern= '-', replacement= '', x= ends[hyphen & ABC.b4.hyphen])
new.names <- paste0(tips, ends)

common <- intersect(new.names, rownames(de.res$methyl))

bp.el <- NULL
for(i in seq_along(common)){
  N <- sum(new.names == common[i])
  tmp.nm <- mirna.g.nm[new.names %in% common[i]]
  tmp <- matrix(c(paste(tmp.nm, 'mirna', sep= '|'),
                   paste(rep(common[i], N), 'methyl', sep= '|')), ncol= 2)
  bp.el <- rbind(bp.el, tmp)
}

bp.methyl.mirna <- igraph::graph_from_edgelist(el= bp.el, directed= FALSE)
V(bp.methyl.mirna)$node_type <- extract_string(V(bp.methyl.mirna)$name, "\\|", 2)
V(bp.methyl.mirna)$name <- extract_string(V(bp.methyl.mirna)$name, "\\|", 1)

vertex_attr(bp.methyl.mirna, "native_nom") <- ifelse(V(bp.methyl.mirna)$node_type == "mirna", "miRID", "Symbol")
id <- which(V(bp.methyl.mirna)$node_type == "mirna")
vertex_attr(bp.methyl.mirna, "miRID", id) <- V(bp.methyl.mirna)$name[id]
id <- which(V(bp.methyl.mirna)$node_type == "methyl")
vertex_attr(bp.methyl.mirna, "Symbol", id) <- V(bp.methyl.mirna)$name[id]
id <- which(V(bp.methyl.mirna)$node_type == "methyl" & V(bp.methyl.mirna)$name %in% rownames(methyl.entrez.dat))
vertex_attr(bp.methyl.mirna, 'EntrezID', id) <- methyl.entrez.dat$EntrezID[match(V(bp.methyl.mirna)$name[id], rownames(methyl.entrez.dat))]

#data list for AMEND
de.res$mrna$node.names <- de.res$mrna$EntrezID
de.res$mirna$node.names <- rownames(de.res$mirna)
de.res$methyl$node.names <- rownames(de.res$methyl)

param.id <- list(c(1,1), c(1,2), c(1,3), c(2,1), c(3,1))[[1]]
data.type <- c('logFC', 'P.Value', 'adj.P.Val')[param.id[1]]
pval.weight <- c(NULL, 'P.Value', 'adj.P.Val')[param.id[2]]

get_seeds <- function(data, data.type, pval.weight = NULL, node.id.col){
  id <- which(!is.na(data[, node.id.col]))
  if(data.type == 'logFC'){
    if(!is.null(pval.weight)){
      tmp <- data[id, data.type] * (1 - data[id, pval.weight])
    }else tmp <- data[id, data.type]
  }
  if(data.type %in% c('P.Value', 'adj.P.Val')){
    tmp <- data[, data.type]
  }
  names(tmp) <- data[id, node.id.col]
  tmp
}

data.list <- lapply(de.res, get_seeds, data.type= data.type, pval.weight= pval.weight, node.id.col= 'node.names')




#brw list for AMEND
de.res <- lapply(de.res, function(x){ x$hr <- rownames(x); x })

get_brw_attr <- function(brw.data, data.type, pval.weight = FALSE, de.data, node.id.col, brw.col){
  id <- which(!is.na(de.data[, node.id.col]))
  feat.nm <- de.data[id, brw.col]
  in.brw <- feat.nm %in% names(brw.data)
  id <- id[in.brw]
  feat.nm <- feat.nm[in.brw]
  tmp.dat <- brw.data[feat.nm]
  if(data.type == 'hr')
    tmp <- unlist(lapply(tmp.dat, function(x) x['hr']))
  if(data.type == '1/hr')
    tmp <- unlist(lapply(tmp.dat, function(x) 1 / x['hr']))
  if(data.type == '|log(hr)|')
    tmp <- unlist(lapply(tmp.dat, function(x) abs(log(x['hr']))))
  if(pval.weight){
    p <- unlist(lapply(tmp.dat, function(x) x['pval']))
    tmp <- tmp * (1 - p)
  }
  names(tmp) <- de.data[id, node.id.col]
  tmp
}

data.type <- c('hr', '1/hr', '|log(hr)|')[1]
pval.weight <- c(TRUE, FALSE)[2]

brw.list <- mapply(FUN= get_brw_attr, brw.data= hr, data.type= data.type,
                   pval.weight= pval.weight, de.data= de.res,
                   node.id.col= 'node.names', brw.col= 'hr')

input_graphs <- list('mirna;mrna'= bp.mirna.mrna,
                     'mirna;methyl'= bp.methyl.mirna,
                     'methyl;mrna'= bp.methyl.mrna,
                     'mrna_string'= g.string,
                     'mrna_path'= g.path,
                     'mrna_aracne'= g.aracne,
                     'mrna_ohmnet'= g.ohmnet)

saveRDS(list(input_graphs= input_graphs,data.list= data.list,brw.list= brw.list),paste0(path.results, "LUAD_net_checkpoint.rds"))

library(igraph)
library(fgsea)
library(DOSE)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(GO.db)
library(ggplot2)
library(dplyr)

path.results <- "Results/"
dir.create(paste0(path.results, "figures"), showWarnings= FALSE)

#load amend results
subnet.brw <- readRDS(paste0(path.results, "LUAD_n50_db.none_brw.rds"))
subnet.nobrw <- readRDS(paste0(path.results, "LUAD_n50_db.none.rds"))
subnet.main <- readRDS(paste0(path.results, "LUAD_n50_IN_agg_brw.rds"))

#ORA universe
chk <- readRDS(paste0(path.results, "LUAD_net_checkpoint.rds"))

universe.entrez <- unique(unlist(lapply(
  grep("^mrna", names(chk$input_graphs), value= TRUE),
  function(nm) {
    g <- chk$input_graphs[[nm]]
    if ("EntrezID" %in% vertex_attr_names(g)) {
      eids <- vertex_attr(g, "EntrezID")
      eids[!is.na(eids) & grepl("^[0-9]+$", eids)]
    } else {
      nms <- V(g)$name
      nms[grepl("^[0-9]+$", nms)]
    }
  }
)))

universe.disease <- unique(unlist(lapply(names(chk$input_graphs), function(nm) {
  g <- chk$input_graphs[[nm]]
  if ("EntrezID" %in% vertex_attr_names(g)) {
    eids <- vertex_attr(g, "EntrezID")
    eids[!is.na(eids) & grepl("^[0-9]+$", eids)]
  } else 
      character(0)
})))

rm(chk)
cat("Universe GO (mRNA):", length(universe.entrez), "\n")
cat("Universe disease (all):", length(universe.disease), "\n")

#extract nodes from module
extract_module_nodes <- function(subnet) {
  if (is_igraph(subnet)) 
    return(V(subnet)$name)
  for (nm in c("module", "subnetwork", "subnet", "graph")) {
    if (nm %in% names(subnet) && is_igraph(subnet[[nm]])) 
      return(V(subnet[[nm]])$name)
  }
  for (nm in names(subnet)) {
    if (is_igraph(subnet[[nm]])) 
      return(V(subnet[[nm]])$name)
  }
  stop("module format not recognized")
}

nodes.brw <- extract_module_nodes(subnet.brw)
nodes.nobrw <- extract_module_nodes(subnet.nobrw)
nodes.main <- extract_module_nodes(subnet.main)

cat("\nRun 1 (BRW):", length(nodes.brw), "noeuds\n")
cat("Run 2 (no BRW):", length(nodes.nobrw), "noeuds\n")
cat("Run 3 (BRW+IN+agg):", length(nodes.main), "noeuds\n")

#extract entrez ids from nodes (only mRNA)
extract_entrez <- function(nodes) {
  mrna_nodes <- nodes[grepl("\\|mrna", nodes)]
  ids <- sub("\\|.*", "", mrna_nodes)
  unique(ids[grepl("^[0-9]+$", ids)])
}

genes.main.entrez <- extract_entrez(nodes.main)
genes.brw.entrez <- extract_entrez(nodes.brw)
genes.nobrw.entrez <- extract_entrez(nodes.nobrw)

cat("\nRun 3 : Entrez IDs mRNA:", length(genes.main.entrez), "\n")
for (typ in c("mrna", "mirna", "methyl")) {
  n <- sum(grepl(paste0("\\|", typ), nodes.main))
  if (n > 0) 
    cat(" ", typ, ":", n, "\n")
}

write.csv(data.frame(gene= genes.main.entrez),paste0(path.results, "module_genes_run3.csv"), row.names=FALSE)
write.csv(data.frame(gene= genes.brw.entrez),paste0(path.results, "module_genes_run1.csv"), row.names=FALSE)
write.csv(data.frame(gene= genes.nobrw.entrez), paste0(path.results, "module_genes_run2.csv"), row.names=FALSE)

#construction gene sets GO
cat("\nConstruction des gene sets GO...\n")
go_map <- suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db,
  keys= keys(org.Hs.eg.db, "ENTREZID"),
  columns= c("GO", "ONTOLOGY"),
  keytype= "ENTREZID"
))
go_map <- go_map[!is.na(go_map$GO) & !is.na(go_map$ONTOLOGY), ]

make_go_sets <- function(ont_code) {
  sub <- go_map[go_map$ONTOLOGY == ont_code, ]
  lapply(split(sub$ENTREZID, sub$GO), unique)
}

go <- list(
  go_BP= make_go_sets("BP"),
  go_MF= make_go_sets("MF"),
  go_CC= make_go_sets("CC")
)
cat("GO sets: BP=", length(go$go_BP), "MF=", length(go$go_MF), "CC=", length(go$go_CC), "\n")

add_go_desc <- function(fora_res) {
  if (nrow(fora_res) == 0)
    return(fora_res)
  terms <- suppressMessages(AnnotationDbi::select(GO.db, keys = fora_res$pathway, columns = "TERM", keytype = "GOID"))
  fora_res$Description <- terms$TERM[match(fora_res$pathway, terms$GOID)]
  fora_res
}

#enrichment GO : BP+MF+CC (fora)
p.adj.cutoff <- 0.01
cat("\n=== GO Enrichment (fora, padj <=", p.adj.cutoff, ") ===\n")

ora.main <- NULL
for (aspect in names(go)) {
  res.tmp <- tryCatch(
    fora(pathways= go[[aspect]], genes= genes.main.entrez, universe= universe.entrez) %>%
      filter(padj<= p.adj.cutoff) %>%
      mutate(
        overlapGenes= sapply(overlapGenes, paste, collapse = ", "),
        GO_aspect=sub("go_", "", aspect)
      ),
    error= function(e) { message("fora error (", aspect, "): ", conditionMessage(e)); NULL }
  )
  if (!is.null(res.tmp) && nrow(res.tmp) > 0) {
    res.tmp <- add_go_desc(res.tmp)
    ora.main <- rbind(ora.main, as.data.frame(res.tmp))
  }
}

if (is.null(ora.main) || nrow(ora.main) == 0) {
  cat("Aucun terme GO significatif\n")
} else {
  ora.main <- ora.main[order(ora.main$padj), ]
  cat("Termes significatifs:", nrow(ora.main),
      "(BP:", sum(ora.main$GO_aspect == "BP"),
      "MF:", sum(ora.main$GO_aspect == "MF"),
      "CC:", sum(ora.main$GO_aspect == "CC"), ")\n")
  print(ora.main[, c("Description", "GO_aspect", "overlap", "size", "padj")])
  write.csv(ora.main, paste0(path.results, "GO_ora_run3.csv"), row.names=FALSE)

  plot.dat <- ora.main
  plot.dat$GeneRatio <- plot.dat$overlap / plot.dat$size
  plot.dat$Description <- factor(plot.dat$Description,
    levels= rev(plot.dat$Description[order(plot.dat$padj)]))
  plot.dat$GO_aspect <- factor(plot.dat$GO_aspect, levels = c("BP", "MF", "CC"))

  p.go <- ggplot(plot.dat, aes(x= GeneRatio, y=Description,color= padj, size= overlap)) +
    geom_point()+
    scale_color_gradient(low = "red2", high = "steelblue", name = "padj", limits = c(0, p.adj.cutoff)) +
    scale_size_continuous(name = "Overlap") +
    scale_y_discrete(expand = expansion(add = 0.8)) +
    facet_grid(GO_aspect ~ ., scales = "free_y", space = "free_y") +
    labs(title = "GO Enrichment : Active Module LUAD (Run 3)", x = "Gene Ratio", y = NULL) +
    theme_bw(base_size = 12) +
    theme(strip.text = element_text(face = "bold"),axis.text.y= element_text(size = 10), panel.spacing = unit(0.8, "lines"))

  ggsave(paste0(path.results, "figures/GO_ora_run3.pdf"), p.go, width=11, height=7)
  ggsave(paste0(path.results, "figures/GO_ora_run3.png"), p.go, width=11, height=7, dpi=150)
}

#enrrichment disease : HDO+ DGN 
cat("\n=== Disease Ontology (HDO, padj <=", p.adj.cutoff, ") ===\n")
edo.main <- tryCatch(
  enrichDO(gene= genes.main.entrez, ont= "HDO", pAdjustMethod= "BH", universe= universe.disease),
  error= function(e) { message("enrichDO error: ", conditionMessage(e));NULL }
)
if (!is.null(edo.main)) {
  edo.main@result <- edo.main@result[!is.na(edo.main@result$p.adjust) &
    edo.main@result$p.adjust <= p.adj.cutoff, ]
  cat("Termes DO:", nrow(edo.main@result), "\n")
  if (nrow(edo.main@result) > 0) {
    print(edo.main@result[, c("Description", "GeneRatio", "p.adjust")])
    saveRDS(edo.main, paste0(path.results, "DO_run3.rds"))
    write.csv(edo.main@result, paste0(path.results, "DO_run3.csv"), row.names=FALSE)
  }
}

cat("\n=== DisGeNET (padj <=", p.adj.cutoff, ") ===\n")
edgn.main <- tryCatch(
  enrichDGN(gene= genes.main.entrez, pAdjustMethod= "BH", universe= universe.disease),
  error= function(e) { message("enrichDGN error: ", conditionMessage(e)); NULL }
)
if (!is.null(edgn.main)) {
  edgn.main@result <- edgn.main@result[!is.na(edgn.main@result$p.adjust) &
    edgn.main@result$p.adjust <= p.adj.cutoff, ]
  cat("Termes DGN:", nrow(edgn.main@result), "\n")
  if (nrow(edgn.main@result) > 0) {
    print(edgn.main@result[, c("Description", "GeneRatio", "p.adjust")])
    saveRDS(edgn.main, paste0(path.results, "DGN_run3.rds"))
    write.csv(edgn.main@result, paste0(path.results, "DGN_run3.csv"), row.names=FALSE)
  }
}

#module comparison (Jaccard)
jaccard <- function(a, b) round(length(intersect(a, b)) / length(union(a, b)), 3)
cat("\nJaccard Run1 vs Run2:", jaccard(genes.brw.entrez, genes.nobrw.entrez), "\n")
cat("Jaccard Run1 vs Run3:", jaccard(genes.brw.entrez, genes.main.entrez), "\n")
cat("Jaccard Run2 vs Run3:", jaccard(genes.nobrw.entrez, genes.main.entrez), "\n")


#main result (run 3): module visualization
g.plot <- if (is_igraph(subnet.main)) {
  subnet.main
} else {
  subnet.main[[which(sapply(subnet.main, is_igraph))[1]]]
}

if (is_igraph(g.plot)) {
  node.ids <- sub("\\|.*", "", V(g.plot)$name)
  node.type <- sub(".*\\|", "", V(g.plot)$name)

  mrna.idx <- which(grepl("mrna", node.type))
  symbols <- suppressMessages(mapIds(org.Hs.eg.db,
    keys= node.ids[mrna.idx], column= "SYMBOL",
    keytype= "ENTREZID", multiVals= "first"))
  labels <- node.ids
  labels[mrna.idx] <- ifelse(!is.na(symbols[node.ids[mrna.idx]]),symbols[node.ids[mrna.idx]], node.ids[mrna.idx])
  node.colors <- ifelse(grepl("mirna", node.type), "steelblue", "tomato")

  set.seed(42)
  lay <- layout_with_fr(g.plot, niter=5000)
  lay <- norm_coords(lay, xmin=-8, xmax=8, ymin=-8, ymax=8)

  write.csv(data.frame(id= V(g.plot)$name, symbol=labels, type= node.type),
    paste0(path.results, "network_nodes_run3.csv"), row.names=FALSE)
  write.csv(igraph::as_data_frame(g.plot, what= "edges"),
    paste0(path.results, "network_edges_run3.csv"), row.names=FALSE)

  do_plot <- function(){
    plot(g.plot,
      layout= lay, rescale= FALSE, xlim= c(-10, 10), ylim= c(-10, 10),
      vertex.color= node.colors, vertex.size= 17,
      vertex.label= labels, vertex.label.cex= 0.45,
      vertex.label.color= "black", vertex.label.font = 2,
      vertex.label.dist = 1.4, vertex.label.degree= -pi / 2,
      edge.width= 0.8, edge.color = "grey70",
      main= "Active module LUAD - Run 3")
    legend("topright", legend = c("mRNA", "miRNA"),
      fill = c("tomato", "steelblue"), title= "Node type", bty = "n")
  }

  tryCatch({
    pdf(paste0(path.results,"figures/network_module_run3.pdf"), width=18, height=18)
    do_plot(); dev.off()
    png(paste0(path.results,"figures/network_module_run3.png"), width=1800, height=1800,res=150)
    do_plot(); dev.off()
    cat("Module saved (PDF + PNG)\n")
  }, error = function(e) {
    try(dev.off(), silent = TRUE)
    cat("Erreur reseau:", conditionMessage(e), "\n")
  })
}

cat("\nFichiers Results/:\n")
cat(list.files(path.results, pattern= "\\.csv$|\\.rds$"), sep= "\n")

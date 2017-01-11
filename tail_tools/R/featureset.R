
#
# A featureset is a list(set=, universe= )
#

get_organism_db <- function(species) {
    if (is.null(species)) 
        return(NULL)

    if (species == "Sc"){
        return(org.Sc.sgd.db::org.Sc.sgd.db)
    } else if (species == "Hs"){
        return(org.Hs.eg.db::org.Hs.eg.db)
    } else if(species == "Ce"){
        return(org.Ce.eg.db::org.Ce.eg.db)
    } else if(species == "Mm"){
        return(org.Mm.eg.db::org.Mm.eg.db)
    }
    
    stop(paste("Unsupported species:", species))
}


.get_genesets <- function(species) {
    org_db <- get_organism_db(species)
    if (is.null(org_db)) return(NULL)
    
    gene_set <- 
        AnnotationDbi::select(
            org_db, 
            AnnotationDbi::keys(org_db, "ENSEMBL"), 
            keytype="ENSEMBL", 
            columns="GOALL") %>%
        select_(name=~ENSEMBL, set_name=~GOALL) %>%
        unique()
    
    set_info <-
        AnnotationDbi::select(
            GO.db::GO.db,
            unique(gene_set$set_name),
            c("TERM","ONTOLOGY")) %>%
        select_(set_name=~GOID, description=~TERM, ontology=~ONTOLOGY)
    
    # Taken from limma "kegga" function
    species.KEGG <- switch(species, Ag = "aga", At = "ath", 
            Bt = "bta", Ce = "cel", Cf = "cfa", Dm = "dme", Dr = "dre", 
            EcK12 = "eco", EcSakai = "ecs", Gg = "gga", Hs = "hsa", 
            Mm = "mmu", Mmu = "mcc", Pf = "pfa", Pt = "ptr", 
            Rn = "rno", Ss = "ssc", Xl = "xla",
            
            Sc = "sce")
    
    if (!is.null(species.KEGG)) {
        keggs <- limma::getGeneKEGGLinks(species.KEGG, convert=TRUE)
        conversion <- AnnotationDbi::select(org_db, unique(keggs$GeneID), 
            keytype="ENTREZID",  columns="ENSEMBL")
        
        keggs <- keggs %>%
            select_(ENTREZID=~GeneID, set_name=~PathwayID) %>%
            inner_join(conversion, "ENTREZID") %>%
            select_(name=~ENSEMBL, ~set_name)
        
        kegg_info <- limma::getKEGGPathwayNames(species.KEGG, remove.qualifier=TRUE) %>%
            select_(set_name=~PathwayID, description=~Description) %>%
            mutate_(ontology=~"KEGG")
        
        gene_set <- rbind(gene_set, keggs)
        set_info <- rbind(set_info, kegg_info)
    }
    
    list(
        gene_set = gene_set, 
        set_info = set_info)
}


get_genesets <- function(species) {
    cached(paste0("cached_genesets_",species), .get_genesets, species, version=2)
}


featureset_test_genesets <- function(fs, species, minimum_set_size=2, ontologies=c("BP","MF","CC","KEGG")) {
    sets <- get_genesets(species)

    set_info <- sets$set_info %>%
         filter_(~ontology %in% ontologies)
    gene_set <- sets$gene_set %>%
         filter_(~set_name %in% set_info$set_name)
    
    universe <- intersect(as.character(fs$universe), gene_set$name)
    set <- intersect(as.character(fs$set), universe)

    gene_set <- gene_set %>%
        filter_(~ name %in% universe) %>%
        group_by_(~ set_name) %>%
        filter_(~ length(name) >= minimum_set_size) %>%
        ungroup() %>%
        as.data.frame
        
    #cat(length(fs$universe), length(universe), length(fs$set), length(set), "\n")

    # kegga is broken if universe supplied for limma 3.4
    result <- limma::kegga(set, #universe, 
        gene.pathway=gene_set, pathway.name=sets$set_info)
    result$set_name <- rownames(result)
    result$FDR <- p.adjust(result$P.DE, method="BH")
    result <- 
        left_join(result, sets$set_info, "set_name") %>%
        select_(~FDR, p_value=~P.DE, size=~N, overlap=~DE, ~set_name, ~ontology, ~description) %>%
        arrange_(~p_value)
    
    list(table=result, set=set, universe=universe)    
}




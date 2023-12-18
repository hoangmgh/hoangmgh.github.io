---
layout: post
title: 📢 Streamline your WGCNA analysis pipeline 
date: 2022-1-12
description: Writing reproducible pipeline for your WGCNA analysis 
tags: math WGCNA systems-biology
categories: Bioinformatics
related_posts: false
thumbnail: /assets/img/network.png
---
As I work with WGCNA to construct  multiple single-cell genes network, I realize that the whole process is quite cumbersome, and the analysis will benefit greatly if we can make reproducible functions for parallel processing of multiple single-cell populations, power detection, module construction, and downstream analysis. Driven by this, I wrote a few functions for my own use. Feel free to use them if you are interested.

The first step in WGCNA network construction is preprocessing of raw reads. Often times for RNA-seq  data, such as scRNA-seq, one simple way  to normalize would be CPM normalization, followed by log1p or log10p. We don't have to provide scaled data into WGCNA's pipeline - it is not necessary, and scaling doesnt change correlation, either. Moreover, there are some internal steps in WGCNA where they would check for missing data or sparseness (in the case that one have many 0s across genes).
```R
preprocess_input=function (meta_file, count_file, cluster_label = "cluster", sample_label = "sample", 
    vars_to_keep = c("sample", "immune_infiltration_level", "Condition"), 
    normalize = FALSE, cluster_to_use = c("all"), normalize_method = "edgeR", 
    remove_genes = Ig_rm, vst = TRUE,remove_low_count=NA
                          ) 
    '''
    parameters:
        meta_file: a metadata file for all population of cells, with
                    cluster_label: the columns correspond to cell-type (cluster)
                    sample_label: the column of sampleID
                    vars_to_keep: are the variables to be retained for further analyses
        normalize:whether to normalize the count data
        normalize_method: using existing package for normalization, such as edgeR or deSeq2
        vst: whether to use variance stabilizing transformation alongside with deseq2 normalization
        cluster_to_use: indicating the cell population to be used, default to all existing populations 
        remove_low_count: remove genes with a total sum of counts surpassing this value, or NA if use all available genes
        removes_genes: remove genes from a given blacklist, or NA to keep the existing set of geness
    output:
        raw_data: keep the raw counts for all populations/clusters
        normalized_data: keep the normalized counts for the given clusters in cluster_to_use
        metadata: the metadata for each population, seperately. The sample names will align with the rownames of the count matrix
 
    '''
{
    print("1. read in metadata:")
    ################################################
    dat_meta <<- data.frame(as.matrix(fread(meta_file), rownames = 1))
    clusters = unique(dat_meta[[cluster_label]])
    if (length(cluster_to_use) == 1 & cluster_to_use[1] == "all") {
        print("use all clusters")
        cluster_to_use = clusters
    }
    else {
        if (any(cluster_to_use %in% dat_meta[[cluster_label]])) {
            print("one of cluster is not in datameta")
        }
        
    }
    print(paste0("using a total of:", length(cluster_to_use), 
        " unique clusters"))
    ################################################
    dat_meta_splitted = lapply(clusters, function(i) {
        dat_meta[which(dat_meta[[cluster_label]] == i), ]
    })
    names(dat_meta_splitted) = clusters
    print("2. sanity check  metadata:")
    if (sum(!table(dat_meta[[cluster_label]], dat_meta[[sample_label]]) <= 
        1)) {
        print("more than 2 identical samples  exist in one cluster")
        break
    }
    if ( !sample_label %in% colnames(dat_meta_splitted)) {
        print("sample label must exist in metadata")
    }
    if ( !all(vars_to_keep %in% colnames(dat_meta_splitted))) {
        print("all variables used  must exist in metadata")
    }
    if (!sample_label %in% vars_to_keep ) {
        print("sample label and cluster label must be in vars_to_keep")
        break
    }
    
    dat =as.matrix(fread(count_file), rownames = 1)
   
    if (length(remove_genes)>1 & !is.na(remove_genes[1])){
        print("3. remove gene in the blacklist:")
        dat = dat[which(!rownames(dat) %in% remove_genes), ]       
    }
      print(paste0("remaining genes", nrow(dat)))
    dat = lapply(clusters, function(i) t(dat[, which(dat_meta[[cluster_label]] == 
        i)]))
                 
    names(dat) = clusters
    if (!is.na(remove_low_count)){
        print("removing genes with low counts overall:")

        keep = lapply(clusters, function(i) colnames(dat[[i]])[which(colSums(dat[[i]],na.rm=T) >=remove_low_count)])
        print(paste0("remaining genes: ", length(keep[[1]])))
    }
    else{
    print(clusters)
    keep = lapply(clusters, function(i) colnames(dat[[i]]))  
    }                      
    names(keep) = clusters
    print("replace NA with 0 for calculating CPM:")
    na_mask = lapply(clusters, function(i) apply(dat[[i]],2,is.na))
    names(na_mask)=clusters
    for (i in clusters){
           dat[[i]][na_mask[[i]]]=0
    }
 
    for (i in seq_along(dat_meta_splitted)) rownames(dat_meta_splitted[[i]]) <- dat_meta_splitted[[i]][[sample_label]]
    for (i in seq_along(dat_meta_splitted)) rownames(dat[[i]]) <- rownames(dat_meta_splitted[[i]])
    if (normalize) {
        if (normalize_method == "deseq2") {
            print("using deseq2 normalization")
            dds = lapply(cluster_to_use, function(i) DESeqDataSetFromMatrix(countData = t(dat[[i]]), 
                colData = data.frame(dat_meta_splitted[[i]]), 
                design = ~1))
            if (vst == TRUE) {
                print("using VST")
                vsd <- lapply(dds, function(i) vst(i, blind = FALSE))
                mat <- lapply(vsd, function(i) t(assay(i)))
                names(mat) = cluster_to_use
                multiExpr = lapply(cluster_to_use, function(i) list(data = mat[[i]][, 
                  keep[[i]]]))
                                   
                names(multiExpr) = cluster_to_use
            }
            else {
                names(dds) = cluster_to_use
                dds = lapply(dds, estimateSizeFactors)
                multiExpr = lapply(cluster_to_use, function(i) list(data = t(log10(counts(dds[[i]], 
                  normalized = T) + 1))))
                names(multiExpr) = cluster_to_use
            }
                                   
        }
                                   
        else if (normalize_method == "edgeR") {
            print("using edgeR normalization:")
            multiExpr = lapply(cluster_to_use, function(i) {
            list(data = t(cpm(t(dat[[i]]), log = T)))})
            names(multiExpr) = cluster_to_use
            for (i in clusters){
               multiExpr[[i]]$data[na_mask[[i]]]=NA
               dat[[i]][na_mask[[i]]]=NA
            }

            multiExpr = lapply(cluster_to_use, function(i) {
                selected=colnames(multiExpr[[i]]$data) %in% keep[[i]]
                return(list(data=multiExpr[[i]]$data[,selected]))
            }
            )
            names(multiExpr) = cluster_to_use

        }        
    }
    else {
        multiExpr = lapply(cluster_to_use, function(i) list(data = dat[[i]][, 
         keep[[i]]]))
        names(multiExpr) = cluster_to_use
    }
    return(list(raw_data=dat, normalized_data=multiExpr, metadata=dat_meta_splitted))
}
```

the next step would be finding the appropriate powers for each population of cell. With this we have the find_power() function:

```R
find_power=function(multiExpr,which_cluster="all",meta,blockSize=30000,verbose=2,corFnc=WGCNA::cor) {
    '''
    multiExpr: a list of normalized datamatrix created from calling preprocess_input()
    which_cluster: which cluster to find power for. default to all
    blockSize: number of genes to be processed in parallel
    corFnc: correlation function to be used, default to WGCNA s cor 
    '''

    if (length(which_cluster)==1 &which_cluster=="all"){
        nSet=length(multiExpr)
        which_cluster=names(multiExpr)
    }
    else {
    nSet=length(which_cluster)
    }
    
    powers <<- c(seq(2,10,by=1), seq(12,20, by=2));
    # Initialize a list to hold the results of scale-free analysis
    powerTables <<- vector(mode = "list", length = nSet);
    # Call the network topology analysis function for each set in turn
    for (set in seq_along(which_cluster)){
      powerTables[[set]] <<-list(data = pickSoftThreshold(multiExpr[[which_cluster[set]]]$data, powerVector=powers,
                                                         verbose = 2,corFnc=corFnc,blockSize=30000)[[2]])
      
    collectGarbage();
    # Plot the results:
    # Will plot these columns of the returned scale free analysis tables
    names(powerTables)=which_cluster
    return(powerTables)
}
```

After calling this function, one would go on to construct the network, the most time-consuming step:


```R
construct_network=function(multiExpr,which_cluster="all",nThreads=30,power_list,minModuleSize = 50, reassignThreshold = 0, 
                           mergeCutHeight = .15,numericLabels = TRUE, pamRespectsDendro = TRUE, saveTOMs = FALSE,
                      verbose = 3,maxBlockSize=30000,deepSplit=4,  checkMissingData = TRUE,
                      corFnc=WGCNA::cor,save_output=NA) {
'''
    multiExpr:the list of normalized matrices from step1 (preprocess_input)
    which_cluster: the list of cluster to use, default to all ("all")
    nThreads: number of threads to run WGCNA
    power_list: the list of power detected from step2 (find_powers()), must match the number of clusters to process
    for other parameters, consult WGCNA manual 
'''
if (length(which_cluster)==1 &which_cluster=="all") cluster_to_use= names(multiExpr)
else cluster_to_use=which_cluster

if (length(cluster_to_use) != length(power_list)){
    print("the length of the power list must match the number of clusters")
}    
names(power_list)=cluster_to_use

net=lapply(cluster_to_use,function(i) {
     blockwiseModules(multiExpr[[i]]$data,nThreads=nThreads, power =power_list[[i]],
                      minModuleSize = minModuleSize, reassignThreshold = reassignThreshold, mergeCutHeight = mergeCutHeight, 
                      numericLabels = numericLabels, pamRespectsDendro = pamRespectsDendro, saveTOMs = saveTOMs,
                      verbose = verbose,maxBlockSize=maxBlockSize,deepSplit=deepSplit,  checkMissingData = checkMissingData,
                      corFnc=corFnc)})
names(net)=names(cluster_to_use)
if (!is.na(save_output)){
saveRDS(net,save_output)
}
return(net)
}
```
The next step would be processing the assignments of genes, calculating module eigengenes, gene module membership, and the correlating module's trend with other variables of interest. The most important parameter here is summary, whether we should take the average of the gene to calculate module eigengenes or use the native WGCNA approach (1st PC):

```R
analyze_network=function (net, multiExpr, meta, use = "pairwise.complete.obs", 
    scale = T, summary = "pc") 
{
    '''
    parameter:
        net: the list of network constructed from the previous step (construct_network)
        multiExpr:the list of normalized matrices from step1 (preprocess_input)
        meta: the metadata obtained from step1
        use: this will be feed into the correlation function. for example use="p" means use the none NA entries
        scale: whether to scale the data before correlation calculation. this is advised for using the mean trend strategy (summary ="mean")
        summary: whether to take the average value of the genes within each module ("mean"), or use WGCNA  module eigengene ("pc") 


    output: 
        cor_df: a dataframe storing  correlation coefficients between each module and variables of interest 
        p_df: a dataframe storing  correlation p-value between each module and variables of interest 
        group_gene: a named lists where each element is a list of genes per module
        geneModuleMembership: the membership of each gene to all modules
        mergedColors: the assignment of each gene to their module (color as the module nmae)

    '''
    mergedColors = lapply(net, function(i) labels2colors(i$colors))
    names(mergedColors) =  names(net)
    geneModuleMembership = lapply(names(net), function(i) as.data.frame(WGCNA::cor(multiExpr[[i]]$data, 
        net[[i]]$MEs, use = use,method="spearman")))
    names(geneModuleMembership) = names(net)
    group_gene = lapply(names(net), function(i) split(colnames(multiExpr[[i]]$data), 
        mergedColors[[i]]))
    MMPvalue = lapply(names(net), function(i) as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership[[i]]), 
        ncol(net[[i]]))))
    print("test")
    names(group_gene) = names(net)
    names(MMPvalue) = names(net)
    names(mergedColors) = names(net)
    groups = names(net)
    cor_result = lapply(names(net), function(i) {
        nGenes = ncol(multiExpr[[i]]$data)
        nSamples = nrow(multiExpr[[i]]$data)
        if (scale) {
            tmp_data = scale(multiExpr[[i]]$data)
        }
        else {
            tmp_data = (multiExpr[[i]]$data[rownames(meta[[i]]), 
                ])
        }
        print(dim(tmp_data))
        if (summary == "mean") 
            MEs = orderMEs(moduleMeangenes(tmp_data, mergedColors[[i]]))
        else if (summary =="pc") MEs = orderMEs(moduleEigengenes(tmp_data, mergedColors[[i]])$eigengenes)
        moduleTraitCor = WGCNA::cor(MEs, meta[[i]], use = use, 
            method = "pearson")
        moduleTraitPvalue = corPvalueStudent(moduleTraitCor, 
            nSamples)
        rownames(moduleTraitCor) = paste0(i, "_", rownames(moduleTraitCor))
        rownames(moduleTraitPvalue) = paste0(i, "_", rownames(moduleTraitPvalue))
        return(list(Cor = moduleTraitCor, p = moduleTraitPvalue))
    })
    cor_df = data.frame(do.call(plyr::rbind.fill, lapply(cor_result, 
        function(i) data.frame(i$Cor, rownames = rownames(i$Cor)))))
    rownames(cor_df) = unlist(sapply(cor_result, function(i) rownames(i$Cor)))
    p_df = data.frame(do.call(plyr::rbind.fill, lapply(cor_result, 
        function(i) data.frame(i$p, rownames = rownames(i$p)))))
    rownames(p_df) = unlist(sapply(cor_result, function(i) rownames(i$p)))
    return(list(cor_df = cor_df, p_df = p_df,group_gene=group_gene,geneModuleMembership=geneModuleMembership,
               mergedColors=mergedColors))
}
```

With these helper functions, the proper way to construct networks across multiple cell populations is like the following:
```R
step1=preprocess_input(meta_file="/projects/Thyroid/Results_nuc/7_NetworkAnalysis/Eco_system/Eco_all_meta.csv",
                          count_file="/projects/Thyroid/Results_nuc/7_NetworkAnalysis/Eco_system/Eco_all_count.csv",
                       vst=T,cluster_label="cluster",sample_label="Channel"
                ,vars_to_keep=c("Channel","immune_infiltration_percent","Condition"),normalize=T,
                          cluster_to_use=c("all"),normalize_method="edgeR",
                        remove_genes=Ig_rm,remove_missing=10)
step1[[3]]$Eco_all$log_immune_infiltration_percent=log10(as.numeric(step1[[3]]$Eco_all$immune_infiltration_percent))

find_power(step1[[2]])

net=construct_network(step1[[2]],save_output="/working-dir/net.RDS")

cor_result=analyze_network(net,step1[[2]],step1[[3]],summary="mean")

```
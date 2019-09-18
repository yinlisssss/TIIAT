OvarianCancerExpr <- system.file("extdata", "sample_input.txt", package="estimate")

OvarianCancerExpr 

filterCommonGenes1 <- function (input.f, output.f, id = c("GeneSymbol", "EntrezID")) 
{
  id <- match.arg(id)
  input.df <- input.f
  merged.df <- merge(common_genes, input.df, by.x = id, by.y = "row.names")
  rownames(merged.df) <- merged.df$GeneSymbol
  merged.df <- merged.df[, -1:-ncol(common_genes)]
  print(sprintf("Merged dataset includes %d genes (%d mismatched).", 
                nrow(merged.df), nrow(common_genes) - nrow(merged.df)))
  outputGCT(merged.df, output.f)
}


www <- list(ACC_unique,
            BLCA_unique,
            BRCA_unique,
            CESC_unique,
            CHOL_unique,
            COAD_unique,
            DLBC_unique,
            ESCA_unique,
            GBM_unique,
            HNSC_unique,
            KICH_unique,
            KIRC_unique,
            KIRP_unique,
            LAML_unique,
            LGG_unique,
            LIHC_unique,
            LUAD_unique,
            LUSC_unique,
            MESO_unique,
            OV_unique,
            PAAD_unique,
            PCPG_unique,
            PRAD_unique,
            READ_unique,
            SARC_unique,
            SKCM_unique,
            STAD_unique,
            TGCT_unique,
            THCA_unique,
            THYM_unique,
            UCEC_unique,
            UCS_unique,
            UVM_unique)

saveRDS(www,'alltumor.RData')

www <- readRDS('alltumor.RData')

purity_tumor <- list()

for (i in 1:33) {
  
  print(paste('trying for',i))
  
  filterCommonGenes1(input.f=www[[i]], output.f=paste0("unique_tumor_purity",i,".gct"), id="GeneSymbol")
  
  estimateScore(paste0("unique_tumor_purity",i,".gct"),output.ds = paste0("unique_tumor_purity",i,".csv"), 
                platform = 'illumina')
  
  aa <- read.csv(paste0("unique_tumor_purity",i,".csv"),sep = '\t',skip = 1)
  
  aa <- as.data.frame(t(aa))
  
  colnames(aa) <- aa[1,]
  
  aa <- aa[c(-1,-2),]
  
  aa$ESTIMATEScore <- as.numeric(aa$ESTIMATEScore)
  
  aa$purity <- cos(0.6049872018 + 0.0001467884*aa$ESTIMATEScore)
  
  aa$NAME <- str_replace_all(aa$NAME,"[.]","-")
  
  row.names(aa) <- aa$NAME
  
  purity_tumor[[i]] <- aa
}













###calculate the top surival gene

###############RFS

OS_cancer <- unique(survival_tcga_rfs_origin$cancer.type.abbreviation)
OS_cancer
OS_cancer <- OS_cancer[c(1:8,10:29)]
OS_cancer
length(OS_cancer)

www <- list(ACC_unique,
            BLCA_unique,
            BRCA_unique,
            CESC_unique,
            CHOL_unique,
            COAD_unique,
            DLBC_unique,
            ESCA_unique,
            HNSC_unique,
            KICH_unique,
            KIRC_unique,
            KIRP_unique,
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
            STAD_unique,
            TGCT_unique,
            THCA_unique,
            UCEC_unique,
            UCS_unique)

final_surivival_top1000 <- list()


for (i in 1:28) {
  
OS_cancer_1 <- filter(survival_tcga_rfs_origin,survival_tcga_rfs_origin$cancer.type.abbreviation==OS_cancer[i])

aaa <- intersect(OS_cancer_1$sample,colnames(www[[i]]))

exp_survi <- www[[i]][,aaa]

exp_survi <- as.data.frame(t(exp_survi))

exp_survi$sample <- row.names(exp_survi)

exp_sur <- merge(survival_tcga_rfs_origin,exp_survi,by='sample')

exp_sur[1:10,1:10]

exp_sur <- exp_sur[,c(-2,-3)]

row.names(exp_sur) <- exp_sur$sample

exp_sur[1:10,1:10]

exp_sur <- exp_sur[,-1]

genes <- colnames(exp_sur)[-c(1:2)]

print('calculate p-value')

system.time(res3 <- mclapply(1:length(genes),mc.cores = 8,function(j){
  group = ifelse(exp_sur[,genes[j]]>median(exp_sur[,genes[j]]),"high","low")
  if(length(table(group))==1) return(NULL)
  surv =as.formula(paste('Surv(DFI.time, DFI)~', "group"))
  data = cbind(exp_sur[,1:2],group)
  data$group <- as.factor(data$group)
  data$group <- relevel(data$group,ref = 'low')
  x = survdiff(surv,data=data)
  p.val = 1 - pchisq(x$chisq, length(x$n) - 1)
  HR = (x$obs[2]/x$exp[2])/(x$obs[1]/x$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/x$exp[2]+1/x$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/x$exp[2]+1/x$exp[1]))
  return(c(genes[j],p.val, HR, up95, low95))
}))

aaaaaaa <- data.frame(do.call(rbind,res3))

colnames(aaaaaaa) <- c('Gene','p.value','HR','95%CI upper','95%CI lower')

for (k in 2:5) {
  aaaaaaa[,k] <- as.numeric(as.vector(aaaaaaa[,k]))
}

wwwwwwww <- filter(aaaaaaa,aaaaaaa$p.value<0.05)

wwwwwwww[wwwwwwww==Inf]<-NA

wwwwwwww <- na.omit(wwwwwwww)

print(paste('done',i))
      
final_surivival_top1000[[i]] <- wwwwwwww

}


final_surivival_signi <- list()

for (i in 1:28) {
  
  nnn_great <- filter(final_surivival_top1000[[i]],final_surivival_top1000[[i]]$`95%CI lower`>1,
                      final_surivival_top1000[[i]]$`95%CI upper`>1)
  
  nnn_lower <- filter(final_surivival_top1000[[i]],final_surivival_top1000[[i]]$`95%CI lower`<1,
                      final_surivival_top1000[[i]]$`95%CI upper`<1)
  
  final_surivival_signi[[i]] <- rbind(nnn_great,nnn_lower)
  
}



saveRDS(final_surivival_signi,'topsurvival_rfs.RData',compress = F)


##




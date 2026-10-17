library(arrow)
set.seed(20261016)
root <- "draft/nature communications/Revision/Revised_code"
input <- file.path(root, "data/Simulated_data_dml.parquet")
base <- read_parquet(input)
if (nrow(base) > 1000) base <- base[sample.int(nrow(base), 500L), , drop = FALSE]
seed_path <- file.path(root, "data/Simulated_data_dml_seed.parquet")
if (!file.exists(seed_path)) write_parquet(base, seed_path)
sim <- base[sample.int(nrow(base), 100000L, replace=TRUE), , drop=FALSE]
sim$ID <- seq_len(nrow(sim))
required <- unique(c("age","sex","race","BMI","total_number_of_COVID_vaccine_doses","poverty_status","health_insurance",names(sim)[79:102],names(sim)[105:113]))
for (nm in required[required %in% names(sim)]) { x<-sim[[nm]]; if(anyNA(x)){ if(is.numeric(x)){v<-median(x,na.rm=TRUE);if(!is.finite(v))v<-0}else{z<-x[!is.na(x)];v<-if(length(z))z[1] else 0}; x[is.na(x)]<-v;sim[[nm]]<-x } }
covid <- as.Date(sim$COVID_19); covid[is.na(covid)] <- as.Date("2021-01-01") + sample(0:1095,sum(is.na(covid)),replace=TRUE); sim$COVID_19<-covid
lc_rows<-sample.int(nrow(sim),20000L);sim$Long_COVID<-as.Date(NA);sim$Long_COVID[lc_rows]<-sim$COVID_19[lc_rows]+sample(91:365,length(lc_rows),replace=TRUE);sim$Long_COVID_binary<-as.integer(!is.na(sim$Long_COVID))
targets<-c("Hypertension","Diabetes_without_complications","Depression","Chronic_lung_disease","Diabetes_with_complications","Kidney_disease","Coronary_artery_disease","Heart_failure","Liver_disease","Cerebrovascular_disease","Peripheral_vascular_disease","Myocardial_infarction","Cardiomyopathy","Dementia","Malignant_tumor","Benign_tumor")
for(nm in targets){ if(nm%in%names(sim)){e<-sample.int(nrow(sim),25000L);x<-rep(as.Date(NA),nrow(sim));x[e]<-sim$COVID_19[e]-sample(30:720,length(e),replace=TRUE);sim[[nm]]<-x} }
sim$Diabetes<-as.integer(!is.na(sim$Diabetes_without_complications)|!is.na(sim$Diabetes_with_complications));sim$Cardiovascular<-as.integer(rowSums(!is.na(sim[c("Cardiomyopathy","Coronary_artery_disease","Heart_failure","Myocardial_infarction","Peripheral_vascular_disease")]))>0);sim$Neuro_psych_disorders<-as.integer(rowSums(!is.na(sim[c("Cerebrovascular_disease","Dementia","Depression")]))>0);sim$Cancer<-as.integer(!is.na(sim$Malignant_tumor)|!is.na(sim$Benign_tumor));sim$Any_comorbidities<-as.integer(rowSums(!is.na(sim[targets]))>0);if("death_date"%in%names(sim))sim$death_date<-as.Date(NA);sim$Dead<-0L
write_parquet(sim,input);write_parquet(sim,file.path(root,"data/Simulated_data_dml_100k.parquet"));cat("Wrote",nrow(sim),"rows and",ncol(sim),"columns\n")

library(tidyverse); library("DECIPHER")




dist_vir <- readDNAStringSet("data/ASVs.fa") %>% 
	OrientNucleotides() %>% 
	DistanceMatrix() %>%
	as_tibble(rownames = "vir_A") %>%
	pivot_longer(names_to = "vir_B",
				 cols = -1,
				 values_to = "distance")

covar_virbact <- read_csv("data/virbact.cov.csv") %>% 
	select(-X1) %>% 
	rename(weight = "cor")

covar_virvir <- read_csv("data/virvir.covar.csv") %>% 
	select(-X1)

vir_all <- covar_virvir %>% 
	mutate(vir_A = pmin(from, to),
		   vir_B = pmax(from, to)) %>%
	left_join(dist_vir) %>% 
	select(vir_A, vir_B, distance, cov = weight)

vir_all %>% write_csv("export/vir_all.csv")






# test <- cor_virbact %>% 
# 	mutate(vir = from %>% str_extract("(?<=vir_ASV_)\\d*"),
# 		   bact = to %>% str_extract("(?<=ASV_)\\d*"))

# test$vir %>% 
# 	map_lgl(~{. %in% test$bact}) %>% 
# 	sum()

library(tidyverse); library("DECIPHER")

dist_vir <- readDNAStringSet("data/ASVs.fa") %>% 
	OrientNucleotides() %>% 
	DistanceMatrix() %>%
	as_tibble(rownames = "vir_A") %>%
	pivot_longer(names_to = "vir_B",
				 cols = -1,
				 values_to = "distance") %>% 
	mutate(across(.cols = 1:2,
				  ~ .x %>% 
				  	str_remove("ASV_") %>% 
				  	as.integer()))

covar_virvir <- read_csv("data/virvir.covar.csv") %>% 
	select(-X1)

vir_all <- covar_virvir %>% 
	mutate(vir_A = pmin(from, to),
		   vir_B = pmax(from, to)) %>%
	mutate(across(-weight, ~str_remove(.x, "ASV_") %>% as.integer())) %>% 
	left_join(dist_vir) %>% 
	select(vir_A, vir_B, distance, cov = weight)

vir_all %>% write_csv("export/vir_all.csv")

covar_virbact <- read_csv("data/virbact.cov.csv") %>% 
	select(-X1) %>% 
	rename(weight = "covar_vir_bact", from = "vir", to = "bact") %>% 
	mutate(vir = vir %>% str_remove("vir_ASV_") %>% as.integer())

vir_cov_dist <- covar_virbact %>% 
	group_by(bact) %>%
	# Restrict the analysis to bact for which we have at least 7 vir covariances
	filter(n() >= 7) %>% 
	group_modify(~{
		print(.y$bact)
		data <- .x %>% 
			# Generate pairs
			tidyr::expand(v1 = vir,v2 = vir) %>% 
			# Eliminate duplicates
			filter(v1 < v2) %>% 
			# Reattach covar with current bact
			left_join(.x, by = c("v1" = "vir")) %>% 
			left_join(.x, by = c("v2" = "vir"), suffix = c("1", "2")) %>% 
			# Attach distance
			left_join(dist_vir, by = c("v1" = "vir_A", "v2" = "vir_B")) %>% 
			mutate(abs_diff = abs(covar_vir_bact1 - covar_vir_bact2))
		
		res <- cor.test(x = data$distance, y = data$abs_diff)
		
		tibble(n_vir = nrow(.x),
			   n_pairs = nrow(data),
			   estimate = res$estimate,
			   p_value = res$p.value,
			   data = list(data))
	}) %>% 
	ungroup()

save.image("vir_cov_done.RData")
vir_cov_dist %>% filter(p_value < .05)

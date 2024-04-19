# Space use

load(file = glue("out/intra_ind_models/size_intra_ind_int_mod_2023-11-20.rdata"))

(n_ind <- length(unique(out$data$ind_f)))
(n <- nrow(out$data))
(n_sp <- length(unique(out$data$species)))


# Niche

load(file = glue("out/intra_ind_models/niche_intra_ind_int_mod_2023-11-20.rdata"))

(n_ind <- length(unique(out$data$ind_f)))
(n <- nrow(out$data))
(n_sp <- length(unique(out$data$scientificname)))

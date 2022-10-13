## code to prepare `DATASET` dataset goes here
m <- 2000 # size for RWE
niter <- 500 # number of replication
thres.psi  <-  sqrt(log(m)) # threshold for ECI psi
alltlocalpar <-round((seq(2,0,length.out=10)),2)

elastic_aipw_ses_l1_list <- lapply(alltlocalpar, function(tlocalpar){
  elastic_list <- sapply(1:niter, function(seed)
  {
    Data.list <- GenData(beta0 = c(0, 1, 1, 1), # for the mu0 function
                              psi0 = c(0, 1, 1), # for the contrast function
                              om.id = 0, es.id = 0,
                              n = 1e5, mean.x = 1,  # setup for the finite population
                              n.t = NULL, equal.weight = c(-1,-1),# for the RCT, use the default sample size
                              m = m, tlocalpar = tlocalpar, # for the RWE
                              seed = seed)
    elasticHTE(Data.list$RT, # RCT
               Data.list$RW, # RWE
               thres.psi = thres.psi,
               fixed = FALSE # adaptive selection strategy
    )
  })
  elastic_list
})

elastic_aipw_ses_l2_list <- lapply(alltlocalpar, function(tlocalpar){
  elastic_list <- sapply(1:niter, function(seed)
  {
    Data.list <- GenData(beta0 = c(0, 1, 1, 1), # for the mu0 function
                         psi0 = c(0, 1, 1), # for the contrast function
                         om.id = 0, es.id = 0,
                         n = 1e5, mean.x = 1,  # setup for the finite population
                         n.t = NULL, equal.weight = c(-2,-2),# for the RCT, use the default sample size
                         m = m, tlocalpar = tlocalpar, # for the RWE
                         seed = seed)
    elasticHTE(Data.list$RT, # RCT
               Data.list$RW, # RWE
               thres.psi = thres.psi,
               fixed = FALSE # adaptive selection strategy
    )
  })
  elastic_list
})

elastic_aipw_ses_l3_list <- lapply(alltlocalpar, function(tlocalpar){
  elastic_list <- sapply(1:niter, function(seed)
  {
    Data.list <- GenData(beta0 = c(0, 1, 1, 1), # for the mu0 function
                         psi0 = c(0, 1, 1), # for the contrast function
                         om.id = 0, es.id = 0,
                         n = 1e5, mean.x = 1,  # setup for the finite population
                         n.t = NULL, equal.weight = c(-3,-3),# for the RCT, use the default sample size
                         m = m, tlocalpar = tlocalpar, # for the RWE
                         seed = seed)
    elasticHTE(Data.list$RT, # RCT
               Data.list$RW, # RWE
               thres.psi = thres.psi,
               fixed = FALSE # adaptive selection strategy
    )
  })
  elastic_list
})

# write out the raw data to Rd
usethis::use_data(elastic_aipw_ses_l1_list,
                  overwrite = TRUE, compress = 'xz')

usethis::use_data(elastic_aipw_ses_l2_list,
                  overwrite = TRUE, compress = 'xz')

usethis::use_data(elastic_aipw_ses_l3_list,
                  overwrite = TRUE, compress = 'xz')

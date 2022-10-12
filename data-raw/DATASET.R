## code to prepare `DATASET` dataset goes here
m <- 2000 # size for RWE
niter <- 2000 # number of replication
thres.psi  <-  sqrt(log(m)) # threshold for ECI psi
alltlocalpar <-round((seq(2,0,length.out=10)),2)

elastic_psi111_lists <- lapply(alltlocalpar, function(tlocalpar){
  elastic_list <- sapply(1:niter, function(seed)
  {
    Data.list <- GenerateData(beta0 = c(0, 1, 1, 1), # for the mu0 function
                              psi0 = c(0, 1, 1), # for the contrast function
                              om.id = 0, es.id = 0,
                              n = 1e5, mean.x = 1,  # setup for the finite population
                              n.t = NULL, # for the RCT, use the default sample size
                              m = m, tlocalpar = tlocalpar, # for the RWE
                              seed = seed)
    elasticHTE(Data.list$RT, # RCT
               Data.list$RW, # RWE
               fixed = FALSE # adaptive selection strategy
    )
  })
  elastic_list
})

# write out the raw data to Rd
usethis::use_data(elastic_psi111_lists,
                  overwrite = TRUE)

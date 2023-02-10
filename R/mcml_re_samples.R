# gen_u_samples <- function(y,Xb,A,D,NN,warmup_iter=100,m=100){
# 
#   data <- list(
#     N = nrow(X),
#     Xb = Xb,
#     y = y,
#     m = nrow(NN),
#     A = A,
#     D= D
#   )
#   
#   res <- rstan::sampling(stanmodels[[fname]],
#                          data=datlist,
#                          chains=chains,
#                          iter = iter_warmup+iter_sampling,
#                          warmup = iter_warmup,
#                          cores = parallel_chains,
#                          refresh = 0)
# 
# 
#   if(!requireNamespace("cmdstanr")){
#     stop("cmdstanr not available")
#   } else {
#     file_type <- mcnr_family(family)
#     model_file <- system.file("stan",
#                               file_type$file,
#                               package = "glmmrMCML",
#                               mustWork = TRUE)
#     mod <- suppressMessages(cmdstanr::cmdstan_model(model_file))
#     data <- list(
#       N = nrow(X),
#       Q = ncol(Z),
#       Xb = drop(as.matrix(X)%*%beta),
#       Z = as.matrix(Z)%*%L,
#       y = y,
#       sigma = sigma,
#       type=as.numeric(file_type$type)
#     )
# 
#     capture.output(fit <- mod$sample(data = data,
#                                      chains = 1,
#                                      iter_warmup = warmup_iter,
#                                      iter_sampling = m,
#                                      refresh = 0),
#                    file=tempfile())
#     dsamps <- fit$draws("gamma",format = "matrix")
#     dsamps <- L%*%t(dsamps)
#     return(dsamps)
#   }
# 
# }



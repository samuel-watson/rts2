test_that("data functions work", {
  b1 <- sf::st_sf(sf::st_sfc(sf::st_polygon(list(cbind(c(0,3,3,0,0),c(0,0,3,3,0))))))
  npoints <- 30
  dp <- data.frame(y=3*sqrt(runif(npoints)),
                   x=3*sqrt(runif(npoints)),
                   date=paste0("2021-01-",sample(11:13,npoints,replace = TRUE)))


  # create a coarse grid over the area
  g1 <- grid$new(poly = b1,cellsize = 0.1)

  expect_s3_class(g1$grid_data,"sf")

  # create the points sf object
  dp <- create_points(dp,pos_vars = c('y','x'),t_var='date')

  expect_s3_class(dp,"sf")

  #create a random covariate over the area to act as population density
  cov1 <- grid$new(poly = b1,cellsize = 0.12)
  cov1$grid_data$cov <- runif(nrow(cov1$grid_data))

  # map the population density to the grid
  g1$add_covariates(cov1$grid_data,
                    zcols="cov",
                    verbose = FALSE)

  expect_s3_class(g1$grid_data,"sf")
  expect_true("cov"%in%colnames(g1$grid_data))

  #aggregate the points to the grid
  g1$points_to_grid(dp, laglength=3)
  expect_s3_class(g1$grid_data,"sf")
  expect_true("t3"%in%colnames(g1$grid_data))
})

test_that("sampler and analysis works", {
  b1 <- sf::st_sf(sf::st_sfc(sf::st_polygon(list(cbind(c(0,3,3,0,0),c(0,0,3,3,0))))))
  npoints <- 30
  dp <- data.frame(y=3*sqrt(runif(npoints)),
                   x=3*sqrt(runif(npoints)),
                   date=paste0("2021-01-",sample(11:13,npoints,replace = TRUE)))


  # create a coarse grid over the area
  g1 <- grid$new(poly = b1,cellsize = 0.1)

  # create the points sf object
  dp <- create_points(dp,pos_vars = c('y','x'),t_var='date')

  #create a random covariate over the area to act as population density
  cov1 <- grid$new(poly = b1,cellsize = 0.12)
  cov1$grid_data$cov <- runif(nrow(cov1$grid_data))

  # map the population density to the grid
  g1$add_covariates(cov1$grid_data,
                    zcols="cov",
                    verbose = FALSE)

  expect_s3_class(g1$grid_data,"sf")

  #aggregate the points to the grid
  g1$points_to_grid(dp, laglength=3)

  g1$priors <- list(
    prior_lscale=c(0,0.5),
    prior_var=c(0,0.5),
    prior_linpred_mean=c(0),
    prior_linpred_sd=c(5)
  )

  # run the model
  res <- g1$lgcp_fit(popdens="cov",
                    priors = priors,
                    iter_warmup = 10,
                    iter_sampling = 10,
                    verbose = FALSE,
                    chains = 1,
                    parallel_chains = 1)

  expect_s4_class(res,"stanfit")

  # extract the predictions
  g1$extract_preds(res,
                   type=c("rr"),
                   popdens="cov")

  expect_s3_class(g1$grid_data,"sf")
  expect_true("rr"%in%colnames(g1$grid_data))

})





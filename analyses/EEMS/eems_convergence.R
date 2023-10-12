#### eems plot ####
library("reemsplots2")
mcmcpaths = paste("chain",1:3,sep="")
plots <-  lapply(mcmcpaths, function(x) make_eems_plots(x,longlat=TRUE) )
names(plots) <- paste("run",1:3,sep="")



#### check the convergence ####
library(dplyr)
library('coda')
library(ggplot2)

### define function here ###
plot_log_posterior <- function(mcmcpath) {
  message("Generate posterior probability trace. ",
          "See plots$pilog01.")
  rleid <- function(x) {
    r <- rle(x)
    rep(seq_along(r$lengths), r$lengths)
  }
  pl_df <- NULL
  for (path in mcmcpath) {
    pl <- read_matrix(file.path(path, "mcmcpilogl.txt"))
    pl_df <- bind_rows(pl_df, as_data_frame(pl) %>% mutate(path))
  }
  pl_df <- pl_df %>%
    setNames(c("pi", "logl", "path")) %>%
    mutate(mcmcpath = factor(rleid(path))) %>%
    group_by(mcmcpath) %>%
    mutate(iter = row_number(), pilogl = pi + logl)
  ggplot(pl_df, aes(x = iter, y = pilogl, color = mcmcpath)) +
    geom_path() +
    labs(x = "MCMC iteration  (after burn-in and thinning)",
         y = "log posterior",
         title = "Have the MCMC chains converged?",
         subtitle = "If not, restart EEMS and/or increase numMCMCIter") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank())
}

read_matrix <- function(file, ncol = 2) {
  stopifnot(file.exists(file))
  matrix(scan(file, what = numeric(), quiet = TRUE),
         ncol = ncol, byrow = TRUE)
}

rleid <- function(x) {
  r <- rle(x)
  rep(seq_along(r$lengths), r$lengths)
}

### process data here ###
pl_df <- NULL

mcmcpaths = paste("chain",1:3,sep="")

for (path in mcmcpaths) {
  pl <- read_matrix(file.path(path, "mcmcpilogl.txt"))
  pl_df <- bind_rows(pl_df, as_data_frame(pl) %>% mutate(path))
}

pl_df <- pl_df %>%
  setNames(c("pi", "logl", "path")) %>%
  mutate(mcmcpath = factor(rleid(path))) %>%
  group_by(mcmcpath) %>%
  mutate(iter = row_number(), pilogl = pi + logl)

### plot ###
ggplot(pl_df, aes(x = iter, y = pilogl, color = mcmcpath)) +
  geom_path() +
  labs(x = "MCMC iteration  (after burn-in and thinning)",
       y = "log posterior",
       title = "Have the MCMC chains converged?",
       subtitle = "If not, restart EEMS and/or increase numMCMCIter") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())


### prepare and run the gelman test ###
step=1000 ### number of observations in one chain, assuming the same number across chains
v1<-as.vector(pl_df$pilogl)

n_chain <- length(v1)/step
start_pos = 1
end_pos = step
for (i in 1:n_chain) {
  assign( paste("chain",i,sep=""),  mcmc(v1[start_pos:end_pos],start=1,end=step,thin=1) )
  start_pos <- start_pos + step
  end_pos <- step *(i+1)
}

combinedchains = mcmc.list( mget(paste("chain",1:n_chain,sep="")) )
plot(combinedchains)
gelman.diag(combinedchains)
gelman.plot(combinedchains)

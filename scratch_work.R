sigma_total <- 4.475 # prop of total: 100%
sigma_error <- 2.12 # prop of total: 47%
sigma_block <- 1.78 # prop of total: 40%
sigma_int <- 0.575 # prop of total: 13%

delta <- 0.40
m <- 4
n <- 10


N <- 2 * m * n
msw <- sigma_error
msb <- sigma_error + 2 * n * sigma_block
msab <- sigma_error + n * sigma_int
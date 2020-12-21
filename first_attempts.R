# need to know:
# msw
# m
# n
# mswb
# mswab
# 
# we're interested in comparing y/s ("d1") (traditional, bottom of page 7) to y/s * sqrt(b) (also bot pag 7) ("d2")
# our variance is the top of page 8 ("v2")
# we think the traditional variance method is sigma_squared_t (sum of all sigma squared) over N - 2 ("v1")
# need m, n, Y, S, MSW, MSB, MSAB
# probably a simulation table
# delta is the mean difference (y) over square root of sigmaTotal, right now we could just keep it fixed at .4 over s
# the delta squared in the variance is just d2 squared
# in the excel spreadsheet did term 1 and term 2
# try to have a basic table by tomorrow
# try varying number in each block and interaction proportion. fix total and error, vary block and interaction proportion
# meeting is tomorrow at 3

msw <- 10; m <- 4; n <- 100; mswb <- 20; mswab <- 30

sigma_squared_e <- msw
sigma_fourth_e <- (2*m*n - 2*m)*msw^2/(2*m*n - 2*m + 2)
sigma_squared_b <- (mswb - msw)/(2*n)
sigma_squared_ab <- (mswab - msw)/n
sigma_squared_t <- sigma_squared_e + sigma_squared_ab +
  sigma_squared_b

e_v <- sigma_squared_e + ( (2*n*(m - 1)*sigma_squared_b + 
                              n*(m - 1)*sigma_squared_ab)/(2*m*n - 2) )

v_v <- ((2*(2*m*n - 2*m)*sigma_squared_e) + 
          (2*(m - 1)*(sigma_squared_e + 2*n*sigma_squared_b)^2) + 
     (2*(m - 1)*(sigma_squared_e + n*sigma_squared_ab)^2))/((2*m*n - 2)^2)

b <- ((2*m*n - 2*m)*msw + (m - 1)*(mswb + mswab))/((2*m*n - 2)*sigma_squared_t)

c <- (1/((2*m*n - 2)^2*(sigma_squared_t)^2)) * (((2*m*n - 2*m)^2*msw^2)/(2*m*n - 2*m + 2) + ((m - 1)^2*mswb^2)/(m + 1) + ((m - 1)^2*mswab^2)/(m + 1))

a <- (2*n + mswab)/(mswb + 2*mswab + (2*n - 3)*msw)

cat("a: ", a, "\nb: ", b, "\nc: ", c, "\ne_v: ", e_v, "\nv_v: ", v_v)

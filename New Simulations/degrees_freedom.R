m=2
n=100
N=2*m*n
var_e<- 0.3
var_inter<- 0
var_block <- 0.7

df_num<- ((N - 2 * m) * var_e + (m - 1) * (var_e + 2 * n * var_block) + (m - 1) * (var_e + n * var_inter))^2

df_den <- ((N - 2 * m) * var_e^2 + (m - 1) * (var_e + 2 * n * var_block)^2 + (m - 1) * (var_e + n * var_inter)^2)

df_total<- df_num/df_den
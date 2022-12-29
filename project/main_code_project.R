# Projet:

rm(list = ls())

# Exercice 1:

# Méthode de Monte Carlo n°1:

# Question 1:

# Soit l'inverse généralisé de la fonction de répartition de Y,

fct_inv_Y <- function(x) {
  return((sqrt(8 * 4 * x) - 17) * (x < (1 / 2)) + (-13 - 2 * log(2 * (1 - x))) * (x > (1 / 2)))
}

# Pour simuler la var Y,

ech_Y <- function(n) {
  return(fct_inv_Y(runif(n)))
}

# Pour simuler X,

ech_X <- function(n) {
  return(rgamma(n, shape = 3, scale = 2))
}

# Pour simuler maintenant Z=(X,Y),

simu <- function(n) {
  return(rbind(ech_X(n), ech_Y(n)))
}

# On pose n=1000.
# On gardera cette affection par la suite, sauf mention contraire explicite.

n <- 1000

# On simule Y et X.

Y <- ech_Y(n)
X <- ech_X(n)

# Question 2:

# On utilise viridisLite pour avoir une représentation graphique plus accessible.
library(viridisLite)
c_pal <- viridis(2)

# On définit la fonction g densité de Y
g <- function(x) {
  return(((x + 17) / 16) * (x > (-17) &
                              x < (-13)) + (exp(-(x + 13) / 2) / 4) * (x > (-13)))
}

par(mfrow = c(1, 1))

# On créé un histogramme pour représenter les réalisations empiriques de Y
hist(
  Y,
  freq = FALSE,
  col = 'grey60',
  border = 'grey70',
  breaks = 30,
  main = 'Histogramme de Y',
  ylim = c(0, 0.3),
  ylab = 'Frequences'
)
# On représente les réalisations théoriques de Y
lines(y <- seq(min(Y), max(Y), 0.01), g(y), col = c_pal[2], lwd = 3)

legend(
  'topright',
  "densité de Y",
  col = c_pal[2],
  lwd = 3,
  box.lty = 0,
  bg = 'gray95',
  inset = .05
)

#Question 3:

# On définit la fonction h utilisée dans l'estimateur MC classique
h <- function(x_1, x_2) {
  return(1 * ((x_1 + x_2) > 0))
}

# Question 4:

# On peut calculer h(X,Y) à partir des réalisations empiriques de X et Y
pn_hat_simu <- function(n) {
  Y <- ech_Y(n)
  X <- ech_X(n)
  return(h(X, Y))
}

pn_hat_ech <- pn_hat_simu(n)
# On peut calculer le résultat de l'estimateur MC classique
pn_hat = mean(pn_hat_ech)

# L'estimateur de Monte Carlo étant sans biais, l'erreur quadratique correspond
# à la variance de l'estimateur

MSE_4 <- var(pn_hat_ech) / n

# Estimateur de la variance de h(X,Y):

estimateur_variance_4 <-
  mean((pn_hat_ech - pn_hat) ** 2) * (n / (n - 1))

# L'intervalle de confiance asymptotique:

Intervalle <-
  pn_hat + c(-1, 1) * qnorm(0.975) * sqrt(estimateur_variance_4 / n)

# Question 5:

# Vérifions si notre estimateur est asymptotiquement normal:

p <- ppoints(500)
# On va simuler 10000 fois l'estimateur pn_hat pour ensuite montrer qu'il suit asymptotiquement une loi normale
echantillon_pn <- colMeans(matrix(pn_hat_simu(1000 * 10000), nrow = 1000))

qqnorm(
  quantile(echantillon_pn, p),
  main = "Normal Q-Q Plot",
  xlab = "Quantiles théoriques",
  ylab = "Quantiles de l'échantillon des estimateurs pn",
  col = c_pal[2]
)
qqline(echantillon_pn, col = 'tomato3', lwd = 2)

# Comme nous pouvons le voir, les points sont alignés sur le graphique, donc notre estimateur est asymptotiquement normal.

# Question 7:

# On peut définir la fonction h associée à l'estimateur MC de la question 6
h_7 <- function(n) {
  return(1 - pgamma(-ech_Y(n), shape = 3, scale = 2))
}

h_tilde <- function(x) {
  return(1 - pgamma(-x, shape = 3, scale = 2))
}

# On calcule le résultat de ce nouvel estimateur MC
ech_7 <- h_7(n)
delta_n_hat <- mean(ech_7)

# Erreur quadratique moyenne:
# L'estimateur étant non biaisé, cette dernière est égale à la variance
MSE_7 <- var(ech_7) / n

# Estimateur de la variance:

estimateur_variance_7 <-
  mean((ech_7 - delta_n_hat) ** 2) * (n / (n - 1))

# Intervalle de confiance asymptotique de niveau 95%:

Intervalle_7 <-
  delta_n_hat + c(-1, 1) * qnorm(0.975) * sqrt(estimateur_variance_7 / n)

# cMéthode de Monte carlo n°2:

# Question 10

# On implémente une fonction qui renvoie l'échantillon permettant de créer
# ensuite l'estimateur de la variable antithétique associé à la fonction A

A <- function(x) {
  return(1 - x)
}

ech_h_antithetique <- function(n) {
  u <- runif(n)
  return((h_7(fct_inv_Y(u)) + h_7(fct_inv_Y(A(
    u
  )))) / 2)
}

ech_h_ant <- ech_h_antithetique(n)

# Résultat de l'estimateur MC de la variable antithétique
delta_A_n_hat <- mean(ech_h_ant)

# Erreur quadratique moyenne associée à ce nouvel estimateur
MSE_anthitetique <- var(ech_h_ant) / n

# Estimateur de la variance:
estimateur_variance_anth <-
  mean((ech_h_ant - delta_A_n_hat) ** 2) * (n / (n - 1))

# Intervalle de confiance associé à ce nouvel estimateur
Intervalle_10 <-
  delta_A_n_hat + c(-1, 1) * qnorm(0.975) * sqrt(estimateur_variance_anth /
                                                   n)

#Méthode de Monte Carlo n°3:

# Question 11
# On commence par regarder quelle corrélation au carré est la plus proche de 1:

q_epsilon <- qgamma(0.6, shape = 3, scale = 2)

h_0_2 <- function(x) {
  return(1 * (x > -q_epsilon))
}

cor_1 <- cor(h_7(Y), Y) ** 2 #h_0_1 est la fonction identité
cor_2 <- cor(h_7(Y), h_0_2(Y)) ** 2

# On observe qu'il s'agit de cor_1: on utilise donc la fonction h_0_1, qui
# est la fonction identité

# Question 12:

# Pour déterminer le coefficient optimal b*, on utilise la stratégie n°1 du cours.

# On commence par choisir pour quel l on implémente la stratégie.
# Pour cela, on regarde à partir de quelle valeur de l b_hat défini dans le cours stagne.


choix_l_opti <- function(k) {
  y <- ech_Y(k)
  return(cumsum((y + (38 / 3)) * (h_tilde(y) - mean(h_tilde(
    y
  )))) / (sum((y + (
    38 / 3
  )) * (y + (
    38 / 3
  )))))
}

plot(choix_l_opti(27),
     col = 'tomato3',
     type = 'b',
     main = 'valeur de b* en fonction de l')

ech_h_controle <- function(n, l = 16) {
  y <- ech_Y(n)
  h <- h_tilde(y)
  b_l_opti_hat <-
    sum((y[1:l] + (38 / 3)) * (h[1:l] - mean(h[1:l]))) / (sum((y[1:l] + (38 /
                                                                           3)) * (y[1:l] + (38 / 3))))
  return(h[n - l:n] - b_l_opti_hat * (y[n - l:n] + (38 / 3)))
}

ech_h_cont <- ech_h_controle(n)

# Résultat de l'estimateur MC de la variable de contrôle
delta_cont_n_hat <- mean(ech_h_cont)

# Erreur quadratique moyenne associée à cet estimateur sans biais
MSE_cont <- var(ech_h_cont) / n

# Estimateur de la variance:
estimateur_variance_cont <-
  mean((ech_h_cont - delta_cont_n_hat) ** 2) * (n / (n - 1))

# Intervalle de confiance associé à cet estimateur de la variable de contrôle
Intervalle_12 <-
  delta_cont_n_hat + c(-1, 1) * qnorm(0.975) * sqrt(estimateur_variance_cont /
                                                      n)

# Question 17:

ech_h_17 <- function(n) {
  sample <- ech_Y(n)
  a <- exp((q_epsilon - 13) / 2)
  sigma <-
    matrix(c(
      47 / 9,
      a - q_epsilon * a / 2,
      a - q_epsilon * a / 2,
      a / 2 - exp(q_epsilon - 13) / 4
    ), nrow = 2)
  fdr <- pgamma(-sample, scale = 2, shape = 3)
  C <- c(cov(fdr, sample), cov(fdr, 1 * (sample > q_epsilon)))
  beta_opti <- solve(sigma) %*% C
  return(1 - fdr - beta_opti[1] * (sample + 38 / 3) - beta_opti[2] *
           (h_0_2(sample) - a / 2))
}

ech_h_question_17 <- ech_h_17(n)

# Résultat de l'estimateur MC de la question 16
delta_n_beta_hat <- mean(ech_h_question_17)

MSE_beta <-
  var(ech_h_question_17) / n

# Estimateur de la variance:
estimateur_variance_17 <-
  mean((ech_h_question_17 - delta_n_beta_hat) ** 2) * (n / (n - 1))

Intervalle_17 <-
  delta_n_beta_hat + c(-1, 1) * qnorm(0.975) * sqrt(estimateur_variance_17 /
                                                      n)

# Question 20:
ech_h_18 <- function(n) {
  k <- 1:n
  return(h_7(fct_inv_Y((k - 1 + runif(
    n
  )) / n)))
}

ech_h_question_18 <- ech_h_18(n)
delta_n_hat_18 <- mean(ech_h_question_18)
MSE_18 <- var(ech_h_question_18) / n
estimateur_variance_18 <-
  mean((ech_h_question_18 - delta_n_hat_18) ** 2) * (n / (n - 1))

Intervalle_18 <-
  delta_n_hat_18 + c(-1, 1) * qnorm(0.975) * sqrt(estimateur_variance_18 /
                                                    n)

# Question 21:

# On calcule ici le coût pour construire chaque échantillon des différentes méthodes,
# car il suffit ensuite de déterminer la moyenne de chaque échantillon pour avoir ensuite
# l'estimateur associé, donc le coût restant pour déterminer chaque estimateur est identique
# et ne va donc pas modifier les efficacités relatives si on les comptabilise.

library(microbenchmark)
microbenchmark(pn_hat_simu(n),
               h_7(n),
               ech_h_antithetique(n),
               ech_h_controle(n),
               ech_h_17(n))

c_pn <- mean(microbenchmark(pn_hat_simu(n))$time)
c_7 <- mean(microbenchmark(h_7(n))$time)
c_anth <- mean(microbenchmark(ech_h_antithetique(n))$time)
c_cont <- mean(microbenchmark(ech_h_controle(n))$time)
c_17 <- mean(microbenchmark(ech_h_17(n))$time)
c_18 <- mean(microbenchmark(ech_h_18(n))$time)

outer(
  c(
    c_pn * estimateur_variance_4,
    c_7 * estimateur_variance_7,
    c_anth * estimateur_variance_anth,
    c_cont * estimateur_variance_cont,
    c_17 * estimateur_variance_17,
    c_18 * estimateur_variance_18
  ),
  c(
    c_pn * estimateur_variance_4,
    c_7 * estimateur_variance_7,
    c_anth * estimateur_variance_anth,
    c_cont * estimateur_variance_cont,
    c_17 * estimateur_variance_17,
    c_18 * estimateur_variance_18
  ),
  '/'
)

# Exercice 2

# On a la densité de psi
psi <- function(x, m, k) {
  return((factorial(m) / (factorial(k - 1) * factorial(m - k))) * G(x) ^ (k -
                                                                            1) * (1 - G(x)) ^ (m - k) * g(x))
}

# On a la fonction de répartition G de l'exercice 1
G <- function(x) {
  return((x > (-17) &
            x < (-13)) * ((1 / 16) * (x / sqrt(2) + 17 / sqrt(2)) ^ 2) + (x >= (-13)) *
           (1 - (1 / 2) * exp((-x - 13) / 2)))
}

# On peut mettre en place un algorithme du rejet pour simuler des réalisations suivant Psi
simulation_rejet_psi <- function(n, m, k) {
  ech <- c()
  while (length(ech) < n) {
    u <- runif(n)
    y <- ech_Y(n)
    ech <-
      c(ech, y[u < ((m - 1) / (k - 1)) ^ (k - 1) * ((m - 1) / (m - k)) ^ (m -
                                                                            k) * G(y) ^ (k - 1) * (1 - G(y)) ^ (m - k)])
  }
  return(ech[1:n])
}

# Question 2
# On définit un échantillon de référence pour la représentation graphique
simulations_exo_2 <- simulation_rejet_psi(1000, 15, 7)

# On représente les réalisations empiriques de notre algorithme du rejet
hist(
  simulations_exo_2,
  freq = FALSE,
  col = 'grey60',
  border = 'grey70',
  breaks = 30,
  main = 'Histogramme de Psi',
  ylim = c(0, 1),
  ylab = 'Frequences'
)

# On représente les réalisations théoriques de Psi
lines(y <-
        seq(min(simulations_exo_2), max(simulations_exo_2), 0.01),
      psi(y, 15, 7),
      col = c_pal[2],
      lwd = 3)

legend(
  'topright',
  "densité de Psi",
  col = c_pal[2],
  lwd = 3,
  box.lty = 0,
  bg = 'gray95',
  inset = .05,
  cex = 0.7
)


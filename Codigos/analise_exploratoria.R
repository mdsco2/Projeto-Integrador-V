#install.packages("moments")
#install.packages("evd")
#install.packages("zoo")
#install.packages("dgof")
#install.packages('nortest')
#install.packages("DescTools")
#install.packages('EnvStats')
#install.packages("distr6")
library('EnvStats')
library("moments")
library("evd")
library("zoo")
library("dgof")
library("tidyverse")
library("ggplot2")
library("ggpubr")
library("nortest")
library("DescTools")
library("EnvStats")
library("distr6")


options(scipen=999)

getwd()

df <- read.csv("dados_versao_final.csv", sep=";")
df

names(df)[1] <- "data_medicao"
names(df)[2] <- "temp_max"
df

sapply(df, class)

summary(df)

df[df == "null"] <- NA

missing <- df[rowSums(is.na(df)) > 0,]

missing

df$temp_max <- as.numeric(as.character(df$temp_max))
df

df$data_medicao <- as.Date(df$data_medicao, "%d/%m/%Y")
df

summary(df)

df <- na.omit(df)

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
getmode(df$temp_max)

skewness(df$temp_max)
kurtosis(df$temp_max)

mean(df$temp_max)
var(df$temp_max)

# Gráfico de Densidade
plot(density(df$temp_max), xlab='Temperatura Máxima', ylab='Densidade', main='Gráfico de Densidade das Temperaturas Máximas de Brasília (1980 a 2021)')
abline(v=mean(df$temp_max), col='red') 

var.test(df$temp_max>26.8, df$temp_max<26.8)

boxplot(df$temp_max, col = "#69b3a2", main = "Temperaturas Máximas Brasília (1980 a 2021)", 
        horizontal = TRUE)

sturges = hist(df$temp_max, col = 'skyblue3', breaks = "Sturges", main = "Temperaturas Máximas Brasília (Sturges)",
     ylab = "Frequência", xlab = "Temperatura") 

fd = hist(df$temp_max, col = 'skyblue3', breaks = "fd", main = "Temperaturas Máximas Brasília (Freedman and Diaconis)",
     ylab = "Frequência", xlab = "Temperatura") 

scott = hist(df$temp_max, col = 'skyblue3', breaks = "scott", main = "Temperaturas Máximas Brasília (Scott)",
     ylab = "Frequência", xlab = "Temperatura") 

nclass.FD(df$temp_max)

nclass.scott(df$temp_max)

nclass.Sturges(df$temp_max)

n = nrow(df)
data <- rnorm(n)
doane_breaks <- (1 + log(n) + log(1 + mean((data - (mean(data))^4)/mean((data - mean(data))^2)^2) * sqrt(n / 6.)))
doane_breaks  

rice_breaks <- (2 * (n^(1/3)))
rice_breaks

larson_breaks <- 1 + (2.2 * log10(n))
larson_breaks

doane = hist(df$temp_max, col = 'skyblue3', breaks = doane_breaks, main = "Temperaturas Máximas Brasília (Doane)",
     ylab = "Frequência", xlab = "Temperatura") 

rice = hist(df$temp_max, col = 'skyblue3', breaks = rice_breaks, main = "Temperaturas Máximas Brasília (Rice)",
     ylab = "Frequência", xlab = "Temperatura")

larson = hist(df$temp_max, col = 'skyblue3', breaks = larson_breaks, main = "Temperaturas Máximas Brasília (Larson)",
     ylab = "Frequência", xlab = "Temperatura")

par(mfrow=c(2,3))

sturges = hist(df$temp_max, col = 'skyblue3', breaks = "Sturges", main = "Sturges",
               ylab = "", xlab = "")

doane = hist(df$temp_max, col = 'skyblue3', breaks = 1 + log(n) + log(1 + mean((data - (mean(data))^4)/mean((data - mean(data))^2)^2) * sqrt(n / 6.)), main = "Doane",
             ylab = "", xlab = "") 

larson = hist(df$temp_max, col = 'skyblue3', breaks = 1 + (2.2 * log10(n)), main ="Larson",
              ylab = "", xlab = "")

scott = hist(df$temp_max, col = 'skyblue3', breaks = "scott", main = "Scott",
             ylab = "", xlab = "") 

rice = hist(df$temp_max, col = 'skyblue3', breaks = 2 * (n^(1/3)), main = "Rice",
            ylab = "", xlab = "")

fd = hist(df$temp_max, col = 'skyblue3', breaks = "fd", main = "Freedman and Diaconis",
          ylab = "", xlab = "")

#Distribuicoes

#Normal:
x = df$temp_max
hist(x, breaks = "fd", freq = FALSE)
curve(dnorm(x,mean(df$temp_max),sd(df$temp_max)), col = "dark green", lwd = 2, add=TRUE)

#Gammma: important due to its relation to exponential and normal distributions. Uses of the gamma function:
#estimation for insurance loss, waiting time. Parameters: scale and shape.
theta <- var(df$temp_max)/mean(df$temp_max)
k <- mean(df$temp_max)/theta
a_gama <- k
b_gama <- 1/theta
x = df$temp_max
hist(x, breaks = "fd", freq = FALSE)
curve(dgamma(x,a_gama,b_gama), col = "red", lwd = 2, add=TRUE)

#Gumbel: used to model the maxima or minima of a number of samples from various different distributions, such as
#maximum level of a river over multiple years, extreme earthquake events, etc. Does not use parameters derived from
#the descriptive statistics. The parameters need to be estimated with Maximum likelihood parameters estimates.
dgumbel <- function(x, a, b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
pgumbel <- function(q, a, b) exp(-exp((a-q)/b))
qgumbel <- function(p, a, b) a-b*log(-log(p))

b_gumbel <- sqrt((6*var(df$temp_max)/pi^2))
a_gumbel <- -b_gumbel*-digamma(1)+mean(df$temp_max)
dgumbel <- function(x, a_gumbel, b_gumbel) 1/b_gumbel*exp((a_gumbel-x)/b_gumbel)*exp(-exp((a_gumbel-x)/b_gumbel))
x = df$temp_max
hist(x, breaks = "fd", freq = FALSE)
curve(dgumbel(x, a_gumbel, b_gumbel), col = "blue", lwd = 2, add=TRUE)

#lognormal: if we apply a log function to a dataset and get a normal distribution, then the data is 
#log normally distributed. Used for income distributions, hospitalization days. Parameters: mean and standard
#deviation on the log scale.
variancia <- log((var(df$temp_max)/exp(2*log(mean(df$temp_max))))+1)
media <- (2*log(mean(df$temp_max))-variancia)/2
x = df$temp_max
hist(x, breaks = "fd", freq = FALSE)
curve(dlnorm(x,media,sqrt(variancia)), col = "black", lwd = 2, add=TRUE)

hist(x, breaks = "fd", freq = FALSE, main = "Possíveis Distribuições para a Temperatura Máxima em Brasília (1980-2021)", ylab = "Densidade", xlab = "Temperatura (°C)")
curve(dgamma(x,a_gama,b_gama), col = "red", lwd = 2, add=TRUE)
curve(dgumbel(x,a_gumbel,b_gumbel), col = "blue", lwd = 2, add=TRUE)
curve(dnorm(x,mean(df$temp_max),sd(df$temp_max)), col = "dark green", lwd = 2, add=TRUE)
curve(dlnorm(x,media,sqrt(variancia)), col = "black", lwd = 2, add=TRUE)
legend(x = "topright",legend = c("Gama", "Gumbel", "Normal", "Log_normal"), col = c("red","blue","dark green","black"), lwd = 2)

# Kolmogorov-Smirnov

ks.test(df$temp_max, "pgamma", a_gama, b_gama, alternative = "greater")
ks.test(df$temp_max, "pgumbel", a_gumbel, b_gumbel, alternative = "greater")
ks.test(df$temp_max, "pnorm", mean(df$temp_max), sd(df$temp_max), alternative = "greater")
ks.test(df$temp_max, "plnorm", media, sqrt(variancia), alternative = "greater")

x_emp <- ecdf(df$temp_max)
x_norm <- ecdf(rnorm(5000,mean(df$temp_max),sd(df$temp_max)))
x_gama <- ecdf(rgamma(5000, a_gama,b_gama))
x_gumbel <- ecdf(rgumbel(5000,a_gumbel,b_gumbel))
x_lnorm <- ecdf(rlnorm(5000,media,sqrt(variancia)))
plot(x_emp, verticals = TRUE, do.points = FALSE, col="blue")
plot(x_norm, verticals = TRUE, do.points = FALSE, col="red", add=TRUE)
plot(x_gama, verticals = TRUE, do.points = FALSE, col="dark green", add=TRUE)
plot(x_gumbel, verticals = TRUE, do.points = FALSE, col="black", add=TRUE)
plot(x_lnorm, verticals = TRUE, do.points = FALSE, col="yellow", add=TRUE)

### Chi-Square
# new_df <- df[which(amostra >= 18.5 & df$temp_max <= 35.0),]
# new_fd = hist(new_amostra, col = 'skyblue3', breaks = "fd", main = "Freedman & Diaconis",
#           ylab = "", xlab = "")
# nclass.FD(new_df$temp_max)
# 
# minV <- 20.5; maxV <- 33.5
# new_amostra <- sapply(amostra, function(y) min(max(y,minV),maxV))
# 
# new_fd$breaks
# new_fd$counts
# 
# breaks_normal <- pnorm(new_fd$breaks, mean(df$temp_max),sd(df$temp_max))
# null.probs_normal <- rollapply(breaks_normal,2,function(x) x[2]-x[1])
# chisq.test(new_fd$counts,p=null.probs_normal,rescale.p = TRUE, simulate.p.value = TRUE)
# 
# breaks_gamma <- pgamma(new_fd$breaks, a_gama,b_gama)
# null.probs_gamma <- rollapply(breaks_gamma,2,function(x) x[2]-x[1])
# chisq.test(new_fd$counts,p=null.probs_gamma,rescale.p = TRUE, simulate.p.value = TRUE)
# 
# breaks_gumbel <- pgumbel(new_fd$breaks, a_gumbel,b_gumbel)
# null.probs_gumbel <- rollapply(breaks_gumbel,2,function(x) x[2]-x[1])
# chisq.test(new_fd$counts,p=null.probs_gumbel,rescale.p = TRUE, simulate.p.value = TRUE)
# 
# breaks_lnorm <- plnorm(new_fd$breaks, media,sqrt(variancia))
# null.probs_lnorm <- rollapply(breaks_lnorm,2,function(x) x[2]-x[1])
# res <- chisq.test(new_fd$counts,p=null.probs_lnorm,rescale.p = TRUE, simulate.p.value = TRUE)
# res$observed
# res$expected


# Anderson Darling
# Pacote nortest
ad.test(df$temp_max)
# Pacote desctools
x = df$temp_max
AndersonDarlingTest(x, "pgamma", a_gama, b_gama)
AndersonDarlingTest(x, "pgumbel", a_gumbel, b_gumbel)
AndersonDarlingTest(x, "pnorm", mean(df$temp_max), sd(df$temp_max))
variancia <- log((var(df$temp_max)/exp(2*log(mean(df$temp_max))))+1)
media <- (2*log(mean(df$temp_max))-variancia)/2
AndersonDarlingTest(x, "plnorm", media, sqrt(variancia))

# Lilliefors
lillie.test(df$temp_max)

# reamostragem monte carlo
ks.test(x, "pgamma", a_gama, b_gama, alternative = "greater", exact = NULL, tol=1e-8, simulate.p.value=TRUE, B=2000)
ks.test(x, "pgumbel", a_gumbel, b_gumbel, alternative = "greater", exact = NULL, tol=1e-8, simulate.p.value=TRUE, B=2000)
ks.test(x, "pnorm", mean(df$temp_max), sd(df$temp_max), alternative = "greater", exact = NULL, tol=1e-8, simulate.p.value=TRUE, B=2000)
ks.test(x, "plnorm", media, sqrt(variancia), alternative = "greater", exact = NULL, tol=1e-8, simulate.p.value=TRUE, B=2000)

# reamostragem bootstrap
pop <- df[,2, drop=FALSE]
colnames(pop) <- ("dados")

set.seed(123)

pop %>%
  ggplot(aes(x = dados)) +
  geom_histogram(alpha = 0.3, fill = "green") +
  labs(
    title = "Distribuição Temperatura Máxima Brasilia (1980-2021)",
    subtitle = paste0("Média: ", round(mean(pop$dados)), 
                      ". Desvio-padrão: ", round(sd(pop$dados)))
  ) +
  theme_light()

# Bootstraping (KS) para Gama
# Criar um vetor para armazenar o valor da reamostragem:
reamostras = rep(NA, 1000)
reamostras_p_value = rep(NA, 1000)
# Simular a reamostragem 1.000 vezes:
for (i in 1:1000) {
  reamostra = sample(pop$dados, 500, replace = T) # Reamostragem com reposição
  reamostras[i] = ks.test(reamostra, "pgamma", a_gama, b_gama, alternative = "greater")$statistic
  reamostras_p_value[i] = ks.test(reamostra, "pgamma", a_gama, b_gama, alternative = "greater")$p.value# Calcular o p_value de cada amostra
}
# Transformar o vetor em formato tidy:
reamostras = data.frame(reamostras,reamostras_p_value)
media = mean(reamostras$reamostras)
p_value = reamostras[which.min(abs(media-reamostras$reamostras)),]
quantile(reamostras$reamostras, c(0.025, 0.975))
# Criar o gráfico de distribuição da estatística KS
gama <- reamostras  %>%
  ggplot(aes(x = reamostras)) +
  geom_histogram(alpha = 0.3, fill = "blue", bins = 50) +
  labs(
    title = "Distribuição KS para Gama - Reamostragem (1.000 de n=500)",
    subtitle = paste0("Média: ", round(media,4),
                      ". p-value: ", round(p_value$reamostras_p_value, 4))
  ) + 
  ylim(0,120) +
  theme_light()
gama

# Bootstraping (KS) para Gumbel
# Criar um vetor para armazenar o valor da reamostragem:
reamostras = rep(NA, 1000)
reamostras_p_value = rep(NA, 1000)
# Simular a reamostragem 1.000 vezes:
for (i in 1:1000) {
  reamostra = sample(pop$dados, 500, replace = T) # Reamostragem com reposição
  reamostras[i] = ks.test(reamostra, "pgumbel", a_gumbel, b_gumbel, alternative = "greater")$statistic
  reamostras_p_value[i] = ks.test(reamostra, "pgumbel", a_gumbel, b_gumbel, alternative = "greater")$p.value# Calcular o p_value de cada amostra
}
# Transformar o vetor em formato tidy:
reamostras = data.frame(reamostras,reamostras_p_value)
media = mean(reamostras$reamostras)
p_value = reamostras[which.min(abs(media-reamostras$reamostras)),]
quantile(reamostras$reamostras, c(0.025, 0.975))
# Criar o gráfico de distribuição da estatística KS
gumbel <- reamostras  %>%
  ggplot(aes(x = reamostras)) +
  geom_histogram(alpha = 0.3, fill = "blue", bins = 50) +
  labs(
    title = "Distribuição KS para Gumbel - Reamostragem (1.000 de n=500)",
    subtitle = paste0("Média: ", round(media,4),
                      ". p-value: ", round(p_value$reamostras_p_value, 4))
  ) + 
  ylim(0,120) +
  theme_light()
gumbel

# Bootstraping (KS) para Normal
# Criar um vetor para armazenar o valor da reamostragem:
reamostras = rep(NA, 1000)
reamostras_p_value = rep(NA, 1000)
# Simular a reamostragem 1.000 vezes:
for (i in 1:1000) {
  reamostra = sample(pop$dados, 500, replace = T) # Reamostragem com reposição
  reamostras[i] = ks.test(reamostra, "pnorm", mean(df$temp_max), sd(df$temp_max), alternative = "greater")$statistic
  reamostras_p_value[i] = ks.test(reamostra, "pnorm", mean(df$temp_max), sd(df$temp_max), alternative = "greater")$p.value# Calcular o p_value de cada amostra
}
# Transformar o vetor em formato tidy:
reamostras = data.frame(reamostras,reamostras_p_value)
media = mean(reamostras$reamostras)
p_value = reamostras[which.min(abs(media-reamostras$reamostras)),]
quantile(reamostras$reamostras, c(0.025, 0.975))
# Criar o gráfico de distribuição da estatística KS
normal <- reamostras  %>%
  ggplot(aes(x = reamostras)) +
  geom_histogram(alpha = 0.3, fill = "blue", bins = 50) +
  labs(
    title = "Distribuição KS para Normal - Reamostragem (1.000 de n=500)",
    subtitle = paste0("Média: ", round(media,4),
                      ". p-value: ", round(p_value$reamostras_p_value, 4))
  ) + 
  ylim(0,120) +
  theme_light()
normal

# Bootstraping (KS) para Log-normal
# Criar um vetor para armazenar o valor da reamostragem:
reamostras = rep(NA, 1000)
reamostras_p_value = rep(NA, 1000)
variancia <- log((var(df$temp_max)/exp(2*log(mean(df$temp_max))))+1)
media <- (2*log(mean(df$temp_max))-variancia)/2
# Simular a reamostragem 1.000 vezes:
for (i in 1:1000) {
  reamostra = sample(pop$dados, 500, replace = T) # Reamostragem com reposição
  reamostras[i] = ks.test(reamostra, "plnorm", media, sqrt(variancia), alternative = "greater")$statistic
  reamostras_p_value[i] = ks.test(reamostra, "plnorm", media, sqrt(variancia), alternative = "greater")$p.value# Calcular o p_value de cada amostra
}
# Transformar o vetor em formato tidy:
reamostras = data.frame(reamostras,reamostras_p_value)
media = mean(reamostras$reamostras)
p_value = reamostras[which.min(abs(media-reamostras$reamostras)),]
quantile(reamostras$reamostras, c(0.025, 0.975))
# Criar o gráfico de distribuição da estatística KS
log_normal <- reamostras  %>%
  ggplot(aes(x = reamostras)) +
  geom_histogram(alpha = 0.3, fill = "blue", bins = 50) +
  labs(
    title = "Distribuição KS para Log-normal - Reamostragem (1.000 de n=500)",
    subtitle = paste0("Média: ", round(media,4),
                      ". p-value: ", round(p_value$reamostras_p_value, 4))
  ) + 
  ylim(0,120) +
  theme_light()
log_normal

ggarrange(log_normal + rremove("xylab"), gama + rremove("xylab"), normal + rremove("xylab"), gumbel + rremove("xylab"), ncol = 2, nrow = 2)

# Bootstraping (Lilliefors) para Normal
# Criar um vetor para armazenar o valor da reamostragem:
reamostras = rep(NA, 1000)
reamostras_p_value = rep(NA, 1000)
# Simular a reamostragem 1.000 vezes:
for (i in 1:1000) {
  reamostra = sample(pop$dados, 500, replace = T) # Reamostragem com reposição
  reamostras[i] = lillie.test(reamostra)$statistic
  reamostras_p_value[i] = lillie.test(reamostra)$p.value# Calcular o p_value de cada amostra
}
# Transformar o vetor em formato tidy:
reamostras = data.frame(reamostras,reamostras_p_value)
media = mean(reamostras$reamostras)
p_value = reamostras[which.min(abs(media-reamostras$reamostras)),]
quantile(reamostras$reamostras, c(0.025, 0.975))
# Criar o gráfico de distribuição da estatística Lilliefors
normal_lilliefors <- reamostras  %>%
  ggplot(aes(x = reamostras)) +
  geom_histogram(alpha = 0.3, fill = "blue", bins = 50) +
  labs(
    title = "Distribuição Lilliefors para Normal - Reamostragem (1.000 de n=500)",
    subtitle = paste0("Média: ", round(media,4),
                      ". p-value: ", round(p_value$reamostras_p_value, 4))
  ) + 
  ylim(0,120) +
  theme_light()
normal_lilliefors


# Bootstraping (AD) para Gama
# Criar um vetor para armazenar o valor da reamostragem:
reamostras = rep(NA, 1000)
reamostras_p_value = rep(NA, 1000)
# Simular a reamostragem 1.000 vezes:
for (i in 1:1000) {
  reamostra = sample(pop$dados, 500, replace = T) # Reamostragem com reposição
  reamostras[i] = AndersonDarlingTest(reamostra, "pgamma", a_gama, b_gama)$statistic
  reamostras_p_value[i] = AndersonDarlingTest(reamostra, "pgamma", a_gama, b_gama)$p.value# Calcular o p_value de cada amostra
}
# Transformar o vetor em formato tidy:
reamostras = data.frame(reamostras,reamostras_p_value)
media = mean(reamostras$reamostras)
p_value = reamostras[which.min(abs(media-reamostras$reamostras)),]
quantile(reamostras$reamostras, c(0.025, 0.975))
# Criar o gráfico de distribuição da estatística AD
gama <- reamostras  %>%
  ggplot(aes(x = reamostras)) +
  geom_histogram(alpha = 0.3, fill = "blue", bins = 50) +
  labs(
    title = "Distribuição AD para Gama - Reamostragem (1.000 de n=500)",
    subtitle = paste0("Média: ", round(media,4),
                      ". p-value: ", round(p_value$reamostras_p_value, 4))
  ) + 
  ylim(0,120) +
  theme_light()
gama

# Bootstraping (AD) para Gumbel
# Criar um vetor para armazenar o valor da reamostragem:
reamostras = rep(NA, 1000)
reamostras_p_value = rep(NA, 1000)
# Simular a reamostragem 1.000 vezes:
for (i in 1:1000) {
  reamostra = sample(pop$dados, 500, replace = T) # Reamostragem com reposição
  reamostras[i] = AndersonDarlingTest(reamostra, "pgumbel", a_gumbel, b_gumbel)$statistic
  reamostras_p_value[i] = AndersonDarlingTest(reamostra, "pgumbel", a_gumbel, b_gumbel)$p.value# Calcular o p_value de cada amostra
}
# Transformar o vetor em formato tidy:
reamostras = data.frame(reamostras,reamostras_p_value)
media = mean(reamostras$reamostras)
p_value = reamostras[which.min(abs(media-reamostras$reamostras)),]
quantile(reamostras$reamostras, c(0.025, 0.975))
# Criar o gráfico de distribuição da estatística AD
gumbel <- reamostras  %>%
  ggplot(aes(x = reamostras)) +
  geom_histogram(alpha = 0.3, fill = "blue", bins = 50) +
  labs(
    title = "Distribuição AD para Gumbel - Reamostragem (1.000 de n=500)",
    subtitle = paste0("Média: ", round(media,4),
                      ". p-value: ", round(p_value$reamostras_p_value, 4))
  ) + 
  ylim(0,120) +
  theme_light()
gumbel

# Bootstraping (AD) para Normal
# Criar um vetor para armazenar o valor da reamostragem:
reamostras = rep(NA, 1000)
reamostras_p_value = rep(NA, 1000)
# Simular a reamostragem 1.000 vezes:
for (i in 1:1000) {
  reamostra = sample(pop$dados, 500, replace = T) # Reamostragem com reposição
  reamostras[i] = AndersonDarlingTest(reamostra, "pnorm", mean(df$temp_max), sd(df$temp_max))$statistic
  reamostras_p_value[i] = AndersonDarlingTest(reamostra, "pnorm", mean(df$temp_max), sd(df$temp_max))$p.value# Calcular o p_value de cada amostra
}
# Transformar o vetor em formato tidy:
reamostras = data.frame(reamostras,reamostras_p_value)
media = mean(reamostras$reamostras)
p_value = reamostras[which.min(abs(media-reamostras$reamostras)),]
quantile(reamostras$reamostras, c(0.025, 0.975))
# Criar o gráfico de distribuição da estatística AD
normal <- reamostras  %>%
  ggplot(aes(x = reamostras)) +
  geom_histogram(alpha = 0.3, fill = "blue", bins = 50) +
  labs(
    title = "Distribuição AD para Normal - Reamostragem (1.000 de n=500)",
    subtitle = paste0("Média: ", round(media,4),
                      ". p-value: ", round(p_value$reamostras_p_value, 4))
  ) + 
  ylim(0,120) +
  theme_light()
normal

# Bootstraping (AD) para Log-normal
# Criar um vetor para armazenar o valor da reamostragem:
reamostras = rep(NA, 1000)
reamostras_p_value = rep(NA, 1000)
variancia <- log((var(df$temp_max)/exp(2*log(mean(df$temp_max))))+1)
media <- (2*log(mean(df$temp_max))-variancia)/2
# Simular a reamostragem 1.000 vezes:
for (i in 1:1000) {
  reamostra = sample(pop$dados, 500, replace = T) # Reamostragem com reposição
  reamostras[i] = AndersonDarlingTest(reamostra, "plnorm", media, sqrt(variancia))$statistic
  reamostras_p_value[i] = AndersonDarlingTest(reamostra, "plnorm", media, sqrt(variancia))$p.value# Calcular o p_value de cada amostra
}
# Transformar o vetor em formato tidy:
reamostras = data.frame(reamostras,reamostras_p_value)
media = mean(reamostras$reamostras)
p_value = reamostras[which.min(abs(media-reamostras$reamostras)),]
quantile(reamostras$reamostras, c(0.025, 0.975))
# Criar o gráfico de distribuição da estatística AD
log_normal <- reamostras  %>%
  ggplot(aes(x = reamostras)) +
  geom_histogram(alpha = 0.3, fill = "blue", bins = 50) +
  labs(
    title = "Distribuição AD para Log-normal - Reamostragem (1.000 de n=500)",
    subtitle = paste0("Média: ", round(media,4),
                      ". p-value: ", round(p_value$reamostras_p_value, 4))
  ) + 
  ylim(0,120) +
  theme_light()
log_normal

ggarrange(normal + rremove("xylab"), gama + rremove("xylab"), log_normal + rremove("xylab"), gumbel + rremove("xylab"), ncol = 2, nrow = 2)
nclass.FD(reamostra)

gofTest(reamostra, test = "chisq", distribution = "gamma", n.classes = 25)

#Chi-squared
gofTest(df$temp_max, test = "chisq", n.classes = 82)
gofTest(df$temp_max, test = "chisq", distribution = "gamma", n.classes = 82)
gofTest(df$temp_max, test = "chisq", distribution = "lnorm", n.classes = 82)

# Bootstraping (Chi-squared) para Gama
# Criar um vetor para armazenar o valor da reamostragem:
reamostras = rep(NA, 1000)
reamostras_p_value = rep(NA, 1000)
# Simular a reamostragem 1.000 vezes:
for (i in 1:1000) {
  reamostra = sample(pop$dados, 500, replace = T) # Reamostragem com reposição
  reamostras[i] = gofTest(reamostra, test = "chisq", distribution = "gamma", n.classes = 10)$statistic
  reamostras_p_value[i] = gofTest(reamostra, test = "chisq", distribution = "gamma", n.classes = 25)$p.value# Calcular o p_value de cada amostra
}
# Transformar o vetor em formato tidy:
reamostras = data.frame(reamostras,reamostras_p_value)
media = mean(reamostras$reamostras)
p_value = reamostras[which.min(abs(media-reamostras$reamostras)),]
quantile(reamostras$reamostras, c(0.025, 0.975))
# Criar o gráfico de distribuição da estatística Chi-squared
gama <- reamostras  %>%
  ggplot(aes(x = reamostras)) +
  geom_histogram(alpha = 0.3, fill = "blue", bins = 50) +
  labs(
    title = "Distribuição Chi-squared para Gama - Reamostragem (1.000 de n=500)",
    subtitle = paste0("Média: ", round(media,4),
                      ". p-value: ", round(p_value$reamostras_p_value, 4))
  ) + 
  ylim(0,120) +
  theme_light()
gama

# Bootstraping (Chi-squared) para Normal
# Criar um vetor para armazenar o valor da reamostragem:
reamostras = rep(NA, 1000)
reamostras_p_value = rep(NA, 1000)
# Simular a reamostragem 1.000 vezes:
for (i in 1:1000) {
  reamostra = sample(pop$dados, 500, replace = T) # Reamostragem com reposição
  reamostras[i] = gofTest(reamostra, test = "chisq", n.classes = 10)$statistic
  reamostras_p_value[i] = gofTest(reamostra, test = "chisq", n.classes = 25)$p.value# Calcular o p_value de cada amostra
}
# Transformar o vetor em formato tidy:
reamostras = data.frame(reamostras,reamostras_p_value)
media = mean(reamostras$reamostras)
p_value = reamostras[which.min(abs(media-reamostras$reamostras)),]
quantile(reamostras$reamostras, c(0.025, 0.975))
# Criar o gráfico de distribuição da estatística Chi-squared
normal <- reamostras  %>%
  ggplot(aes(x = reamostras)) +
  geom_histogram(alpha = 0.3, fill = "blue", bins = 50) +
  labs(
    title = "Distribuição Chi-squared para Normal - Reamostragem (1.000 de n=500)",
    subtitle = paste0("Média: ", round(media,4),
                      ". p-value: ", round(p_value$reamostras_p_value, 4))
  ) + 
  ylim(0,120) +
  theme_light()
normal

# Bootstraping (Chi-squared) para Log-normal
# Criar um vetor para armazenar o valor da reamostragem:
reamostras = rep(NA, 1000)
reamostras_p_value = rep(NA, 1000)
variancia <- log((var(df$temp_max)/exp(2*log(mean(df$temp_max))))+1)
media <- (2*log(mean(df$temp_max))-variancia)/2
# Simular a reamostragem 1.000 vezes:
for (i in 1:1000) {
  reamostra = sample(pop$dados, 500, replace = T) # Reamostragem com reposição
  reamostras[i] = gofTest(reamostra, test = "chisq", distribution = "lnorm", n.classes = 10)$statistic
  reamostras_p_value[i] =gofTest(reamostra, test = "chisq", distribution = "lnorm", n.classes = 10)$p.value# Calcular o p_value de cada amostra
}
# Transformar o vetor em formato tidy:
reamostras = data.frame(reamostras,reamostras_p_value)
media = mean(reamostras$reamostras)
p_value = reamostras[which.min(abs(media-reamostras$reamostras)),]
quantile(reamostras$reamostras, c(0.025, 0.975))
# Criar o gráfico de distribuição da estatística AD
log_normal <- reamostras  %>%
  ggplot(aes(x = reamostras)) +
  geom_histogram(alpha = 0.3, fill = "blue", bins = 50) +
  labs(
    title = "Distribuição AD para Log-normal - Reamostragem (1.000 de n=500)",
    subtitle = paste0("Média: ", round(media,4),
                      ". p-value: ", round(p_value$reamostras_p_value, 4))
  ) + 
  ylim(0,120) +
  theme_light()
log_normal

ggarrange(gama + rremove("xylab"), log_normal + rremove("xylab"), normal + rremove("xylab"), ncol = 2, nrow = 2)
nclass.FD(reamostra)

# Log-normal Distribution
variancia <- log((var(df$temp_max)/exp(2*log(mean(df$temp_max))))+1)
media <- (2*log(mean(df$temp_max))-variancia)/2
variancia
media

set.seed(42)
lnormal_dist <- rlnorm(n=length(df$temp_max), media, sqrt(variancia))
lnormal_dist

x <- lnormal_dist

hist(lnormal_dist, breaks = "fd", freq = FALSE, main = "Distribuição Log-normal", ylab = "Densidade", xlab = "Temperatura (°C)")
curve(dlnorm(x, media, sqrt(variancia)), col = "black", lwd = 2, add=TRUE)

# Função geratriz de momento
x <- lnormal_dist

variancia_ln <- log((var(x)/exp(2*log(mean(x))))+1)
media_ln <- (2*log(mean(x))-variancia_ln)/2

# Momento 1:
# From library "moments"
momento1 <- moment(x, order=1)
momento1

moment(x, order=2)

# Momento 2:
# From library "moments"
momento2 <- moment(x, central=TRUE, order=2)
momento2

# Momento 3:
# From library "moments"
momento3 <- moment(x, central=TRUE, order=3)
momento3

# Momento 4:
# From library "moments"
momento4 <- moment(x, central=TRUE, order=4)
momento4

# Gráfico de Densidade
plot(density(x), xlab='Temperatura Máxima', ylab='Densidade', main='Gráfico de Densidade da Distribuição Log-normal')
abline(v=mean(x), col='red') 

4.209049 / 6.081948 ^ (3/2) #Assimetria da Log-normal
116.8326 / 6.081948 ^ 2 # Curtose da Log-normal

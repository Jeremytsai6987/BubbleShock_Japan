getwd()
library(data.table)
library(readxl)
Path = "/Users/caiyawei/Desktop/VAR"
setwd(Path)
source("/Users/caiyawei/Desktop/VAR/VAR_functions.r")  
inv_tol = 1e-20
options(warn=-1)  
options(scipen=999) 
library(texreg)

### DATA ###
file = "/Users/caiyawei/Desktop/Japan Data/SVAR_forR.xlsx"
data = read_excel(file)
data = data[-c(107),]

###Order1 Stock -> House -> Macro###
### Employment ###
#data$date <- as.Date(data$date)

By <- data %>% 
  select(stock_price, property_price, employment) %>% 
  as.matrix

### RAW DATA ###
dim(By)
kk <- dim(By)[2]
draw_raw <- function(vlb_name = ""){
  vlb_name <- sym(vlb_name)
  plot.return <- data %>%
    ggplot()+
    geom_line(aes(year, !!vlb_name))
  return(plot.return)
}

p1 <- draw_raw("stock_price")
p2 <- draw_raw("property_price")
p3 <- draw_raw("employment")

multiplot(p1, p2, p3)

### Setting Model ###
VAR.P = 4                       
CONST = TRUE                    
Y     = VAR.Y(By, VAR.P)        
X     = VAR.X(By, VAR.P)        

hrz=39 # the length of response

### Reduced Form VAR ###
(Coef.OLS    = VAR.OLS(Y, X, CONST)                  )
(Sigma.OLS   = VAR.Sigma.OLS(Y, X, Coef.OLS, CONST)  )
(Sigma.MLE   = VAR.Sigma.MLE(Y, X, Coef.OLS, CONST))

VAR.P = 4 # By AIC

### SVAR ###
### 6-variable model

Amat = diag(3)
# Identification Conditions

Amat[2,1]  = NA; 
Amat[3,1]  = NA; Amat[3,2]  = NA; 


Bmat = diag(3)
diag(Bmat) = NA

Amat;Bmat

### Estimate A^hat and B^hat ###
C.Prime <- chol(Sigma.OLS)
C <- t(C.Prime)
C

### Solving system of linear equations ###
B0 <- diag(diag(C), ncol = 3, nrow = 3)
B0

A0 <- B0 %*% solve(C)
A0

SVAR_AB_est <- list("A0.svar" = A0, "B0.svar" = B0)

### IRF Without Bootstrap CI ###
### IRF
SVAR_AB_IRF <- VAR.svarirf.AB(By, VAR.P, Amat, Bmat, h = hrz, CONST, SVAR_AB_est = SVAR_AB_est)


# 6*6 time series
df_IRF_plot <- matrix(NA, hrz+1, kk^2) #%>% as.tibble() ## hrz+1
#dim(df_IRF_plot)
h <- 0 
for(period in SVAR_AB_IRF){
  k <- 0 
  h <- h+1 
  for(j in 1:kk){
    for(i in 1:kk){
      k <- k+1 
      df_IRF_plot[h,k] <- period[i,j]
    }
  }
}
df_IRF_plot <- df_IRF_plot %>% as_tibble()

kk*1:kk

#output entire table
IRF_TABLE <- df_IRF_plot[,kk*1:kk] %>% select(stock_price = 1,
                                              property_price = 2,
                                              employment = 3)
write.table(IRF_TABLE, file = "/Users/caiyawei/Desktop/IRF_TABLE.csv", sep = ",", row.names = FALSE)

draw_IRF <- function(df = df_IRF_plot, V1 = 1){
  V1 <- paste0("V", V1) %>% sym()
  plot.return <- ggplot(df) + 
    geom_line(aes(x = 1:nrow(df), y = !!V1))
  return(plot.return)
}

for(i in 1:(kk^2)){
  assign(paste0("p",i), draw_IRF(df_IRF_plot, i))
}

multiplot(p1,p2,p3,p4,p5,p6,
          p7,p8,p9,cols = 3)

### IRF (Bootstrap C.I.) ###
lower = 0.025                                        # 95% CI
upper = 1-lower
kk = ncol(By)
ddY = VAR.ddY(By, VAR.P)
ddX = VAR.ddX(By, VAR.P)

# dim(ddY); dim(ddX)

T   = nrow(ddY)
T.total= nrow(By)
Ik  = diag(rep(1, kk))
# 16 coef if 4 variables; 55 coef if 5 variables
Coef = t(VAR.EbyE(ddY, ddX, CONST)$ddA)              
# residuals
U    = VAR.EbyE(ddY, ddX, CONST)$ddU
BSigma.u = VAR.ddSigma.OLS(ddY, ddX, CONST)
if(CONST == TRUE){
  const = Coef[, ncol(Coef)]
  Coef.noc= Coef[,-ncol(Coef)]                      
}else{
  const = matrix(0, kk, 1)
  Coef.noc = Coef
}

Theta.unit= VAR.Theta(Coef, h, BSigma.u, CONST)$unit 
Theta.std = VAR.Theta(Coef, h, BSigma.u, CONST)$std  

# dm.U <- U-mean(U)
dm.U <- U

N = 2000 
Theta.unit.sim = vector("list", N)
Theta.std.sim  = vector("list", N)

# check dimension
print("check dimensionality")
dim(ddX); dim(Coef.noc); dim(dm.U)


df_IRF.sim <- array(NA, c(hrz+1,kk^2,N)) #dimensions are: Time Period, Number of shock interacts with variables, page (number of Bootstrap resamplings)
counter <- 1
while(TRUE){
  
  #cat("Now, there are ", counter-1, " sets of resamples.\n")
  Y.sim = matrix(0, nrow = T.total, ncol = kk)          # Y.sim = 0 #pseudo time series
  Y.sim[c(1:VAR.P),] = By[c(1:VAR.P), ] #initial values
  
  boot.number = sample(c(1:T), replace = TRUE)      
  U.sim = dm.U[boot.number,]
  
  # predicted values given the above initial values
  last.y= c(t(By[VAR.P:1,]))
  for(ii in 1:T){
    last.y = last.y[1:(kk*VAR.P)]
    Y.sim[ii+VAR.P, ] = Coef.noc %*% last.y + const + U.sim[ii,]     
    last.y = c(Y.sim[ii+VAR.P,], last.y)
  }
  
  #   Y.sim[-c(1:VAR.P),] <- matrix(const, nrow = T.total-VAR.P, ncol = kk, byrow = T) + ddX %*% t(Coef.noc) + U.sim
  
  
  #`Y.sim` is the pseudo time series
  
  ### SVAR.sim Start ###
  
  Y_pseudo     = VAR.Y(Y.sim, VAR.P)        
  X_pseudo     = VAR.X(Y.sim, VAR.P)        
  Coef.OLS_pseudo    = VAR.OLS(Y_pseudo, X_pseudo, CONST)
  Sigma.OLS_pseudo   = VAR.Sigma.OLS(Y_pseudo, X_pseudo, Coef.OLS_pseudo, CONST)
  C.Prime_pseudo <- chol(Sigma.OLS_pseudo)
  C_pseudo <- t(C.Prime_pseudo)
  B0_pseudo <- diag(diag(C_pseudo), ncol = kk, nrow = kk)
  A0_pseudo <- B0_pseudo %*% solve(C_pseudo)
  SVAR_AB_est.sim <- list("A0.svar" = A0_pseudo, "B0.svar" = B0_pseudo)
  SVAR_AB_IRF.sim <- VAR.svarirf.AB(Y.sim, VAR.P, Amat, Bmat, h = hrz, CONST, SVAR_AB_est = SVAR_AB_est.sim)
  
  # 5*5 time series
  df_IRF_plot.sim <- matrix(NA, hrz+1, kk^2) #%>% as.tibble()
  # df_IRF.sim <- array(1:(120*25*N), c(120,25,N))
  # df_IRF.sim[2,1,1] # slicing
  
  h <- 0 # h period of IRF
  for(period in SVAR_AB_IRF.sim){
    k <- 0 
    h <- h+1 # h = period of IRF
    for(j in 1:kk){
      for(i in 1:kk){
        k <- k+1 
        df_IRF_plot.sim[h,k] <- period[i,j]
      }
    }
  }
  df_IRF.sim[,,counter] <- df_IRF_plot.sim
  ### SVAR.sim Ends ###
  if(counter>=N){
    break
  }
  counter <- counter+1
}
# Save
saveRDS(df_IRF.sim, file = "df_IRF.sim.rds")

df_IRF.sim <- read_rds("df_IRF.sim.rds")

head(df_IRF.sim[,,1000])
print(sum(is.na(df_IRF.sim)))

# Draw IRF & Bootstrap C.I.
df_IRF_plot.BS.L <- matrix(NA, nrow = hrz+1, ncol = kk^2)
df_IRF_plot.BS.U <- matrix(NA, nrow = hrz+1, ncol = kk^2)
df_IRF_plot.BS.Median <- matrix(NA, nrow = hrz+1, ncol = kk^2)
df_IRF_plot.BS.Mean <- matrix(NA, nrow = hrz+1, ncol = kk^2)
for(col in 1:(kk^2)){
  for(row in 1:(hrz+1) ){
    df_IRF_plot.BS.L[row,col] <- quantile(df_IRF.sim[row,col,], probs = 0.025)
    df_IRF_plot.BS.U[row,col] <- quantile(df_IRF.sim[row,col,], probs = 0.975)
    df_IRF_plot.BS.Median[row,col] <- quantile(df_IRF.sim[row,col,], probs = 0.5)
    df_IRF_plot.BS.Mean[row,col] <- mean(df_IRF.sim[row,col,])
  }
}

df_IRF_plot.BS.L <- df_IRF_plot.BS.L %>% as_tibble()
df_IRF_plot.BS.U <- df_IRF_plot.BS.U %>% as_tibble()
df_IRF_plot.BS.Median <- df_IRF_plot.BS.Median %>% as_tibble()
df_IRF_plot.BS.Mean <- df_IRF_plot.BS.Mean %>% as_tibble()

ind <- 0
for(i in 1:kk){
  for(j in 1:kk){
    ind <- ind+1
    nam <- paste("shock", j, "y", i, sep = '')
    assign(nam, bind_cols(df_IRF_plot.BS.L[ind], df_IRF_plot.BS.U[ind],
                          df_IRF_plot.BS.Median[ind], df_IRF_plot.BS.Mean[ind],
                          df_IRF_plot[ind]))
    evalStr <- paste0("colnames(", nam, ") <- c('Lower', 'Upper', 'Median', 'Mean', 'Actual')")
    eval(parse(text=evalStr))
    evalStr <- paste0("p", ind, " <- ", "ggplot(",nam,") +geom_hline(yintercept=0, color = 'grey')+ geom_line(aes(x = 1:nrow(", nam, "), y = Lower), linetype = 'dashed', col='red')+geom_line(aes(x = 1:nrow(", nam, "), y = Upper), linetype = 'dashed', col='red')+geom_line(aes(x = 1:nrow(", nam, "), y = Median), col = 'Blue')")
    eval(parse(text=evalStr))
  }
}  

Text_Size_Theme = theme(
  axis.title.x = element_text(size = 12),
  axis.text.x = element_text(size = 12),
  axis.title.y = element_text(size = 12),
  axis.title = element_text(size = 12),
  plot.title = element_text(size=12))

add_layout <- function(p = p1,
                       response_of = "this variable", react_to = "the shock"){
  title_text <- paste("Response of", response_of, "to", react_to, sep = ' ')
  plot.return <- p + labs(x = "Period", y = "", title = title_text) + Text_Size_Theme
  return(plot.return)
}

vlb_list = c("Stock Price",
             "House Price",
             "Employment"
)
p1 <- add_layout(p1, vlb_list[1], "Stock shock")
p2 <- add_layout(p2, vlb_list[2], "Stock shock")
p3 <- add_layout(p3, vlb_list[3], "Stock shock")

p4 <- add_layout(p4, vlb_list[1], "Property Price shock")
p5 <- add_layout(p5, vlb_list[2], "Property Price shock")
p6 <- add_layout(p6, vlb_list[3], "Property Price shock")

p7 <- add_layout(p7, vlb_list[1], "Employment shock")
p8 <- add_layout(p8, vlb_list[2], "Employment shock")
p9 <- add_layout(p9, vlb_list[3], "Employment shock")

multiplot(p1,p2,p3,p4,p5,p6,
          p7,p8,p9,cols = 3)

multiplot(p1,p2,p3, cols = 1)
ggsave(filename = "/Users/caiyawei/Desktop/VAR/IRF_shock1.png", 
       plot = multiplot(p1,p2,p3, cols = 1),
       width = 30, height = 20, units = "cm",
       device = "png")

multiplot(p4,p5,p6, cols = 1)
ggsave(filename = "/Users/caiyawei/Desktop/VAR/IRF_shock2.png", 
       plot = multiplot(p4,p5,p6, cols = 1),
       width = 30, height = 20, units = "cm",
       device = "png")

multiplot(p7,p8,p9, cols = 1)
ggsave(filename = "/Users/caiyawei/Desktop/VAR/IRF_shock3.png", 
       plot = multiplot(p7,p8,p9, cols = 1),
       width = 30, height = 20, units = "cm",
       device = "png")

### Obtain the exogenous shock###
Coef.noc
Coef.OLS
const

# check the order of variables is correct (stock shock) (employment)
lm(stock_price ~ lag(stock_price)+lag(property_price)+lag(employment)+
     +lag(stock_price, n=2)+lag(property_price, n=2)+lag(employment, n=2)+lag(stock_price, n=3)+lag(property_price, n=3)+lag(employment, n=3)+lag(stock_price, n=4)+lag(property_price, n=4)+lag(employment, n=4),
   data = data)

get_residual <- function(vlb_order = 1){
  X.Mat <- By %>% as_tibble() %>%
    mutate(V1 = lag(stock_price),
           V2 = lag(property_price),
           V3 = lag(employment),
           V4 = lag(stock_price, n = 2),
           V5 = lag(property_price, n = 2),
           V6 = lag(employment, n = 2),
           V7 = lag(stock_price, n = 3),
           V8 = lag(property_price, n = 3),
           V9 = lag(employment, n = 3),
           V10 = lag(stock_price, n = 4),
           V11 = lag(property_price, n = 4),
           V12 = lag(employment, n = 4),
           V13 = rep(1, nrow(By))
    ) %>%
    select(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13
    ) %>%
    drop_na() %>%
    as.matrix
  Y.Mat <- By[(VAR.P+1):nrow(By),vlb_order] %>% as.matrix
  Y.Mat_hat <- X.Mat %*% (Coef.OLS[vlb_order, ] %>% as.matrix)
  residual <- Y.Mat - Y.Mat_hat
  return(residual)
}

# resid1 <- get_residual(1)
for(i in 1:3){
  assign(paste0("resid_",i), get_residual(i))
}

resid.Mat <- cbind(resid_1, resid_2, resid_3)
error.Mat <- (solve(B0) %*% A0 %*% t(resid.Mat)) %>% t()
error.Mat[,1] %>% hist

# Extract the stock shock

error.Mat %>% as_tibble() %>%
  select(e1 = 1, e2 = 2, e3 = 3) %>%
  mutate(no = row_number()) %>%
  ggplot()+
  geom_line(aes(no, e1))

data.err <- data %>% 
  bind_cols(e1 = c(rep(NA, VAR.P), error.Mat[,1])) %>%
  mutate(e1_pos = fifelse(e1 > 0, e1, 0),
         e1_neg = fifelse(e1 < 0, e1, 0))
data.err %>% head(10)

data.err %>% ggplot()+
  geom_line(aes(year, e1))

# check the order of variables is correct (house shock) (employment)
lm(property_price ~ lag(stock_price)+lag(property_price)+lag(employment)+
     +lag(stock_price, n=2)+lag(property_price, n=2)+lag(employment, n=2)+lag(stock_price, n=3)+lag(property_price, n=3)+lag(employment, n=3)+lag(stock_price, n=4)+lag(property_price, n=4)+lag(employment, n=4),
   data = data)

get_residual <- function(vlb_order = 1){
  X.Mat <- By %>% as_tibble() %>%
    mutate(V1 = lag(stock_price),
           V2 = lag(property_price),
           V3 = lag(employment),
           V4 = lag(stock_price, n = 2),
           V5 = lag(property_price, n = 2),
           V6 = lag(employment, n = 2),
           V7 = lag(stock_price, n = 3),
           V8 = lag(property_price, n = 3),
           V9 = lag(employment, n = 3),
           V10 = lag(stock_price, n = 4),
           V11 = lag(property_price, n = 4),
           V12 = lag(employment, n = 4),
           V13 = rep(1, nrow(By))
    ) %>%
    select(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13
    ) %>%
    drop_na() %>%
    as.matrix
  Y.Mat <- By[(VAR.P+1):nrow(By),vlb_order] %>% as.matrix
  Y.Mat_hat <- X.Mat %*% (Coef.OLS[vlb_order, ] %>% as.matrix)
  residual <- Y.Mat - Y.Mat_hat
  return(residual)
}

# resid1 <- get_residual(1)
for(i in 1:3){
  assign(paste0("resid_",i), get_residual(i))
}

resid.Mat <- cbind(resid_1, resid_2, resid_3)
error.Mat <- (solve(B0) %*% A0 %*% t(resid.Mat)) %>% t()
error.Mat[,2] %>% hist

# Extract the house shock

error.Mat %>% as_tibble() %>%
  select(e1 = 1, e2 = 2, e3 = 3) %>%
  mutate(no = row_number()) %>%
  ggplot()+
  geom_line(aes(no, e2))

data.err2 <- data %>% 
  bind_cols(e2 = c(rep(NA, VAR.P), error.Mat[,2])) %>%
  mutate(e2_pos = fifelse(e2 > 0, e2, 0),
         e2_neg = fifelse(e2 < 0, e2, 0))
data.err2 %>% head(10)

data.err2 %>% ggplot()+
  geom_line(aes(year, e2))

#data.err["Period"] = data$year - 7
data.err["e2"] = data.err2$e2
### Regression ###
#reg1 <- lm(em_d ~ Period + e1 + e2 + interest_rate_d, data = data.err)
#summary(reg1)

###衝擊也減去1991Q4###
data.err["e2"] = data.err["e2"] - data.err$e2[6]
data.err["e1"] = data.err["e1"] - data.err$e1[6]
#reg1 <- lm(em_d ~ Period + e1 + e2 + interest_rate_d, data = data.err)
reg1 <- lm(em_d ~ e1 + e2 + interest_rate_d, data = data.err)
summary(reg1)


### DATA CPI ###
file = "/Users/caiyawei/Desktop/Japan Data/SVAR_forR.xlsx"
data_cpi = read_excel(file)
data_cpi = data_cpi[-c(107),]

By <- data_cpi %>% 
  select(stock_price, property_price, cpi) %>% 
  as.matrix

### RAW DATA ###
dim(By)
kk <- dim(By)[2]
draw_raw <- function(vlb_name = ""){
  vlb_name <- sym(vlb_name)
  plot.return <- data_cpi %>%
    ggplot()+
    geom_line(aes(year, !!vlb_name))
  return(plot.return)
}

p1 <- draw_raw("stock_price")
p2 <- draw_raw("property_price")
p3 <- draw_raw("cpi")

multiplot(p1, p2, p3)

### Setting Model ###
VAR.P = 4                       
CONST = TRUE                    
Y     = VAR.Y(By, VAR.P)      
X     = VAR.X(By, VAR.P)        

hrz=39 # the length of response

### Reduced Form VAR ###
(Coef.OLS    = VAR.OLS(Y, X, CONST)                  )
(Sigma.OLS   = VAR.Sigma.OLS(Y, X, Coef.OLS, CONST)  )
(Sigma.MLE   = VAR.Sigma.MLE(Y, X, Coef.OLS, CONST))

VAR.P = 4 # By AIC

### SVAR ###
### 6-variable model

Amat = diag(3)
# Identification Conditions

Amat[2,1]  = NA; 
Amat[3,1]  = NA; Amat[3,2]  = NA; 


Bmat = diag(3)
diag(Bmat) = NA

Amat;Bmat

### Estimate A^hat and B^hat ###
C.Prime <- chol(Sigma.OLS)
C <- t(C.Prime)
C

### Solving system of linear equations ###
B0 <- diag(diag(C), ncol = 3, nrow = 3)
B0

A0 <- B0 %*% solve(C)
A0

SVAR_AB_est <- list("A0.svar" = A0, "B0.svar" = B0)

### IRF Without Bootstrap CI ###
### IRF
SVAR_AB_IRF <- VAR.svarirf.AB(By, VAR.P, Amat, Bmat, h = hrz, CONST, SVAR_AB_est = SVAR_AB_est)


df_IRF_plot <- matrix(NA, hrz+1, kk^2) #%>% as.tibble() ## hrz+1
#dim(df_IRF_plot)
h <- 0 
for(period in SVAR_AB_IRF){
  k <- 0 
  h <- h+1 
  for(j in 1:kk){
    for(i in 1:kk){
      k <- k+1 
      df_IRF_plot[h,k] <- period[i,j]
    }
  }
}
df_IRF_plot <- df_IRF_plot %>% as_tibble()

kk*1:kk

#output entire table
IRF_TABLE <- df_IRF_plot[,kk*1:kk] %>% select(stock_price = 1,
                                              property_price = 2,
                                              cpi = 3)
write.table(IRF_TABLE, file = "/Users/caiyawei/Desktop/IRF_TABLE3.csv", sep = ",", row.names = FALSE)

draw_IRF <- function(df = df_IRF_plot, V1 = 1){
  V1 <- paste0("V", V1) %>% sym()
  plot.return <- ggplot(df) + 
    geom_line(aes(x = 1:nrow(df), y = !!V1))
  return(plot.return)
}

for(i in 1:(kk^2)){
  assign(paste0("p",i), draw_IRF(df_IRF_plot, i))
}

multiplot(p1,p2,p3,p4,p5,p6,
          p7,p8,p9,cols = 3)

### IRF (Bootstrap C.I.) ###
lower = 0.025                                        # 95% CI
upper = 1-lower
kk = ncol(By)
ddY = VAR.ddY(By, VAR.P)
ddX = VAR.ddX(By, VAR.P)

# dim(ddY); dim(ddX)

T   = nrow(ddY)
T.total= nrow(By)
Ik  = diag(rep(1, kk))
# 16 coef if 4 variables; 55 coef if 5 variables
Coef = t(VAR.EbyE(ddY, ddX, CONST)$ddA)              
# residuals
U    = VAR.EbyE(ddY, ddX, CONST)$ddU
BSigma.u = VAR.ddSigma.OLS(ddY, ddX, CONST)
if(CONST == TRUE){
  const = Coef[, ncol(Coef)]
  Coef.noc= Coef[,-ncol(Coef)]                      
}else{
  const = matrix(0, kk, 1)
  Coef.noc = Coef
}

Theta.unit= VAR.Theta(Coef, h, BSigma.u, CONST)$unit # estimate Theta.unit
Theta.std = VAR.Theta(Coef, h, BSigma.u, CONST)$std  # estimate Theta.std

# dm.U <- U-mean(U)
dm.U <- U

N = 2000 
Theta.unit.sim = vector("list", N)
Theta.std.sim  = vector("list", N)

# check dimension
print("check dimensionality")
dim(ddX); dim(Coef.noc); dim(dm.U)


df_IRF.sim <- array(NA, c(hrz+1,kk^2,N)) #dimensions are: Time Period, Number of shock interacts with variables, page (number of Bootstrap resamplings)
counter <- 1
while(TRUE){
  
  #cat("Now, there are ", counter-1, " sets of resamples.\n")
  Y.sim = matrix(0, nrow = T.total, ncol = kk)          # Y.sim = 0 #pseudo time series
  Y.sim[c(1:VAR.P),] = By[c(1:VAR.P), ] #initial values
  
  boot.number = sample(c(1:T), replace = TRUE)      
  U.sim = dm.U[boot.number,]
  
  # predicted values given the above initial values
  last.y= c(t(By[VAR.P:1,]))
  for(ii in 1:T){
    last.y = last.y[1:(kk*VAR.P)]
    Y.sim[ii+VAR.P, ] = Coef.noc %*% last.y + const + U.sim[ii,]      
    last.y = c(Y.sim[ii+VAR.P,], last.y)
  }
  
  #   Y.sim[-c(1:VAR.P),] <- matrix(const, nrow = T.total-VAR.P, ncol = kk, byrow = T) + ddX %*% t(Coef.noc) + U.sim
  
  
  #`Y.sim` is the pseudo time series
  
  ### SVAR.sim Start ###
  
  Y_pseudo     = VAR.Y(Y.sim, VAR.P)        
  X_pseudo     = VAR.X(Y.sim, VAR.P)        
  Coef.OLS_pseudo    = VAR.OLS(Y_pseudo, X_pseudo, CONST)
  Sigma.OLS_pseudo   = VAR.Sigma.OLS(Y_pseudo, X_pseudo, Coef.OLS_pseudo, CONST)
  C.Prime_pseudo <- chol(Sigma.OLS_pseudo)
  C_pseudo <- t(C.Prime_pseudo)
  B0_pseudo <- diag(diag(C_pseudo), ncol = kk, nrow = kk)
  A0_pseudo <- B0_pseudo %*% solve(C_pseudo)
  SVAR_AB_est.sim <- list("A0.svar" = A0_pseudo, "B0.svar" = B0_pseudo)
  SVAR_AB_IRF.sim <- VAR.svarirf.AB(Y.sim, VAR.P, Amat, Bmat, h = hrz, CONST, SVAR_AB_est = SVAR_AB_est.sim)
  
  df_IRF_plot.sim <- matrix(NA, hrz+1, kk^2) #%>% as.tibble()
  # df_IRF.sim <- array(1:(120*25*N), c(120,25,N))
  # df_IRF.sim[2,1,1] # slicing
  
  h <- 0 
  for(period in SVAR_AB_IRF.sim){
    k <- 0 
    h <- h+1 
    for(j in 1:kk){
      for(i in 1:kk){
        k <- k+1 
        df_IRF_plot.sim[h,k] <- period[i,j]
      }
    }
  }
  df_IRF.sim[,,counter] <- df_IRF_plot.sim
  ### SVAR.sim Ends ###
  if(counter>=N){
    break
  }
  counter <- counter+1
}
# Save
saveRDS(df_IRF.sim, file = "df_IRF.sim.rds")

df_IRF.sim <- read_rds("df_IRF.sim.rds")

head(df_IRF.sim[,,1000])
print(sum(is.na(df_IRF.sim)))

# IRF & Bootstrap C.I.
df_IRF_plot.BS.L <- matrix(NA, nrow = hrz+1, ncol = kk^2)
df_IRF_plot.BS.U <- matrix(NA, nrow = hrz+1, ncol = kk^2)
df_IRF_plot.BS.Median <- matrix(NA, nrow = hrz+1, ncol = kk^2)
df_IRF_plot.BS.Mean <- matrix(NA, nrow = hrz+1, ncol = kk^2)
for(col in 1:(kk^2)){
  for(row in 1:(hrz+1) ){
    df_IRF_plot.BS.L[row,col] <- quantile(df_IRF.sim[row,col,], probs = 0.025)
    df_IRF_plot.BS.U[row,col] <- quantile(df_IRF.sim[row,col,], probs = 0.975)
    df_IRF_plot.BS.Median[row,col] <- quantile(df_IRF.sim[row,col,], probs = 0.5)
    df_IRF_plot.BS.Mean[row,col] <- mean(df_IRF.sim[row,col,])
  }
}

df_IRF_plot.BS.L <- df_IRF_plot.BS.L %>% as_tibble()
df_IRF_plot.BS.U <- df_IRF_plot.BS.U %>% as_tibble()
df_IRF_plot.BS.Median <- df_IRF_plot.BS.Median %>% as_tibble()
df_IRF_plot.BS.Mean <- df_IRF_plot.BS.Mean %>% as_tibble()

ind <- 0
for(i in 1:kk){
  for(j in 1:kk){
    ind <- ind+1
    nam <- paste("shock", j, "y", i, sep = '')
    assign(nam, bind_cols(df_IRF_plot.BS.L[ind], df_IRF_plot.BS.U[ind],
                          df_IRF_plot.BS.Median[ind], df_IRF_plot.BS.Mean[ind],
                          df_IRF_plot[ind]))
    evalStr <- paste0("colnames(", nam, ") <- c('Lower', 'Upper', 'Median', 'Mean', 'Actual')")
    eval(parse(text=evalStr))
    evalStr <- paste0("p", ind, " <- ", "ggplot(",nam,") +geom_hline(yintercept=0, color = 'grey')+ geom_line(aes(x = 1:nrow(", nam, "), y = Lower), linetype = 'dashed', col='red')+geom_line(aes(x = 1:nrow(", nam, "), y = Upper), linetype = 'dashed', col='red')+geom_line(aes(x = 1:nrow(", nam, "), y = Median), col = 'Blue')")
    eval(parse(text=evalStr))
  }
}  

Text_Size_Theme = theme(
  axis.title.x = element_text(size = 12),
  axis.text.x = element_text(size = 12),
  axis.title.y = element_text(size = 12),
  axis.title = element_text(size = 12),
  plot.title = element_text(size=12))

add_layout <- function(p = p1,
                       response_of = "this variable", react_to = "the shock"){
  title_text <- paste("Response of", response_of, "to", react_to, sep = ' ')
  plot.return <- p + labs(x = "Period", y = "", title = title_text) + Text_Size_Theme
  return(plot.return)
}

vlb_list = c("Stock Price",
             "House Price",
             "CPI"
)
p1 <- add_layout(p1, vlb_list[1], "Stock shock")
p2 <- add_layout(p2, vlb_list[2], "Stock shock")
p3 <- add_layout(p3, vlb_list[3], "Stock shock")

p4 <- add_layout(p4, vlb_list[1], "Property Price shock")
p5 <- add_layout(p5, vlb_list[2], "Property Price shock")
p6 <- add_layout(p6, vlb_list[3], "Property Price shock")

p7 <- add_layout(p7, vlb_list[1], "CPI shock")
p8 <- add_layout(p8, vlb_list[2], "CPI shock")
p9 <- add_layout(p9, vlb_list[3], "CPI shock")

multiplot(p1,p2,p3,p4,p5,p6,
          p7,p8,p9,cols = 3)

multiplot(p1,p2,p3, cols = 1)
ggsave(filename = "/Users/caiyawei/Desktop/VAR/IRF_shock_cpi1.png", 
       plot = multiplot(p1,p2,p3, cols = 1),
       width = 30, height = 20, units = "cm",
       device = "png")

multiplot(p4,p5,p6, cols = 1)
ggsave(filename = "/Users/caiyawei/Desktop/VAR/IRF_shock_cpi2.png", 
       plot = multiplot(p4,p5,p6, cols = 1),
       width = 30, height = 20, units = "cm",
       device = "png")

multiplot(p7,p8,p9, cols = 1)
ggsave(filename = "/Users/caiyawei/Desktop/VAR/IRF_shock_cpi3.png", 
       plot = multiplot(p7,p8,p9, cols = 1),
       width = 30, height = 20, units = "cm",
       device = "png")

### Obtain the exogenous shock###
Coef.noc
Coef.OLS
const

# check the order of variables is correct (stock shock) (cpi)
lm(stock_price ~ lag(stock_price)+lag(property_price)+lag(cpi)+
     +lag(stock_price, n=2)+lag(property_price, n=2)+lag(cpi, n=2)+lag(stock_price, n=3)+lag(property_price, n=3)+lag(cpi, n=3)+lag(stock_price, n=4)+lag(property_price, n=4)+lag(cpi, n=4),
   data = data_cpi)

get_residual <- function(vlb_order = 1){
  X.Mat <- By %>% as_tibble() %>%
    mutate(V1 = lag(stock_price),
           V2 = lag(property_price),
           V3 = lag(cpi),
           V4 = lag(stock_price, n = 2),
           V5 = lag(property_price, n = 2),
           V6 = lag(cpi, n = 2),
           V7 = lag(stock_price, n = 3),
           V8 = lag(property_price, n = 3),
           V9 = lag(cpi, n = 3),
           V10 = lag(stock_price, n = 4),
           V11 = lag(property_price, n = 4),
           V12 = lag(cpi, n = 4),
           V13 = rep(1, nrow(By))
    ) %>%
    select(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13
    ) %>%
    drop_na() %>%
    as.matrix
  Y.Mat <- By[(VAR.P+1):nrow(By),vlb_order] %>% as.matrix
  Y.Mat_hat <- X.Mat %*% (Coef.OLS[vlb_order, ] %>% as.matrix)
  residual <- Y.Mat - Y.Mat_hat
  return(residual)
}

# resid1 <- get_residual(1)
for(i in 1:3){
  assign(paste0("resid_",i), get_residual(i))
}

resid.Mat <- cbind(resid_1, resid_2, resid_3)
error.Mat <- (solve(B0) %*% A0 %*% t(resid.Mat)) %>% t()
error.Mat[,1] %>% hist

# Extract the stock shock

error.Mat %>% as_tibble() %>%
  select(e1 = 1, e2 = 2, e3 = 3) %>%
  mutate(no = row_number()) %>%
  ggplot()+
  geom_line(aes(no, e1))

data.err <- data_cpi %>% 
  bind_cols(e1 = c(rep(NA, VAR.P), error.Mat[,1])) %>%
  mutate(e1_pos = fifelse(e1 > 0, e1, 0),
         e1_neg = fifelse(e1 < 0, e1, 0))
data.err %>% head(10)

data.err %>% ggplot()+
  geom_line(aes(year, e1))

# check the order of variables is correct (house shock) (cpi)
lm(property_price ~ lag(stock_price)+lag(property_price)+lag(cpi)+
     +lag(stock_price, n=2)+lag(property_price, n=2)+lag(cpi, n=2)+lag(stock_price, n=3)+lag(property_price, n=3)+lag(cpi, n=3)+lag(stock_price, n=4)+lag(property_price, n=4)+lag(cpi, n=4),
   data = data_cpi)

get_residual <- function(vlb_order = 1){
  X.Mat <- By %>% as_tibble() %>%
    mutate(V1 = lag(stock_price),
           V2 = lag(property_price),
           V3 = lag(cpi),
           V4 = lag(stock_price, n = 2),
           V5 = lag(property_price, n = 2),
           V6 = lag(cpi, n = 2),
           V7 = lag(stock_price, n = 3),
           V8 = lag(property_price, n = 3),
           V9 = lag(cpi, n = 3),
           V10 = lag(stock_price, n = 4),
           V11 = lag(property_price, n = 4),
           V12 = lag(cpi, n = 4),
           V13 = rep(1, nrow(By))
    ) %>%
    select(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13
    ) %>%
    drop_na() %>%
    as.matrix
  Y.Mat <- By[(VAR.P+1):nrow(By),vlb_order] %>% as.matrix
  Y.Mat_hat <- X.Mat %*% (Coef.OLS[vlb_order, ] %>% as.matrix)
  residual <- Y.Mat - Y.Mat_hat
  return(residual)
}

# resid1 <- get_residual(1)
for(i in 1:3){
  assign(paste0("resid_",i), get_residual(i))
}

resid.Mat <- cbind(resid_1, resid_2, resid_3)
error.Mat <- (solve(B0) %*% A0 %*% t(resid.Mat)) %>% t()
error.Mat[,2] %>% hist

# Extract the house shock

error.Mat %>% as_tibble() %>%
  select(e1 = 1, e2 = 2, e3 = 3) %>%
  mutate(no = row_number()) %>%
  ggplot()+
  geom_line(aes(no, e2))

data.err2 <- data_cpi %>% 
  bind_cols(e2 = c(rep(NA, VAR.P), error.Mat[,2])) %>%
  mutate(e2_pos = fifelse(e2 > 0, e2, 0),
         e2_neg = fifelse(e2 < 0, e2, 0))
data.err2 %>% head(10)

data.err2 %>% ggplot()+
  geom_line(aes(year, e2))


data.err["e2"] = data.err2$e2
### Regression ###
reg3 <- lm(cpi_d ~ Period + e1 + e2 + interest_rate_d, data = data.err)
summary(reg3)

data.err["e2"] = data.err["e2"] - data.err$e2[6]
data.err["e1"] = data.err["e1"] - data.err$e1[6]
reg3 <- lm(cpi_d ~ Period + e1 + e2 + interest_rate_d, data = data.err)
reg3 <- lm(cpi_d ~ e1 + e2 + interest_rate_d, data = data.err)
reg3 <- lm(cpi_d ~ interest_rate_d, data = data.err)
reg3 <- lm(cpi ~ interest_rate, data = data.err)
summary(reg3)








### DATA GDP ###
file = "/Users/caiyawei/Desktop/Japan Data/SVAR_forR.xlsx"
data_gdp = read_excel(file)
data_gdp = data_gdp[-c(107),]



By <- data_gdp %>% 
  select(stock_price, property_price, gdp) %>% 
  as.matrix

### RAW DATA ###
dim(By)
kk <- dim(By)[2]
draw_raw <- function(vlb_name = ""){
  vlb_name <- sym(vlb_name)
  plot.return <- data_gdp %>%
    ggplot()+
    geom_line(aes(year, !!vlb_name))
  return(plot.return)
}

p1 <- draw_raw("stock_price")
p2 <- draw_raw("property_price")
p3 <- draw_raw("gdp")

multiplot(p1, p2, p3)

### Setting Model ###
VAR.P = 4                       
CONST = TRUE                    
Y     = VAR.Y(By, VAR.P)        
X     = VAR.X(By, VAR.P)        

hrz=39 # the length of response

### Reduced Form VAR ###
(Coef.OLS    = VAR.OLS(Y, X, CONST)                  )
(Sigma.OLS   = VAR.Sigma.OLS(Y, X, Coef.OLS, CONST)  )
(Sigma.MLE   = VAR.Sigma.MLE(Y, X, Coef.OLS, CONST))

VAR.P = 4 # By AIC

### SVAR ###
### 6-variable model

Amat = diag(3)
# Identification Conditions

Amat[2,1]  = NA; 
Amat[3,1]  = NA; Amat[3,2]  = NA; 


Bmat = diag(3)
diag(Bmat) = NA

Amat;Bmat

### Estimate A^hat and B^hat ###
C.Prime <- chol(Sigma.OLS)
C <- t(C.Prime)
C

### Solving system of linear equations ###
B0 <- diag(diag(C), ncol = 3, nrow = 3)
B0

A0 <- B0 %*% solve(C)
A0

SVAR_AB_est <- list("A0.svar" = A0, "B0.svar" = B0)

### IRF Without Bootstrap CI ###
### IRF
SVAR_AB_IRF <- VAR.svarirf.AB(By, VAR.P, Amat, Bmat, h = hrz, CONST, SVAR_AB_est = SVAR_AB_est)


df_IRF_plot <- matrix(NA, hrz+1, kk^2) #%>% as.tibble() ## hrz+1
#dim(df_IRF_plot)
h <- 0 
for(period in SVAR_AB_IRF){
  k <- 0 
  h <- h+1 
  for(j in 1:kk){
    for(i in 1:kk){
      k <- k+1 
      df_IRF_plot[h,k] <- period[i,j]
    }
  }
}
df_IRF_plot <- df_IRF_plot %>% as_tibble()

kk*1:kk

#output entire table
IRF_TABLE <- df_IRF_plot[,kk*1:kk] %>% select(stock_price = 1,
                                              property_price = 2,
                                              gdp = 3)
write.table(IRF_TABLE, file = "/Users/caiyawei/Desktop/IRF_TABLE2.csv", sep = ",", row.names = FALSE)

draw_IRF <- function(df = df_IRF_plot, V1 = 1){
  V1 <- paste0("V", V1) %>% sym()
  plot.return <- ggplot(df) + 
    geom_line(aes(x = 1:nrow(df), y = !!V1))
  return(plot.return)
}

for(i in 1:(kk^2)){
  assign(paste0("p",i), draw_IRF(df_IRF_plot, i))
}

multiplot(p1,p2,p3,p4,p5,p6,
          p7,p8,p9,cols = 3)

### IRF (Bootstrap C.I.) ###
lower = 0.025                                        # 95% CI
upper = 1-lower
kk = ncol(By)
ddY = VAR.ddY(By, VAR.P)
ddX = VAR.ddX(By, VAR.P)

# dim(ddY); dim(ddX)

T   = nrow(ddY)
T.total= nrow(By)
Ik  = diag(rep(1, kk))
# 16 coef if 4 variables; 55 coef if 5 variables
Coef = t(VAR.EbyE(ddY, ddX, CONST)$ddA)          
# residuals
U    = VAR.EbyE(ddY, ddX, CONST)$ddU
BSigma.u = VAR.ddSigma.OLS(ddY, ddX, CONST)
if(CONST == TRUE){
  const = Coef[, ncol(Coef)]
  Coef.noc= Coef[,-ncol(Coef)]                   
}else{
  const = matrix(0, kk, 1)
  Coef.noc = Coef
}

Theta.unit= VAR.Theta(Coef, h, BSigma.u, CONST)$unit 
Theta.std = VAR.Theta(Coef, h, BSigma.u, CONST)$std  

# dm.U <- U-mean(U)
dm.U <- U

N = 2000 
Theta.unit.sim = vector("list", N)
Theta.std.sim  = vector("list", N)

# check dimension
print("check dimensionality")
dim(ddX); dim(Coef.noc); dim(dm.U)


df_IRF.sim <- array(NA, c(hrz+1,kk^2,N)) #dimensions are: Time Period, Number of shock interacts with variables, page (number of Bootstrap resamplings)
counter <- 1
while(TRUE){
  
  #cat("Now, there are ", counter-1, " sets of resamples.\n")
  Y.sim = matrix(0, nrow = T.total, ncol = kk)          # Y.sim = 0 #pseudo time series
  Y.sim[c(1:VAR.P),] = By[c(1:VAR.P), ] #initial values
  
  boot.number = sample(c(1:T), replace = TRUE)      
  U.sim = dm.U[boot.number,]
  
  # predicted values given the above initial values
  last.y= c(t(By[VAR.P:1,]))
  for(ii in 1:T){
    last.y = last.y[1:(kk*VAR.P)]
    Y.sim[ii+VAR.P, ] = Coef.noc %*% last.y + const + U.sim[ii,]      
    last.y = c(Y.sim[ii+VAR.P,], last.y)
  }
  
  #   Y.sim[-c(1:VAR.P),] <- matrix(const, nrow = T.total-VAR.P, ncol = kk, byrow = T) + ddX %*% t(Coef.noc) + U.sim
  
  
  #`Y.sim` is the pseudo time series
  
  ### SVAR.sim Start ###
  
  Y_pseudo     = VAR.Y(Y.sim, VAR.P)        
  X_pseudo     = VAR.X(Y.sim, VAR.P)        
  Coef.OLS_pseudo    = VAR.OLS(Y_pseudo, X_pseudo, CONST)
  Sigma.OLS_pseudo   = VAR.Sigma.OLS(Y_pseudo, X_pseudo, Coef.OLS_pseudo, CONST)
  C.Prime_pseudo <- chol(Sigma.OLS_pseudo)
  C_pseudo <- t(C.Prime_pseudo)
  B0_pseudo <- diag(diag(C_pseudo), ncol = kk, nrow = kk)
  A0_pseudo <- B0_pseudo %*% solve(C_pseudo)
  SVAR_AB_est.sim <- list("A0.svar" = A0_pseudo, "B0.svar" = B0_pseudo)
  SVAR_AB_IRF.sim <- VAR.svarirf.AB(Y.sim, VAR.P, Amat, Bmat, h = hrz, CONST, SVAR_AB_est = SVAR_AB_est.sim)
  
  df_IRF_plot.sim <- matrix(NA, hrz+1, kk^2) #%>% as.tibble()
  # df_IRF.sim <- array(1:(120*25*N), c(120,25,N))
  # df_IRF.sim[2,1,1] # slicing
  
  h <- 0 
  for(period in SVAR_AB_IRF.sim){
    k <- 0 
    h <- h+1 
    for(j in 1:kk){
      for(i in 1:kk){
        k <- k+1 
        df_IRF_plot.sim[h,k] <- period[i,j]
      }
    }
  }
  df_IRF.sim[,,counter] <- df_IRF_plot.sim
  ### SVAR.sim Ends ###
  if(counter>=N){
    break
  }
  counter <- counter+1
}
# Save
saveRDS(df_IRF.sim, file = "df_IRF.sim.rds")

df_IRF.sim <- read_rds("df_IRF.sim.rds")

head(df_IRF.sim[,,1000])
print(sum(is.na(df_IRF.sim)))

# Draw IRF & Bootstrap C.I.
df_IRF_plot.BS.L <- matrix(NA, nrow = hrz+1, ncol = kk^2)
df_IRF_plot.BS.U <- matrix(NA, nrow = hrz+1, ncol = kk^2)
df_IRF_plot.BS.Median <- matrix(NA, nrow = hrz+1, ncol = kk^2)
df_IRF_plot.BS.Mean <- matrix(NA, nrow = hrz+1, ncol = kk^2)
for(col in 1:(kk^2)){
  for(row in 1:(hrz+1) ){
    df_IRF_plot.BS.L[row,col] <- quantile(df_IRF.sim[row,col,], probs = 0.025)
    df_IRF_plot.BS.U[row,col] <- quantile(df_IRF.sim[row,col,], probs = 0.975)
    df_IRF_plot.BS.Median[row,col] <- quantile(df_IRF.sim[row,col,], probs = 0.5)
    df_IRF_plot.BS.Mean[row,col] <- mean(df_IRF.sim[row,col,])
  }
}

df_IRF_plot.BS.L <- df_IRF_plot.BS.L %>% as_tibble()
df_IRF_plot.BS.U <- df_IRF_plot.BS.U %>% as_tibble()
df_IRF_plot.BS.Median <- df_IRF_plot.BS.Median %>% as_tibble()
df_IRF_plot.BS.Mean <- df_IRF_plot.BS.Mean %>% as_tibble()

ind <- 0
for(i in 1:kk){
  for(j in 1:kk){
    ind <- ind+1
    nam <- paste("shock", j, "y", i, sep = '')
    assign(nam, bind_cols(df_IRF_plot.BS.L[ind], df_IRF_plot.BS.U[ind],
                          df_IRF_plot.BS.Median[ind], df_IRF_plot.BS.Mean[ind],
                          df_IRF_plot[ind]))
    evalStr <- paste0("colnames(", nam, ") <- c('Lower', 'Upper', 'Median', 'Mean', 'Actual')")
    eval(parse(text=evalStr))
    evalStr <- paste0("p", ind, " <- ", "ggplot(",nam,") +geom_hline(yintercept=0, color = 'grey')+ geom_line(aes(x = 1:nrow(", nam, "), y = Lower), linetype = 'dashed', col='red')+geom_line(aes(x = 1:nrow(", nam, "), y = Upper), linetype = 'dashed', col='red')+geom_line(aes(x = 1:nrow(", nam, "), y = Median), col = 'Blue')")
    eval(parse(text=evalStr))
  }
}  

Text_Size_Theme = theme(
  axis.title.x = element_text(size = 12),
  axis.text.x = element_text(size = 12),
  axis.title.y = element_text(size = 12),
  axis.title = element_text(size = 12),
  plot.title = element_text(size=12))

add_layout <- function(p = p1,
                       response_of = "this variable", react_to = "the shock"){
  title_text <- paste("Response of", response_of, "to", react_to, sep = ' ')
  plot.return <- p + labs(x = "Period", y = "", title = title_text) + Text_Size_Theme
  return(plot.return)
}

vlb_list = c("Stock Price",
             "House Price",
             "GDP"
)
p1 <- add_layout(p1, vlb_list[1], "Stock shock")
p2 <- add_layout(p2, vlb_list[2], "Stock shock")
p3 <- add_layout(p3, vlb_list[3], "Stock shock")

p4 <- add_layout(p4, vlb_list[1], "Property Price shock")
p5 <- add_layout(p5, vlb_list[2], "Property Price shock")
p6 <- add_layout(p6, vlb_list[3], "Property Price shock")

p7 <- add_layout(p7, vlb_list[1], "GDP shock")
p8 <- add_layout(p8, vlb_list[2], "GDP shock")
p9 <- add_layout(p9, vlb_list[3], "GDP shock")

multiplot(p1,p2,p3,p4,p5,p6,
          p7,p8,p9,cols = 3)

multiplot(p1,p2,p3, cols = 1)
ggsave(filename = "/Users/caiyawei/Desktop/VAR/IRF_shock_gdp1.png", 
       plot = multiplot(p1,p2,p3, cols = 1),
       width = 30, height = 20, units = "cm",
       device = "png")

multiplot(p4,p5,p6, cols = 1)
ggsave(filename = "/Users/caiyawei/Desktop/VAR/IRF_shock_gdp2.png", 
       plot = multiplot(p4,p5,p6, cols = 1),
       width = 30, height = 20, units = "cm",
       device = "png")

multiplot(p7,p8,p9, cols = 1)
ggsave(filename = "/Users/caiyawei/Desktop/VAR/IRF_shock_gdp33.png", 
       plot = multiplot(p7,p8,p9, cols = 1),
       width = 30, height = 20, units = "cm",
       device = "png")

### Obtain the exogenous shock###
Coef.noc
Coef.OLS
const

# check the order of variables is correct (stock shock) (gdp)
lm(stock_price ~ lag(stock_price)+lag(property_price)+lag(gdp)+
     +lag(stock_price, n=2)+lag(property_price, n=2)+lag(gdp, n=2)+lag(stock_price, n=3)+lag(property_price, n=3)+lag(gdp, n=3)+lag(stock_price, n=4)+lag(property_price, n=4)+lag(gdp, n=4),
   data = data_gdp)

get_residual <- function(vlb_order = 1){
  X.Mat <- By %>% as_tibble() %>%
    mutate(V1 = lag(stock_price),
           V2 = lag(property_price),
           V3 = lag(gdp),
           V4 = lag(stock_price, n = 2),
           V5 = lag(property_price, n = 2),
           V6 = lag(gdp, n = 2),
           V7 = lag(stock_price, n = 3),
           V8 = lag(property_price, n = 3),
           V9 = lag(gdp, n = 3),
           V10 = lag(stock_price, n = 4),
           V11 = lag(property_price, n = 4),
           V12 = lag(gdp, n = 4),
           V13 = rep(1, nrow(By))
    ) %>%
    select(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13
    ) %>%
    drop_na() %>%
    as.matrix
  Y.Mat <- By[(VAR.P+1):nrow(By),vlb_order] %>% as.matrix
  Y.Mat_hat <- X.Mat %*% (Coef.OLS[vlb_order, ] %>% as.matrix)
  residual <- Y.Mat - Y.Mat_hat
  return(residual)
}

# resid1 <- get_residual(1)
for(i in 1:3){
  assign(paste0("resid_",i), get_residual(i))
}

resid.Mat <- cbind(resid_1, resid_2, resid_3)
error.Mat <- (solve(B0) %*% A0 %*% t(resid.Mat)) %>% t()
error.Mat[,1] %>% hist

# Extract the stock shock

error.Mat %>% as_tibble() %>%
  select(e1 = 1, e2 = 2, e3 = 3) %>%
  mutate(no = row_number()) %>%
  ggplot()+
  geom_line(aes(no, e1))

data.err <- data_gdp %>% 
  bind_cols(e1 = c(rep(NA, VAR.P), error.Mat[,1])) %>%
  mutate(e1_pos = fifelse(e1 > 0, e1, 0),
         e1_neg = fifelse(e1 < 0, e1, 0))
data.err %>% head(10)

data.err %>% ggplot()+
  geom_line(aes(year, e1))

# check the order of variables is correct (house shock) (gdp)
lm(property_price ~ lag(stock_price)+lag(property_price)+lag(gdp)+
     +lag(stock_price, n=2)+lag(property_price, n=2)+lag(gdp, n=2)+lag(stock_price, n=3)+lag(property_price, n=3)+lag(gdp, n=3)+lag(stock_price, n=4)+lag(property_price, n=4)+lag(gdp, n=4),
   data = data_gdp)

get_residual <- function(vlb_order = 1){
  X.Mat <- By %>% as_tibble() %>%
    mutate(V1 = lag(stock_price),
           V2 = lag(property_price),
           V3 = lag(gdp),
           V4 = lag(stock_price, n = 2),
           V5 = lag(property_price, n = 2),
           V6 = lag(gdp, n = 2),
           V7 = lag(stock_price, n = 3),
           V8 = lag(property_price, n = 3),
           V9 = lag(gdp, n = 3),
           V10 = lag(stock_price, n = 4),
           V11 = lag(property_price, n = 4),
           V12 = lag(gdp, n = 4),
           V13 = rep(1, nrow(By))
    ) %>%
    select(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13
    ) %>%
    drop_na() %>%
    as.matrix
  Y.Mat <- By[(VAR.P+1):nrow(By),vlb_order] %>% as.matrix
  Y.Mat_hat <- X.Mat %*% (Coef.OLS[vlb_order, ] %>% as.matrix)
  residual <- Y.Mat - Y.Mat_hat
  return(residual)
}

# resid1 <- get_residual(1)
for(i in 1:3){
  assign(paste0("resid_",i), get_residual(i))
}

resid.Mat <- cbind(resid_1, resid_2, resid_3)
error.Mat <- (solve(B0) %*% A0 %*% t(resid.Mat)) %>% t()
error.Mat[,2] %>% hist

# Extract the house shock

error.Mat %>% as_tibble() %>%
  select(e1 = 1, e2 = 2, e3 = 3) %>%
  mutate(no = row_number()) %>%
  ggplot()+
  geom_line(aes(no, e2))

data.err2 <- data_gdp %>% 
  bind_cols(e2 = c(rep(NA, VAR.P), error.Mat[,2])) %>%
  mutate(e2_pos = fifelse(e2 > 0, e2, 0),
         e2_neg = fifelse(e2 < 0, e2, 0))
data.err2 %>% head(10)

data.err2 %>% ggplot()+
  geom_line(aes(year, e2))


data.err["e2"] = data.err2$e2
### Regression ###
reg2 <- lm(gdp_d~ Period + e1 + e2 + interest_rate_d, data = data.err)
summary(reg2)

### substract 1994Q4###
data.err["e2"] = data.err["e2"] - data.err$e2[6]
data.err["e1"] = data.err["e1"] - data.err$e1[6]
reg2 <- lm(gdp_d~ Period + e1 + e2 + interest_rate_d, data = data.err)
summary(reg2)

texreg(list(reg1, reg2, reg3), caption = "Three SVAR models.", digits = 5)




###Order2 House -> Stock -> Macro###
file = "/Users/caiyawei/Desktop/Japan Data/SVAR_forR.xlsx"
data = read_excel(file)
data = data[-c(107),]

By <- data %>% 
  select(property_price, stock_price, employment) %>% 
  as.matrix

### RAW DATA ###
dim(By)
kk <- dim(By)[2]
draw_raw <- function(vlb_name = ""){
  vlb_name <- sym(vlb_name)
  plot.return <- data %>%
    ggplot()+
    geom_line(aes(year, !!vlb_name))
  return(plot.return)
}

p1 <- draw_raw("property_price")
p2 <- draw_raw("stock_price")
p3 <- draw_raw("employment")

multiplot(p1, p2, p3)

### Setting Model ###
VAR.P = 4                       # 最大的落後項數
CONST = TRUE                    # 是否有常數項
Y     = VAR.Y(By, VAR.P)        # 設定 Y
X     = VAR.X(By, VAR.P)        # 設定 X

hrz=39 # the length of response

### Reduced Form VAR ###
(Coef.OLS    = VAR.OLS(Y, X, CONST)                  )
(Sigma.OLS   = VAR.Sigma.OLS(Y, X, Coef.OLS, CONST)  )
(Sigma.MLE   = VAR.Sigma.MLE(Y, X, Coef.OLS, CONST))

VAR.P = 4 # By AIC

### SVAR ###
### 6-variable model

Amat = diag(3)
# Identification Conditions

Amat[2,1]  = NA; 
Amat[3,1]  = NA; Amat[3,2]  = NA; 


Bmat = diag(3)
diag(Bmat) = NA

Amat;Bmat

### Estimate A^hat and B^hat ###
C.Prime <- chol(Sigma.OLS)
C <- t(C.Prime)
C

### Solving system of linear equations ###
B0 <- diag(diag(C), ncol = 3, nrow = 3)
B0

A0 <- B0 %*% solve(C)
A0

SVAR_AB_est <- list("A0.svar" = A0, "B0.svar" = B0)

### IRF Without Bootstrap CI ###
### IRF
SVAR_AB_IRF <- VAR.svarirf.AB(By, VAR.P, Amat, Bmat, h = hrz, CONST, SVAR_AB_est = SVAR_AB_est)


# 6*6個圖的time series
df_IRF_plot <- matrix(NA, hrz+1, kk^2) #%>% as.tibble() ## hrz+1
#dim(df_IRF_plot)
h <- 0 # h表示第幾期的IRF
for(period in SVAR_AB_IRF){
  k <- 0 # k表示把5*5的矩陣攤平到25個col的df時，要攤到第幾個columns上
  h <- h+1 # h表示第幾期的IRF
  for(j in 1:kk){
    for(i in 1:kk){
      k <- k+1 # k表示把5*5的矩陣攤平到25個col的df時，要攤到第幾個columns上
      df_IRF_plot[h,k] <- period[i,j]
    }
  }
}
df_IRF_plot <- df_IRF_plot %>% as_tibble()

kk*1:kk

#output entire table
IRF_TABLE <- df_IRF_plot[,kk*1:kk] %>% select(property_price = 1,
                                              stock_price = 2,
                                              employment = 3)
write.table(IRF_TABLE, file = "/Users/caiyawei/Desktop/IRF_TABLE_order2.csv", sep = ",", row.names = FALSE)

draw_IRF <- function(df = df_IRF_plot, V1 = 1){
  V1 <- paste0("V", V1) %>% sym()
  plot.return <- ggplot(df) + 
    geom_line(aes(x = 1:nrow(df), y = !!V1))
  return(plot.return)
}

for(i in 1:(kk^2)){
  assign(paste0("p",i), draw_IRF(df_IRF_plot, i))
}

multiplot(p1,p2,p3,p4,p5,p6,
          p7,p8,p9,cols = 3)

### IRF (Bootstrap C.I.) ###
lower = 0.025                                        # 控制成 95% CI
upper = 1-lower
kk = ncol(By)
ddY = VAR.ddY(By, VAR.P)
ddX = VAR.ddX(By, VAR.P)

# dim(ddY); dim(ddX)

T   = nrow(ddY)
T.total= nrow(By)
Ik  = diag(rep(1, kk))
# 16 coef if 4 variables; 55 coef if 5 variables
Coef = t(VAR.EbyE(ddY, ddX, CONST)$ddA)              # Step 1 估計模型
# residuals
U    = VAR.EbyE(ddY, ddX, CONST)$ddU
BSigma.u = VAR.ddSigma.OLS(ddY, ddX, CONST)
if(CONST == TRUE){
  const = Coef[, ncol(Coef)]
  Coef.noc= Coef[,-ncol(Coef)]                      # 刪掉 const
}else{
  const = matrix(0, kk, 1)
  Coef.noc = Coef
}

Theta.unit= VAR.Theta(Coef, h, BSigma.u, CONST)$unit # 估算 Theta.unit
Theta.std = VAR.Theta(Coef, h, BSigma.u, CONST)$std  # 估算 Theta.std

# dm.U <- U-mean(U)
dm.U <- U

N = 2000 #重抽次數
Theta.unit.sim = vector("list", N)
Theta.std.sim  = vector("list", N)

# check dimension
print("check dimensionality")
dim(ddX); dim(Coef.noc); dim(dm.U)


# 存N次重抽的IRF
df_IRF.sim <- array(NA, c(hrz+1,kk^2,N)) #dimensions are: Time Period, Number of shock interacts with variables, page (number of Bootstrap resamplings)
counter <- 1
while(TRUE){
  
  #cat("Now, there are ", counter-1, " sets of resamples.\n")
  Y.sim = matrix(0, nrow = T.total, ncol = kk)          # Y.sim = 0 #pseudo time series
  Y.sim[c(1:VAR.P),] = By[c(1:VAR.P), ] #initial values
  
  boot.number = sample(c(1:T), replace = TRUE)      # Step 3 取出放回
  U.sim = dm.U[boot.number,]
  
  # predicted values given the above initial values
  last.y= c(t(By[VAR.P:1,]))
  for(ii in 1:T){
    last.y = last.y[1:(kk*VAR.P)]
    Y.sim[ii+VAR.P, ] = Coef.noc %*% last.y + const + U.sim[ii,]      # Step 4 模擬資料
    last.y = c(Y.sim[ii+VAR.P,], last.y)
  }
  
  #   Y.sim[-c(1:VAR.P),] <- matrix(const, nrow = T.total-VAR.P, ncol = kk, byrow = T) + ddX %*% t(Coef.noc) + U.sim
  
  
  #`Y.sim` is the pseudo time series
  # Step 5 重新估算SVAR
  
  ### SVAR.sim Start ###
  
  Y_pseudo     = VAR.Y(Y.sim, VAR.P)        # 設定 Y
  X_pseudo     = VAR.X(Y.sim, VAR.P)        # 設定 X
  Coef.OLS_pseudo    = VAR.OLS(Y_pseudo, X_pseudo, CONST)
  Sigma.OLS_pseudo   = VAR.Sigma.OLS(Y_pseudo, X_pseudo, Coef.OLS_pseudo, CONST)
  C.Prime_pseudo <- chol(Sigma.OLS_pseudo)
  C_pseudo <- t(C.Prime_pseudo)
  B0_pseudo <- diag(diag(C_pseudo), ncol = kk, nrow = kk)
  A0_pseudo <- B0_pseudo %*% solve(C_pseudo)
  SVAR_AB_est.sim <- list("A0.svar" = A0_pseudo, "B0.svar" = B0_pseudo)
  SVAR_AB_IRF.sim <- VAR.svarirf.AB(Y.sim, VAR.P, Amat, Bmat, h = hrz, CONST, SVAR_AB_est = SVAR_AB_est.sim)
  
  # 5*5個圖的time series
  df_IRF_plot.sim <- matrix(NA, hrz+1, kk^2) #%>% as.tibble()
  # df_IRF.sim <- array(1:(120*25*N), c(120,25,N))
  # df_IRF.sim[2,1,1] # slicing
  
  h <- 0 # h表示第幾期的IRF
  for(period in SVAR_AB_IRF.sim){
    k <- 0 # k表示把5*5的矩陣攤平到25個col的df時，要攤到第幾個columns上
    h <- h+1 # h表示第幾期的IRF
    for(j in 1:kk){
      for(i in 1:kk){
        k <- k+1 # k表示把5*5的矩陣攤平到25個col的df時，要攤到第幾個columns上
        df_IRF_plot.sim[h,k] <- period[i,j]
      }
    }
  }
  # 把這一次重抽得到的IRF append進`df_IRF.sim`中
  df_IRF.sim[,,counter] <- df_IRF_plot.sim
  ### SVAR.sim Ends ###
  if(counter>=N){
    break
  }
  counter <- counter+1
}
# Save
saveRDS(df_IRF.sim, file = "df_IRF.sim.rds")

df_IRF.sim <- read_rds("df_IRF.sim.rds")

# 看某一頁
head(df_IRF.sim[,,1000])
print(sum(is.na(df_IRF.sim)))

# 畫IRF & Bootstrap C.I.
df_IRF_plot.BS.L <- matrix(NA, nrow = hrz+1, ncol = kk^2)
df_IRF_plot.BS.U <- matrix(NA, nrow = hrz+1, ncol = kk^2)
df_IRF_plot.BS.Median <- matrix(NA, nrow = hrz+1, ncol = kk^2)
df_IRF_plot.BS.Mean <- matrix(NA, nrow = hrz+1, ncol = kk^2)
for(col in 1:(kk^2)){
  for(row in 1:(hrz+1) ){
    df_IRF_plot.BS.L[row,col] <- quantile(df_IRF.sim[row,col,], probs = 0.025)
    df_IRF_plot.BS.U[row,col] <- quantile(df_IRF.sim[row,col,], probs = 0.975)
    df_IRF_plot.BS.Median[row,col] <- quantile(df_IRF.sim[row,col,], probs = 0.5)
    df_IRF_plot.BS.Mean[row,col] <- mean(df_IRF.sim[row,col,])
  }
}

df_IRF_plot.BS.L <- df_IRF_plot.BS.L %>% as_tibble()
df_IRF_plot.BS.U <- df_IRF_plot.BS.U %>% as_tibble()
df_IRF_plot.BS.Median <- df_IRF_plot.BS.Median %>% as_tibble()
df_IRF_plot.BS.Mean <- df_IRF_plot.BS.Mean %>% as_tibble()

ind <- 0
for(i in 1:kk){
  for(j in 1:kk){
    ind <- ind+1
    nam <- paste("shock", j, "y", i, sep = '')
    assign(nam, bind_cols(df_IRF_plot.BS.L[ind], df_IRF_plot.BS.U[ind],
                          df_IRF_plot.BS.Median[ind], df_IRF_plot.BS.Mean[ind],
                          df_IRF_plot[ind]))
    # 改名
    evalStr <- paste0("colnames(", nam, ") <- c('Lower', 'Upper', 'Median', 'Mean', 'Actual')")
    eval(parse(text=evalStr))
    # 圖層
    evalStr <- paste0("p", ind, " <- ", "ggplot(",nam,") +geom_hline(yintercept=0, color = 'grey')+ geom_line(aes(x = 1:nrow(", nam, "), y = Lower), linetype = 'dashed', col='red')+geom_line(aes(x = 1:nrow(", nam, "), y = Upper), linetype = 'dashed', col='red')+geom_line(aes(x = 1:nrow(", nam, "), y = Median), col = 'Blue')")
    eval(parse(text=evalStr))
  }
}  

Text_Size_Theme = theme(
  axis.title.x = element_text(size = 12),
  axis.text.x = element_text(size = 12),
  axis.title.y = element_text(size = 12),
  axis.title = element_text(size = 12),
  plot.title = element_text(size=12))

add_layout <- function(p = p1,
                       response_of = "this variable", react_to = "the shock"){
  title_text <- paste("Response of", response_of, "to", react_to, sep = ' ')
  plot.return <- p + labs(x = "Period", y = "", title = title_text) + Text_Size_Theme
  return(plot.return)
}

vlb_list = c("House Price",
             "Stock Price",
             "Employment"
)
p1 <- add_layout(p1, vlb_list[1], "Property Price shock")
p2 <- add_layout(p2, vlb_list[2], "Property Price shock")
p3 <- add_layout(p3, vlb_list[3], "Property Price shock")

p4 <- add_layout(p4, vlb_list[1], "Stock shock")
p5 <- add_layout(p5, vlb_list[2], "Stock shock")
p6 <- add_layout(p6, vlb_list[3], "Stock shock")

p7 <- add_layout(p7, vlb_list[1], "Employment shock")
p8 <- add_layout(p8, vlb_list[2], "Employment shock")
p9 <- add_layout(p9, vlb_list[3], "Employment shock")

multiplot(p1,p2,p3,p4,p5,p6,
          p7,p8,p9,cols = 3)

multiplot(p1,p2,p3, cols = 1)
ggsave(filename = "/Users/caiyawei/Desktop/VAR/IRF_shock1.png", 
       plot = multiplot(p1,p2,p3, cols = 1),
       width = 30, height = 20, units = "cm",
       device = "png")

multiplot(p4,p5,p6, cols = 1)
ggsave(filename = "/Users/caiyawei/Desktop/VAR/IRF_shock2.png", 
       plot = multiplot(p4,p5,p6, cols = 1),
       width = 30, height = 20, units = "cm",
       device = "png")

multiplot(p7,p8,p9, cols = 1)
ggsave(filename = "/Users/caiyawei/Desktop/VAR/IRF_shock3.png", 
       plot = multiplot(p7,p8,p9, cols = 1),
       width = 30, height = 20, units = "cm",
       device = "png")

### Obtain the exogenous shock###
Coef.noc
Coef.OLS
const

# check the order of variables is correct (stock shock) (employment)
lm(stock_price ~ lag(property_price)+lag(stock_price)+lag(employment)+
     +lag(property_price, n=2)+lag(stock_price, n=2)+lag(employment, n=2)+lag(property_price, n=3)+lag(stock_price, n=3)+lag(employment, n=3)+lag(property_price, n=4)+lag(stock_price, n=4)+lag(employment, n=4),
   data = data)

get_residual <- function(vlb_order = 1){
  X.Mat <- By %>% as_tibble() %>%
    mutate(V1 = lag(property_price),
           V2 = lag(stock_price),
           V3 = lag(employment),
           V4 = lag(property_price, n = 2),
           V5 = lag(stock_price, n = 2),
           V6 = lag(employment, n = 2),
           V7 = lag(property_price, n = 3),
           V8 = lag(stock_price, n = 3),
           V9 = lag(employment, n = 3),
           V10 = lag(property_price, n = 4),
           V11 = lag(stock_price, n = 4),
           V12 = lag(employment, n = 4),
           V13 = rep(1, nrow(By))
    ) %>%
    select(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13
    ) %>%
    drop_na() %>%
    as.matrix
  Y.Mat <- By[(VAR.P+1):nrow(By),vlb_order] %>% as.matrix
  Y.Mat_hat <- X.Mat %*% (Coef.OLS[vlb_order, ] %>% as.matrix)
  residual <- Y.Mat - Y.Mat_hat
  return(residual)
}

# resid1 <- get_residual(1)
for(i in 1:3){
  assign(paste0("resid_",i), get_residual(i))
}

resid.Mat <- cbind(resid_1, resid_2, resid_3)
error.Mat <- (solve(B0) %*% A0 %*% t(resid.Mat)) %>% t()
error.Mat[,1] %>% hist

# Extract the stock shock

error.Mat %>% as_tibble() %>%
  select(e1 = 1, e2 = 2, e3 = 3) %>%
  mutate(no = row_number()) %>%
  ggplot()+
  geom_line(aes(no, e1))

data.err <- data %>% 
  bind_cols(e1 = c(rep(NA, VAR.P), error.Mat[,1])) %>%
  mutate(e1_pos = fifelse(e1 > 0, e1, 0),
         e1_neg = fifelse(e1 < 0, e1, 0))
data.err %>% head(10)

data.err %>% ggplot()+
  geom_line(aes(year, e1))

# check the order of variables is correct (house shock) (employment)
lm(property_price ~ lag(property_price)+lag(stock_price)+lag(employment)+
     +lag(property_price, n=2)+lag(stock_price, n=2)+lag(employment, n=2)+lag(property_price, n=3)+lag(stock_price, n=3)+lag(employment, n=3)+lag(property_price, n=4)+lag(stock_price, n=4)+lag(employment, n=4),
   data = data)

get_residual <- function(vlb_order = 1){
  X.Mat <- By %>% as_tibble() %>%
    mutate(V1 = lag(property_price),
           V2 = lag(stock_price),
           V3 = lag(employment),
           V4 = lag(property_price, n = 2),
           V5 = lag(stock_price, n = 2),
           V6 = lag(employment, n = 2),
           V7 = lag(property_price, n = 3),
           V8 = lag(stock_price, n = 3),
           V9 = lag(employment, n = 3),
           V10 = lag(property_price, n = 4),
           V11 = lag(stock_price, n = 4),
           V12 = lag(employment, n = 4),
           V13 = rep(1, nrow(By))
    ) %>%
    select(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13
    ) %>%
    drop_na() %>%
    as.matrix
  Y.Mat <- By[(VAR.P+1):nrow(By),vlb_order] %>% as.matrix
  Y.Mat_hat <- X.Mat %*% (Coef.OLS[vlb_order, ] %>% as.matrix)
  residual <- Y.Mat - Y.Mat_hat
  return(residual)
}

# resid1 <- get_residual(1)
for(i in 1:3){
  assign(paste0("resid_",i), get_residual(i))
}

resid.Mat <- cbind(resid_1, resid_2, resid_3)
error.Mat <- (solve(B0) %*% A0 %*% t(resid.Mat)) %>% t()
error.Mat[,2] %>% hist

# Extract the house shock

error.Mat %>% as_tibble() %>%
  select(e1 = 1, e2 = 2, e3 = 3) %>%
  mutate(no = row_number()) %>%
  ggplot()+
  geom_line(aes(no, e2))

data.err2 <- data %>% 
  bind_cols(e2 = c(rep(NA, VAR.P), error.Mat[,2])) %>%
  mutate(e2_pos = fifelse(e2 > 0, e2, 0),
         e2_neg = fifelse(e2 < 0, e2, 0))
data.err2 %>% head(10)

data.err2 %>% ggplot()+
  geom_line(aes(year, e2))


data.err["e2"] = data.err2$e2
### Regression ###
reg1 <- lm(em_d ~ Period + e1 + e2 + interest_rate_d, data = data.err)
summary(reg1)

###衝擊減掉1994Q4###
data.err["e2"] = data.err["e2"] - data.err$e2[6]
data.err["e1"] = data.err["e1"] - data.err$e1[6]
reg1 <- lm(em_d ~ Period + e1 + e2 + interest_rate_d, data = data.err)
summary(reg1)




### DATA CPI ###
#file = "/Users/caiyawei/Desktop/Japan Data/paper_data_quart.xlsx"
file = "/Users/caiyawei/Desktop/Japan Data/SVAR_forR.xlsx"
data = read_excel(file)
data = data[-c(107),]
data_cpi = data

By <- data_cpi %>% 
  select(property_price, stock_price, cpi) %>% 
  as.matrix

### RAW DATA ###
dim(By)
kk <- dim(By)[2]
draw_raw <- function(vlb_name = ""){
  vlb_name <- sym(vlb_name)
  plot.return <- data_cpi %>%
    ggplot()+
    geom_line(aes(year, !!vlb_name))
  return(plot.return)
}

p1 <- draw_raw("property_price")
p2 <- draw_raw("stock_price")
p3 <- draw_raw("cpi")

multiplot(p1, p2, p3)

### Setting Model ###
#----- 模型設定 -----#
VAR.P = 4                       # 最大的落後項數
CONST = TRUE                    # 是否有常數項
Y     = VAR.Y(By, VAR.P)        # 設定 Y
X     = VAR.X(By, VAR.P)        # 設定 X

hrz=39 # the length of response

### Reduced Form VAR ###
(Coef.OLS    = VAR.OLS(Y, X, CONST)                  )
(Sigma.OLS   = VAR.Sigma.OLS(Y, X, Coef.OLS, CONST)  )
(Sigma.MLE   = VAR.Sigma.MLE(Y, X, Coef.OLS, CONST))

VAR.P = 4 # By AIC

### SVAR ###
### 6-variable model

Amat = diag(3)
# Identification Conditions

Amat[2,1]  = NA; 
Amat[3,1]  = NA; Amat[3,2]  = NA; 


Bmat = diag(3)
diag(Bmat) = NA

Amat;Bmat

### Estimate A^hat and B^hat ###
C.Prime <- chol(Sigma.OLS)
C <- t(C.Prime)
C

### Solving system of linear equations ###
B0 <- diag(diag(C), ncol = 3, nrow = 3)
B0

A0 <- B0 %*% solve(C)
A0

SVAR_AB_est <- list("A0.svar" = A0, "B0.svar" = B0)

### IRF Without Bootstrap CI ###
### IRF
SVAR_AB_IRF <- VAR.svarirf.AB(By, VAR.P, Amat, Bmat, h = hrz, CONST, SVAR_AB_est = SVAR_AB_est)


# 6*6個圖的time series
df_IRF_plot <- matrix(NA, hrz+1, kk^2) #%>% as.tibble() ## hrz+1
#dim(df_IRF_plot)
h <- 0 # h表示第幾期的IRF
for(period in SVAR_AB_IRF){
  k <- 0 # k表示把5*5的矩陣攤平到25個col的df時，要攤到第幾個columns上
  h <- h+1 # h表示第幾期的IRF
  for(j in 1:kk){
    for(i in 1:kk){
      k <- k+1 # k表示把5*5的矩陣攤平到25個col的df時，要攤到第幾個columns上
      df_IRF_plot[h,k] <- period[i,j]
    }
  }
}
df_IRF_plot <- df_IRF_plot %>% as_tibble()

kk*1:kk

#output entire table
IRF_TABLE <- df_IRF_plot[,kk*1:kk] %>% select(property_price = 1,
                                              stock_price = 2,
                                              cpi = 3)
write.table(IRF_TABLE, file = "/Users/caiyawei/Desktop/IRF_TABLE3_order2.csv", sep = ",", row.names = FALSE)

draw_IRF <- function(df = df_IRF_plot, V1 = 1){
  V1 <- paste0("V", V1) %>% sym()
  plot.return <- ggplot(df) + 
    geom_line(aes(x = 1:nrow(df), y = !!V1))
  return(plot.return)
}

for(i in 1:(kk^2)){
  assign(paste0("p",i), draw_IRF(df_IRF_plot, i))
}

multiplot(p1,p2,p3,p4,p5,p6,
          p7,p8,p9,cols = 3)

### IRF (Bootstrap C.I.) ###
lower = 0.025                                        # 控制成 95% CI
upper = 1-lower
kk = ncol(By)
ddY = VAR.ddY(By, VAR.P)
ddX = VAR.ddX(By, VAR.P)

# dim(ddY); dim(ddX)

T   = nrow(ddY)
T.total= nrow(By)
Ik  = diag(rep(1, kk))
# 16 coef if 4 variables; 55 coef if 5 variables
Coef = t(VAR.EbyE(ddY, ddX, CONST)$ddA)              # Step 1 估計模型
# residuals
U    = VAR.EbyE(ddY, ddX, CONST)$ddU
BSigma.u = VAR.ddSigma.OLS(ddY, ddX, CONST)
if(CONST == TRUE){
  const = Coef[, ncol(Coef)]
  Coef.noc= Coef[,-ncol(Coef)]                      # 刪掉 const
}else{
  const = matrix(0, kk, 1)
  Coef.noc = Coef
}

Theta.unit= VAR.Theta(Coef, h, BSigma.u, CONST)$unit # 估算 Theta.unit
Theta.std = VAR.Theta(Coef, h, BSigma.u, CONST)$std  # 估算 Theta.std

# dm.U <- U-mean(U)
dm.U <- U

N = 2000 #重抽次數
Theta.unit.sim = vector("list", N)
Theta.std.sim  = vector("list", N)

# check dimension
print("check dimensionality")
dim(ddX); dim(Coef.noc); dim(dm.U)


# 存N次重抽的IRF
df_IRF.sim <- array(NA, c(hrz+1,kk^2,N)) #dimensions are: Time Period, Number of shock interacts with variables, page (number of Bootstrap resamplings)
counter <- 1
while(TRUE){
  
  #cat("Now, there are ", counter-1, " sets of resamples.\n")
  Y.sim = matrix(0, nrow = T.total, ncol = kk)          # Y.sim = 0 #pseudo time series
  Y.sim[c(1:VAR.P),] = By[c(1:VAR.P), ] #initial values
  
  boot.number = sample(c(1:T), replace = TRUE)      # Step 3 取出放回
  U.sim = dm.U[boot.number,]
  
  # predicted values given the above initial values
  last.y= c(t(By[VAR.P:1,]))
  for(ii in 1:T){
    last.y = last.y[1:(kk*VAR.P)]
    Y.sim[ii+VAR.P, ] = Coef.noc %*% last.y + const + U.sim[ii,]      # Step 4 模擬資料
    last.y = c(Y.sim[ii+VAR.P,], last.y)
  }
  
  #   Y.sim[-c(1:VAR.P),] <- matrix(const, nrow = T.total-VAR.P, ncol = kk, byrow = T) + ddX %*% t(Coef.noc) + U.sim
  
  
  #`Y.sim` is the pseudo time series
  # Step 5 重新估算SVAR
  
  ### SVAR.sim Start ###
  
  Y_pseudo     = VAR.Y(Y.sim, VAR.P)        # 設定 Y
  X_pseudo     = VAR.X(Y.sim, VAR.P)        # 設定 X
  Coef.OLS_pseudo    = VAR.OLS(Y_pseudo, X_pseudo, CONST)
  Sigma.OLS_pseudo   = VAR.Sigma.OLS(Y_pseudo, X_pseudo, Coef.OLS_pseudo, CONST)
  C.Prime_pseudo <- chol(Sigma.OLS_pseudo)
  C_pseudo <- t(C.Prime_pseudo)
  B0_pseudo <- diag(diag(C_pseudo), ncol = kk, nrow = kk)
  A0_pseudo <- B0_pseudo %*% solve(C_pseudo)
  SVAR_AB_est.sim <- list("A0.svar" = A0_pseudo, "B0.svar" = B0_pseudo)
  SVAR_AB_IRF.sim <- VAR.svarirf.AB(Y.sim, VAR.P, Amat, Bmat, h = hrz, CONST, SVAR_AB_est = SVAR_AB_est.sim)
  
  # 5*5個圖的time series
  df_IRF_plot.sim <- matrix(NA, hrz+1, kk^2) #%>% as.tibble()
  # df_IRF.sim <- array(1:(120*25*N), c(120,25,N))
  # df_IRF.sim[2,1,1] # slicing
  
  h <- 0 # h表示第幾期的IRF
  for(period in SVAR_AB_IRF.sim){
    k <- 0 # k表示把5*5的矩陣攤平到25個col的df時，要攤到第幾個columns上
    h <- h+1 # h表示第幾期的IRF
    for(j in 1:kk){
      for(i in 1:kk){
        k <- k+1 # k表示把5*5的矩陣攤平到25個col的df時，要攤到第幾個columns上
        df_IRF_plot.sim[h,k] <- period[i,j]
      }
    }
  }
  # 把這一次重抽得到的IRF append進`df_IRF.sim`中
  df_IRF.sim[,,counter] <- df_IRF_plot.sim
  ### SVAR.sim Ends ###
  if(counter>=N){
    break
  }
  counter <- counter+1
}
# Save
saveRDS(df_IRF.sim, file = "df_IRF.sim.rds")

df_IRF.sim <- read_rds("df_IRF.sim.rds")

# 看某一頁
head(df_IRF.sim[,,1000])
print(sum(is.na(df_IRF.sim)))

# 畫IRF & Bootstrap C.I.
df_IRF_plot.BS.L <- matrix(NA, nrow = hrz+1, ncol = kk^2)
df_IRF_plot.BS.U <- matrix(NA, nrow = hrz+1, ncol = kk^2)
df_IRF_plot.BS.Median <- matrix(NA, nrow = hrz+1, ncol = kk^2)
df_IRF_plot.BS.Mean <- matrix(NA, nrow = hrz+1, ncol = kk^2)
for(col in 1:(kk^2)){
  for(row in 1:(hrz+1) ){
    df_IRF_plot.BS.L[row,col] <- quantile(df_IRF.sim[row,col,], probs = 0.025)
    df_IRF_plot.BS.U[row,col] <- quantile(df_IRF.sim[row,col,], probs = 0.975)
    df_IRF_plot.BS.Median[row,col] <- quantile(df_IRF.sim[row,col,], probs = 0.5)
    df_IRF_plot.BS.Mean[row,col] <- mean(df_IRF.sim[row,col,])
  }
}

df_IRF_plot.BS.L <- df_IRF_plot.BS.L %>% as_tibble()
df_IRF_plot.BS.U <- df_IRF_plot.BS.U %>% as_tibble()
df_IRF_plot.BS.Median <- df_IRF_plot.BS.Median %>% as_tibble()
df_IRF_plot.BS.Mean <- df_IRF_plot.BS.Mean %>% as_tibble()

ind <- 0
for(i in 1:kk){
  for(j in 1:kk){
    ind <- ind+1
    nam <- paste("shock", j, "y", i, sep = '')
    assign(nam, bind_cols(df_IRF_plot.BS.L[ind], df_IRF_plot.BS.U[ind],
                          df_IRF_plot.BS.Median[ind], df_IRF_plot.BS.Mean[ind],
                          df_IRF_plot[ind]))
    # 改名
    evalStr <- paste0("colnames(", nam, ") <- c('Lower', 'Upper', 'Median', 'Mean', 'Actual')")
    eval(parse(text=evalStr))
    # 圖層
    evalStr <- paste0("p", ind, " <- ", "ggplot(",nam,") +geom_hline(yintercept=0, color = 'grey')+ geom_line(aes(x = 1:nrow(", nam, "), y = Lower), linetype = 'dashed', col='red')+geom_line(aes(x = 1:nrow(", nam, "), y = Upper), linetype = 'dashed', col='red')+geom_line(aes(x = 1:nrow(", nam, "), y = Median), col = 'Blue')")
    eval(parse(text=evalStr))
  }
}  

Text_Size_Theme = theme(
  axis.title.x = element_text(size = 12),
  axis.text.x = element_text(size = 12),
  axis.title.y = element_text(size = 12),
  axis.title = element_text(size = 12),
  plot.title = element_text(size=12))

add_layout <- function(p = p1,
                       response_of = "this variable", react_to = "the shock"){
  title_text <- paste("Response of", response_of, "to", react_to, sep = ' ')
  plot.return <- p + labs(x = "Period", y = "", title = title_text) + Text_Size_Theme
  return(plot.return)
}

vlb_list = c("House Price",
             "Stock Price",
             "CPI"
)
p1 <- add_layout(p1, vlb_list[1], "Property Price shock")
p2 <- add_layout(p2, vlb_list[2], "Property Price shock")
p3 <- add_layout(p3, vlb_list[3], "Property Price shock")

p4 <- add_layout(p4, vlb_list[1], "Stock shock")
p5 <- add_layout(p5, vlb_list[2], "Stock shock")
p6 <- add_layout(p6, vlb_list[3], "Stock shock")

p7 <- add_layout(p7, vlb_list[1], "CPI shock")
p8 <- add_layout(p8, vlb_list[2], "CPI shock")
p9 <- add_layout(p9, vlb_list[3], "CPI shock")

multiplot(p1,p2,p3,p4,p5,p6,
          p7,p8,p9,cols = 3)

multiplot(p1,p2,p3, cols = 1)
ggsave(filename = "/Users/caiyawei/Desktop/VAR/IRF_shock_cpi1.png", 
       plot = multiplot(p1,p2,p3, cols = 1),
       width = 30, height = 20, units = "cm",
       device = "png")

multiplot(p4,p5,p6, cols = 1)
ggsave(filename = "/Users/caiyawei/Desktop/VAR/IRF_shock_cpi2.png", 
       plot = multiplot(p4,p5,p6, cols = 1),
       width = 30, height = 20, units = "cm",
       device = "png")

multiplot(p7,p8,p9, cols = 1)
ggsave(filename = "/Users/caiyawei/Desktop/VAR/IRF_shock_cpi3.png", 
       plot = multiplot(p7,p8,p9, cols = 1),
       width = 30, height = 20, units = "cm",
       device = "png")

### Obtain the exogenous shock###
Coef.noc
Coef.OLS
const

# check the order of variables is correct (stock shock) (cpi)
lm(stock_price ~ lag(property_price)+lag(stock_price)+lag(cpi)+
     +lag(property_price, n=2)+lag(stock_price, n=2)+lag(cpi, n=2)+lag(property_price, n=3)+lag(stock_price, n=3)+lag(cpi, n=3)+lag(property_price, n=4)+lag(stock_price, n=4)+lag(cpi, n=4),
   data = data_cpi)

get_residual <- function(vlb_order = 1){
  X.Mat <- By %>% as_tibble() %>%
    mutate(V1 = lag(property_price),
           V2 = lag(stock_price),
           V3 = lag(cpi),
           V4 = lag(property_price, n = 2),
           V5 = lag(stock_price, n = 2),
           V6 = lag(cpi, n = 2),
           V7 = lag(property_price, n = 3),
           V8 = lag(stock_price, n = 3),
           V9 = lag(cpi, n = 3),
           V10 = lag(property_price, n = 4),
           V11 = lag(stock_price, n = 4),
           V12 = lag(cpi, n = 4),
           V13 = rep(1, nrow(By))
    ) %>%
    select(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13
    ) %>%
    drop_na() %>%
    as.matrix
  Y.Mat <- By[(VAR.P+1):nrow(By),vlb_order] %>% as.matrix
  Y.Mat_hat <- X.Mat %*% (Coef.OLS[vlb_order, ] %>% as.matrix)
  residual <- Y.Mat - Y.Mat_hat
  return(residual)
}

# resid1 <- get_residual(1)
for(i in 1:3){
  assign(paste0("resid_",i), get_residual(i))
}

resid.Mat <- cbind(resid_1, resid_2, resid_3)
error.Mat <- (solve(B0) %*% A0 %*% t(resid.Mat)) %>% t()
error.Mat[,1] %>% hist

# Extract the stock shock

error.Mat %>% as_tibble() %>%
  select(e1 = 1, e2 = 2, e3 = 3) %>%
  mutate(no = row_number()) %>%
  ggplot()+
  geom_line(aes(no, e1))

data.err <- data_cpi %>% 
  bind_cols(e1 = c(rep(NA, VAR.P), error.Mat[,1])) %>%
  mutate(e1_pos = fifelse(e1 > 0, e1, 0),
         e1_neg = fifelse(e1 < 0, e1, 0))
data.err %>% head(10)

data.err %>% ggplot()+
  geom_line(aes(year, e1))

# check the order of variables is correct (house shock) (cpi)
lm(property_price ~ lag(property_price)+lag(stock_price)+lag(cpi)+
     +lag(property_price, n=2)+lag(stock_price, n=2)+lag(cpi, n=2)+lag(property_price, n=3)+lag(stock_price, n=3)+lag(cpi, n=3)+lag(property_price, n=4)+lag(stock_price, n=4)+lag(cpi, n=4),
   data = data_cpi)

get_residual <- function(vlb_order = 1){
  X.Mat <- By %>% as_tibble() %>%
    mutate(V1 = lag(property_price),
           V2 = lag(stock_price),
           V3 = lag(cpi),
           V4 = lag(property_price, n = 2),
           V5 = lag(stock_price, n = 2),
           V6 = lag(cpi, n = 2),
           V7 = lag(property_price, n = 3),
           V8 = lag(stock_price, n = 3),
           V9 = lag(cpi, n = 3),
           V10 = lag(property_price, n = 4),
           V11 = lag(stock_price, n = 4),
           V12 = lag(cpi, n = 4),
           V13 = rep(1, nrow(By))
    ) %>%
    select(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13
    ) %>%
    drop_na() %>%
    as.matrix
  Y.Mat <- By[(VAR.P+1):nrow(By),vlb_order] %>% as.matrix
  Y.Mat_hat <- X.Mat %*% (Coef.OLS[vlb_order, ] %>% as.matrix)
  residual <- Y.Mat - Y.Mat_hat
  return(residual)
}

# resid1 <- get_residual(1)
for(i in 1:3){
  assign(paste0("resid_",i), get_residual(i))
}

resid.Mat <- cbind(resid_1, resid_2, resid_3)
error.Mat <- (solve(B0) %*% A0 %*% t(resid.Mat)) %>% t()
error.Mat[,2] %>% hist

# Extract the house shock

error.Mat %>% as_tibble() %>%
  select(e1 = 1, e2 = 2, e3 = 3) %>%
  mutate(no = row_number()) %>%
  ggplot()+
  geom_line(aes(no, e2))

data.err2 <- data_cpi %>% 
  bind_cols(e2 = c(rep(NA, VAR.P), error.Mat[,2])) %>%
  mutate(e2_pos = fifelse(e2 > 0, e2, 0),
         e2_neg = fifelse(e2 < 0, e2, 0))
data.err2 %>% head(10)

data.err2 %>% ggplot()+
  geom_line(aes(year, e2))


data.err["e2"] = data.err2$e2
### Regression ###
reg3 <- lm(cpi_d ~ Period + e1 + e2 + interest_rate_d, data = data.err)
summary(reg3)
###衝擊減掉1991Q4###
data.err["e2"] = data.err["e2"] - data.err$e2[6]
data.err["e1"] = data.err["e1"] - data.err$e1[6]
reg3 <- lm(cpi_d ~ Period + e1 + e2 + interest_rate_d, data = data.err)
summary(reg3)







### DATA GDP ###
data = read_excel(file)
data = data[-c(107),]
data_gdp = data

By <- data_gdp %>% 
  select(property_price, stock_price, gdp) %>% 
  as.matrix

### RAW DATA ###
dim(By)
kk <- dim(By)[2]
draw_raw <- function(vlb_name = ""){
  vlb_name <- sym(vlb_name)
  plot.return <- data_gdp %>%
    ggplot()+
    geom_line(aes(year, !!vlb_name))
  return(plot.return)
}

p1 <- draw_raw("property_price")
p2 <- draw_raw("stock_price")
p3 <- draw_raw("gdp")

multiplot(p1, p2, p3)

### Setting Model ###
#----- 模型設定 -----#
VAR.P = 4                       # 最大的落後項數
CONST = TRUE                    # 是否有常數項
Y     = VAR.Y(By, VAR.P)        # 設定 Y
X     = VAR.X(By, VAR.P)        # 設定 X

hrz=39 # the length of response

### Reduced Form VAR ###
(Coef.OLS    = VAR.OLS(Y, X, CONST)                  )
(Sigma.OLS   = VAR.Sigma.OLS(Y, X, Coef.OLS, CONST)  )
(Sigma.MLE   = VAR.Sigma.MLE(Y, X, Coef.OLS, CONST))

VAR.P = 4 # By AIC

### SVAR ###
### 6-variable model

Amat = diag(3)
# Identification Conditions

Amat[2,1]  = NA; 
Amat[3,1]  = NA; Amat[3,2]  = NA; 


Bmat = diag(3)
diag(Bmat) = NA

Amat;Bmat

### Estimate A^hat and B^hat ###
C.Prime <- chol(Sigma.OLS)
C <- t(C.Prime)
C

### Solving system of linear equations ###
B0 <- diag(diag(C), ncol = 3, nrow = 3)
B0

A0 <- B0 %*% solve(C)
A0

SVAR_AB_est <- list("A0.svar" = A0, "B0.svar" = B0)

### IRF Without Bootstrap CI ###
### IRF
SVAR_AB_IRF <- VAR.svarirf.AB(By, VAR.P, Amat, Bmat, h = hrz, CONST, SVAR_AB_est = SVAR_AB_est)


# 6*6個圖的time series
df_IRF_plot <- matrix(NA, hrz+1, kk^2) #%>% as.tibble() ## hrz+1
#dim(df_IRF_plot)
h <- 0 # h表示第幾期的IRF
for(period in SVAR_AB_IRF){
  k <- 0 # k表示把5*5的矩陣攤平到25個col的df時，要攤到第幾個columns上
  h <- h+1 # h表示第幾期的IRF
  for(j in 1:kk){
    for(i in 1:kk){
      k <- k+1 # k表示把5*5的矩陣攤平到25個col的df時，要攤到第幾個columns上
      df_IRF_plot[h,k] <- period[i,j]
    }
  }
}
df_IRF_plot <- df_IRF_plot %>% as_tibble()

kk*1:kk

#output entire table
IRF_TABLE <- df_IRF_plot[,kk*1:kk] %>% select(property_price = 1,
                                              stock_price = 2,
                                              gdp = 3)
write.table(IRF_TABLE, file = "/Users/caiyawei/Desktop/IRF_TABLE2.csv", sep = ",", row.names = FALSE)

draw_IRF <- function(df = df_IRF_plot, V1 = 1){
  V1 <- paste0("V", V1) %>% sym()
  plot.return <- ggplot(df) + 
    geom_line(aes(x = 1:nrow(df), y = !!V1))
  return(plot.return)
}

for(i in 1:(kk^2)){
  assign(paste0("p",i), draw_IRF(df_IRF_plot, i))
}

multiplot(p1,p2,p3,p4,p5,p6,
          p7,p8,p9,cols = 3)

### IRF (Bootstrap C.I.) ###
lower = 0.025                                        # 控制成 95% CI
upper = 1-lower
kk = ncol(By)
ddY = VAR.ddY(By, VAR.P)
ddX = VAR.ddX(By, VAR.P)

# dim(ddY); dim(ddX)

T   = nrow(ddY)
T.total= nrow(By)
Ik  = diag(rep(1, kk))
# 16 coef if 4 variables; 55 coef if 5 variables
Coef = t(VAR.EbyE(ddY, ddX, CONST)$ddA)              # Step 1 估計模型
# residuals
U    = VAR.EbyE(ddY, ddX, CONST)$ddU
BSigma.u = VAR.ddSigma.OLS(ddY, ddX, CONST)
if(CONST == TRUE){
  const = Coef[, ncol(Coef)]
  Coef.noc= Coef[,-ncol(Coef)]                      # 刪掉 const
}else{
  const = matrix(0, kk, 1)
  Coef.noc = Coef
}

Theta.unit= VAR.Theta(Coef, h, BSigma.u, CONST)$unit # 估算 Theta.unit
Theta.std = VAR.Theta(Coef, h, BSigma.u, CONST)$std  # 估算 Theta.std

# dm.U <- U-mean(U)
dm.U <- U

N = 2000 #重抽次數
Theta.unit.sim = vector("list", N)
Theta.std.sim  = vector("list", N)

# check dimension
print("check dimensionality")
dim(ddX); dim(Coef.noc); dim(dm.U)


# 存N次重抽的IRF
df_IRF.sim <- array(NA, c(hrz+1,kk^2,N)) #dimensions are: Time Period, Number of shock interacts with variables, page (number of Bootstrap resamplings)
counter <- 1
while(TRUE){
  
  #cat("Now, there are ", counter-1, " sets of resamples.\n")
  Y.sim = matrix(0, nrow = T.total, ncol = kk)          # Y.sim = 0 #pseudo time series
  Y.sim[c(1:VAR.P),] = By[c(1:VAR.P), ] #initial values
  
  boot.number = sample(c(1:T), replace = TRUE)      # Step 3 取出放回
  U.sim = dm.U[boot.number,]
  
  # predicted values given the above initial values
  last.y= c(t(By[VAR.P:1,]))
  for(ii in 1:T){
    last.y = last.y[1:(kk*VAR.P)]
    Y.sim[ii+VAR.P, ] = Coef.noc %*% last.y + const + U.sim[ii,]      # Step 4 模擬資料
    last.y = c(Y.sim[ii+VAR.P,], last.y)
  }
  
  #   Y.sim[-c(1:VAR.P),] <- matrix(const, nrow = T.total-VAR.P, ncol = kk, byrow = T) + ddX %*% t(Coef.noc) + U.sim
  
  
  #`Y.sim` is the pseudo time series
  # Step 5 重新估算SVAR
  
  ### SVAR.sim Start ###
  
  Y_pseudo     = VAR.Y(Y.sim, VAR.P)        # 設定 Y
  X_pseudo     = VAR.X(Y.sim, VAR.P)        # 設定 X
  Coef.OLS_pseudo    = VAR.OLS(Y_pseudo, X_pseudo, CONST)
  Sigma.OLS_pseudo   = VAR.Sigma.OLS(Y_pseudo, X_pseudo, Coef.OLS_pseudo, CONST)
  C.Prime_pseudo <- chol(Sigma.OLS_pseudo)
  C_pseudo <- t(C.Prime_pseudo)
  B0_pseudo <- diag(diag(C_pseudo), ncol = kk, nrow = kk)
  A0_pseudo <- B0_pseudo %*% solve(C_pseudo)
  SVAR_AB_est.sim <- list("A0.svar" = A0_pseudo, "B0.svar" = B0_pseudo)
  SVAR_AB_IRF.sim <- VAR.svarirf.AB(Y.sim, VAR.P, Amat, Bmat, h = hrz, CONST, SVAR_AB_est = SVAR_AB_est.sim)
  
  # 5*5個圖的time series
  df_IRF_plot.sim <- matrix(NA, hrz+1, kk^2) #%>% as.tibble()
  # df_IRF.sim <- array(1:(120*25*N), c(120,25,N))
  # df_IRF.sim[2,1,1] # slicing
  
  h <- 0 # h表示第幾期的IRF
  for(period in SVAR_AB_IRF.sim){
    k <- 0 # k表示把5*5的矩陣攤平到25個col的df時，要攤到第幾個columns上
    h <- h+1 # h表示第幾期的IRF
    for(j in 1:kk){
      for(i in 1:kk){
        k <- k+1 # k表示把5*5的矩陣攤平到25個col的df時，要攤到第幾個columns上
        df_IRF_plot.sim[h,k] <- period[i,j]
      }
    }
  }
  # 把這一次重抽得到的IRF append進`df_IRF.sim`中
  df_IRF.sim[,,counter] <- df_IRF_plot.sim
  ### SVAR.sim Ends ###
  if(counter>=N){
    break
  }
  counter <- counter+1
}
# Save
saveRDS(df_IRF.sim, file = "df_IRF.sim.rds")

df_IRF.sim <- read_rds("df_IRF.sim.rds")

# 看某一頁
head(df_IRF.sim[,,1000])
print(sum(is.na(df_IRF.sim)))

# 畫IRF & Bootstrap C.I.
df_IRF_plot.BS.L <- matrix(NA, nrow = hrz+1, ncol = kk^2)
df_IRF_plot.BS.U <- matrix(NA, nrow = hrz+1, ncol = kk^2)
df_IRF_plot.BS.Median <- matrix(NA, nrow = hrz+1, ncol = kk^2)
df_IRF_plot.BS.Mean <- matrix(NA, nrow = hrz+1, ncol = kk^2)
for(col in 1:(kk^2)){
  for(row in 1:(hrz+1) ){
    df_IRF_plot.BS.L[row,col] <- quantile(df_IRF.sim[row,col,], probs = 0.025)
    df_IRF_plot.BS.U[row,col] <- quantile(df_IRF.sim[row,col,], probs = 0.975)
    df_IRF_plot.BS.Median[row,col] <- quantile(df_IRF.sim[row,col,], probs = 0.5)
    df_IRF_plot.BS.Mean[row,col] <- mean(df_IRF.sim[row,col,])
  }
}

df_IRF_plot.BS.L <- df_IRF_plot.BS.L %>% as_tibble()
df_IRF_plot.BS.U <- df_IRF_plot.BS.U %>% as_tibble()
df_IRF_plot.BS.Median <- df_IRF_plot.BS.Median %>% as_tibble()
df_IRF_plot.BS.Mean <- df_IRF_plot.BS.Mean %>% as_tibble()

ind <- 0
for(i in 1:kk){
  for(j in 1:kk){
    ind <- ind+1
    nam <- paste("shock", j, "y", i, sep = '')
    assign(nam, bind_cols(df_IRF_plot.BS.L[ind], df_IRF_plot.BS.U[ind],
                          df_IRF_plot.BS.Median[ind], df_IRF_plot.BS.Mean[ind],
                          df_IRF_plot[ind]))
    # 改名
    evalStr <- paste0("colnames(", nam, ") <- c('Lower', 'Upper', 'Median', 'Mean', 'Actual')")
    eval(parse(text=evalStr))
    # 圖層
    evalStr <- paste0("p", ind, " <- ", "ggplot(",nam,") +geom_hline(yintercept=0, color = 'grey')+ geom_line(aes(x = 1:nrow(", nam, "), y = Lower), linetype = 'dashed', col='red')+geom_line(aes(x = 1:nrow(", nam, "), y = Upper), linetype = 'dashed', col='red')+geom_line(aes(x = 1:nrow(", nam, "), y = Median), col = 'Blue')")
    eval(parse(text=evalStr))
  }
}  

Text_Size_Theme = theme(
  axis.title.x = element_text(size = 12),
  axis.text.x = element_text(size = 12),
  axis.title.y = element_text(size = 12),
  axis.title = element_text(size = 12),
  plot.title = element_text(size=12))

add_layout <- function(p = p1,
                       response_of = "this variable", react_to = "the shock"){
  title_text <- paste("Response of", response_of, "to", react_to, sep = ' ')
  plot.return <- p + labs(x = "Period", y = "", title = title_text) + Text_Size_Theme
  return(plot.return)
}

vlb_list = c("House Price",
             "Stock Price",
             "GDP"
)
p1 <- add_layout(p1, vlb_list[1], "Property Price shock")
p2 <- add_layout(p2, vlb_list[2], "Property Price shock")
p3 <- add_layout(p3, vlb_list[3], "Property Price shock")

p4 <- add_layout(p4, vlb_list[1], "Stock shock")
p5 <- add_layout(p5, vlb_list[2], "Stock shock")
p6 <- add_layout(p6, vlb_list[3], "Stock shock")

p7 <- add_layout(p7, vlb_list[1], "GDP shock")
p8 <- add_layout(p8, vlb_list[2], "GDP shock")
p9 <- add_layout(p9, vlb_list[3], "GDP shock")

multiplot(p1,p2,p3,p4,p5,p6,
          p7,p8,p9,cols = 3)

multiplot(p1,p2,p3, cols = 1)
ggsave(filename = "/Users/caiyawei/Desktop/VAR/IRF_shock_gdp1.png", 
       plot = multiplot(p1,p2,p3, cols = 1),
       width = 30, height = 20, units = "cm",
       device = "png")

multiplot(p4,p5,p6, cols = 1)
ggsave(filename = "/Users/caiyawei/Desktop/VAR/IRF_shock_gdp2.png", 
       plot = multiplot(p4,p5,p6, cols = 1),
       width = 30, height = 20, units = "cm",
       device = "png")

multiplot(p3,p6,p9, cols = 1)
ggsave(filename = "/Users/caiyawei/Desktop/VAR/IRF_shock_gdp33.png", 
       plot = multiplot(p7,p8,p9, cols = 1),
       width = 30, height = 20, units = "cm",
       device = "png")

### Obtain the exogenous shock###
Coef.noc
Coef.OLS
const

# check the order of variables is correct (stock shock) (gdp)
lm(stock_price ~ lag(property_price)+lag(stock_price)+lag(gdp)+
     +lag(property_price, n=2)+lag(stock_price, n=2)+lag(gdp, n=2)+lag(property_price, n=3)+lag(stock_price, n=3)+lag(gdp, n=3)+lag(property_price, n=4)+lag(stock_price, n=4)+lag(gdp, n=4),
   data = data_gdp)

get_residual <- function(vlb_order = 1){
  X.Mat <- By %>% as_tibble() %>%
    mutate(V1 = lag(property_price),
           V2 = lag(stock_price),
           V3 = lag(gdp),
           V4 = lag(property_price, n = 2),
           V5 = lag(stock_price, n = 2),
           V6 = lag(gdp, n = 2),
           V7 = lag(property_price, n = 3),
           V8 = lag(stock_price, n = 3),
           V9 = lag(gdp, n = 3),
           V10 = lag(property_price, n = 4),
           V11 = lag(stock_price, n = 4),
           V12 = lag(gdp, n = 4),
           V13 = rep(1, nrow(By))
    ) %>%
    select(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13
    ) %>%
    drop_na() %>%
    as.matrix
  Y.Mat <- By[(VAR.P+1):nrow(By),vlb_order] %>% as.matrix
  Y.Mat_hat <- X.Mat %*% (Coef.OLS[vlb_order, ] %>% as.matrix)
  residual <- Y.Mat - Y.Mat_hat
  return(residual)
}

# resid1 <- get_residual(1)
for(i in 1:3){
  assign(paste0("resid_",i), get_residual(i))
}

resid.Mat <- cbind(resid_1, resid_2, resid_3)
error.Mat <- (solve(B0) %*% A0 %*% t(resid.Mat)) %>% t()
error.Mat[,1] %>% hist

# Extract the stock shock

error.Mat %>% as_tibble() %>%
  select(e1 = 1, e2 = 2, e3 = 3) %>%
  mutate(no = row_number()) %>%
  ggplot()+
  geom_line(aes(no, e1))

data.err <- data_gdp %>% 
  bind_cols(e1 = c(rep(NA, VAR.P), error.Mat[,1])) %>%
  mutate(e1_pos = fifelse(e1 > 0, e1, 0),
         e1_neg = fifelse(e1 < 0, e1, 0))
data.err %>% head(10)

data.err %>% ggplot()+
  geom_line(aes(year, e1))

# check the order of variables is correct (house shock) (gdp)
lm(property_price ~ lag(property_price)+lag(stock_price)+lag(gdp)+
     +lag(property_price, n=2)+lag(stock_price, n=2)+lag(gdp, n=2)+lag(property_price, n=3)+lag(stock_price, n=3)+lag(gdp, n=3)+lag(property_price, n=4)+lag(stock_price, n=4)+lag(gdp, n=4),
   data = data_gdp)

get_residual <- function(vlb_order = 1){
  X.Mat <- By %>% as_tibble() %>%
    mutate(V1 = lag(property_price),
           V2 = lag(stock_price),
           V3 = lag(gdp),
           V4 = lag(property_price, n = 2),
           V5 = lag(stock_price, n = 2),
           V6 = lag(gdp, n = 2),
           V7 = lag(property_price, n = 3),
           V8 = lag(stock_price, n = 3),
           V9 = lag(gdp, n = 3),
           V10 = lag(property_price, n = 4),
           V11 = lag(stock_price, n = 4),
           V12 = lag(gdp, n = 4),
           V13 = rep(1, nrow(By))
    ) %>%
    select(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13
    ) %>%
    drop_na() %>%
    as.matrix
  Y.Mat <- By[(VAR.P+1):nrow(By),vlb_order] %>% as.matrix
  Y.Mat_hat <- X.Mat %*% (Coef.OLS[vlb_order, ] %>% as.matrix)
  residual <- Y.Mat - Y.Mat_hat
  return(residual)
}

# resid1 <- get_residual(1)
for(i in 1:3){
  assign(paste0("resid_",i), get_residual(i))
}

resid.Mat <- cbind(resid_1, resid_2, resid_3)
error.Mat <- (solve(B0) %*% A0 %*% t(resid.Mat)) %>% t()
error.Mat[,2] %>% hist

# Extract the house shock

error.Mat %>% as_tibble() %>%
  select(e1 = 1, e2 = 2, e3 = 3) %>%
  mutate(no = row_number()) %>%
  ggplot()+
  geom_line(aes(no, e2))

data.err2 <- data_gdp %>% 
  bind_cols(e2 = c(rep(NA, VAR.P), error.Mat[,2])) %>%
  mutate(e2_pos = fifelse(e2 > 0, e2, 0),
         e2_neg = fifelse(e2 < 0, e2, 0))
data.err2 %>% head(10)

data.err2 %>% ggplot()+
  geom_line(aes(year...1, e2))


data.err["e2"] = data.err2$e2
### Regression ###
#reg2 <- lm(gdp_d ~ Period + e1 + e2 + interest_rate_d, data = data.err)
#summary(reg2)


### 衝擊減掉1991Q4 ###
data.err["e2"] = data.err["e2"] - data.err$e2[6]
data.err["e1"] = data.err["e1"] - data.err$e1[6]
reg2 <- lm(gdp_d ~ Period + e1 + e2 + interest_rate_d, data = data.err)
summary(reg2)

texreg(list(reg1, reg2, reg3), caption = "Three SVAR models.", digits = 5)




##### VAR result in Order1 EM #####
file = "/Users/caiyawei/Desktop/Japan Data/SVAR_forR.xlsx"
data = read_excel(file)
data = data[-c(107),]

###Order1 Stock -> House -> Macro###
### Employment ###
#data$date <- as.Date(data$date)

By <- data %>% 
  select(stock_price, property_price, employment) %>% 
  as.matrix

### RAW DATA ###
dim(By)
kk <- dim(By)[2]
draw_raw <- function(vlb_name = ""){
  vlb_name <- sym(vlb_name)
  plot.return <- data %>%
    ggplot()+
    geom_line(aes(year, !!vlb_name))
  return(plot.return)
}

p1 <- draw_raw("stock_price")
p2 <- draw_raw("property_price")
p3 <- draw_raw("employment")

multiplot(p1, p2, p3)

### Setting Model ###
#----- 模型設定 -----#
VAR.P = 4                       # 最大的落後項數
CONST = TRUE                    # 是否有常數項
Y     = VAR.Y(By, VAR.P)        # 設定 Y
X     = VAR.X(By, VAR.P)        # 設定 X

hrz=39 # the length of response

### Reduced Form VAR ###
(Coef.OLS    = VAR.OLS(Y, X, CONST)                  )
(Sigma.OLS   = VAR.Sigma.OLS(Y, X, Coef.OLS, CONST)  )
(Sigma.MLE   = VAR.Sigma.MLE(Y, X, Coef.OLS, CONST))

VAR.P = 4 # By AIC

### SVAR ###
### 6-variable model

Amat = diag(3)
# Identification Conditions

Bmat = diag(3)


Amat;Bmat

### Estimate A^hat and B^hat ###
#C.Prime <- chol(Sigma.OLS)
#C <- t(C.Prime)
#C

### Solving system of linear equations ###
B0 <-Bmat


A0 <-Amat


SVAR_AB_est <- list("A0.svar" = A0, "B0.svar" = B0)

### IRF Without Bootstrap CI ###
### IRF
SVAR_AB_IRF <- VAR.svarirf.AB(By, VAR.P, Amat, Bmat, h = hrz, CONST, SVAR_AB_est = SVAR_AB_est)


# 6*6個圖的time series
df_IRF_plot <- matrix(NA, hrz+1, kk^2) #%>% as.tibble() ## hrz+1
#dim(df_IRF_plot)
h <- 0 # h表示第幾期的IRF
for(period in SVAR_AB_IRF){
  k <- 0 # k表示把5*5的矩陣攤平到25個col的df時，要攤到第幾個columns上
  h <- h+1 # h表示第幾期的IRF
  for(j in 1:kk){
    for(i in 1:kk){
      k <- k+1 # k表示把5*5的矩陣攤平到25個col的df時，要攤到第幾個columns上
      df_IRF_plot[h,k] <- period[i,j]
    }
  }
}
df_IRF_plot <- df_IRF_plot %>% as_tibble()

kk*1:kk



#output entire table
IRF_TABLE <- df_IRF_plot[,kk*1:kk] %>% select(stock_price = 1,
                                              property_price = 2,
                                              employment = 3)
write.table(IRF_TABLE, file = "/Users/caiyawei/Desktop/IRF_TABLE.csv", sep = ",", row.names = FALSE)

draw_IRF <- function(df = df_IRF_plot, V1 = 1){
  V1 <- paste0("V", V1) %>% sym()
  plot.return <- ggplot(df) + 
    geom_line(aes(x = 1:nrow(df), y = !!V1))
  return(plot.return)
}

for(i in 1:(kk^2)){
  assign(paste0("p",i), draw_IRF(df_IRF_plot, i))
}

multiplot(p1,p2,p3,p4,p5,p6,
          p7,p8,p9,cols = 3)

### IRF (Bootstrap C.I.) ###
lower = 0.025                                        # 控制成 95% CI
upper = 1-lower
kk = ncol(By)
ddY = VAR.ddY(By, VAR.P)
ddX = VAR.ddX(By, VAR.P)

# dim(ddY); dim(ddX)

T   = nrow(ddY)
T.total= nrow(By)
Ik  = diag(rep(1, kk))
# 16 coef if 4 variables; 55 coef if 5 variables
Coef = t(VAR.EbyE(ddY, ddX, CONST)$ddA)              # Step 1 估計模型
# residuals
U    = VAR.EbyE(ddY, ddX, CONST)$ddU
BSigma.u = VAR.ddSigma.OLS(ddY, ddX, CONST)
if(CONST == TRUE){
  const = Coef[, ncol(Coef)]
  Coef.noc= Coef[,-ncol(Coef)]                      # 刪掉 const
}else{
  const = matrix(0, kk, 1)
  Coef.noc = Coef
}

Theta.unit= VAR.Theta(Coef, h, BSigma.u, CONST)$unit # 估算 Theta.unit
Theta.std = VAR.Theta(Coef, h, BSigma.u, CONST)$std  # 估算 Theta.std

# dm.U <- U-mean(U)
dm.U <- U

N = 2000 #重抽次數
Theta.unit.sim = vector("list", N)
Theta.std.sim  = vector("list", N)

# check dimension
print("check dimensionality")
dim(ddX); dim(Coef.noc); dim(dm.U)


# 存N次重抽的IRF
df_IRF.sim <- array(NA, c(hrz+1,kk^2,N)) #dimensions are: Time Period, Number of shock interacts with variables, page (number of Bootstrap resamplings)
counter <- 1
while(TRUE){
  
  #cat("Now, there are ", counter-1, " sets of resamples.\n")
  Y.sim = matrix(0, nrow = T.total, ncol = kk)          # Y.sim = 0 #pseudo time series
  Y.sim[c(1:VAR.P),] = By[c(1:VAR.P), ] #initial values
  
  boot.number = sample(c(1:T), replace = TRUE)      # Step 3 取出放回
  U.sim = dm.U[boot.number,]
  
  # predicted values given the above initial values
  last.y= c(t(By[VAR.P:1,]))
  for(ii in 1:T){
    last.y = last.y[1:(kk*VAR.P)]
    Y.sim[ii+VAR.P, ] = Coef.noc %*% last.y + const + U.sim[ii,]      # Step 4 模擬資料
    last.y = c(Y.sim[ii+VAR.P,], last.y)
  }
  
  #   Y.sim[-c(1:VAR.P),] <- matrix(const, nrow = T.total-VAR.P, ncol = kk, byrow = T) + ddX %*% t(Coef.noc) + U.sim
  
  
  #`Y.sim` is the pseudo time series
  # Step 5 重新估算SVAR
  
  ### SVAR.sim Start ###
  
  Y_pseudo     = VAR.Y(Y.sim, VAR.P)        # 設定 Y
  X_pseudo     = VAR.X(Y.sim, VAR.P)        # 設定 X
  Coef.OLS_pseudo    = VAR.OLS(Y_pseudo, X_pseudo, CONST)
  Sigma.OLS_pseudo   = VAR.Sigma.OLS(Y_pseudo, X_pseudo, Coef.OLS_pseudo, CONST)
  C.Prime_pseudo <- chol(Sigma.OLS_pseudo)
  C_pseudo <- t(C.Prime_pseudo)
  B0_pseudo <- diag(diag(C_pseudo), ncol = kk, nrow = kk)
  A0_pseudo <- B0_pseudo %*% solve(C_pseudo)
  SVAR_AB_est.sim <- list("A0.svar" = A0_pseudo, "B0.svar" = B0_pseudo)
  SVAR_AB_IRF.sim <- VAR.svarirf.AB(Y.sim, VAR.P, Amat, Bmat, h = hrz, CONST, SVAR_AB_est = SVAR_AB_est.sim)
  
  # 5*5個圖的time series
  df_IRF_plot.sim <- matrix(NA, hrz+1, kk^2) #%>% as.tibble()
  # df_IRF.sim <- array(1:(120*25*N), c(120,25,N))
  # df_IRF.sim[2,1,1] # slicing
  
  h <- 0 # h表示第幾期的IRF
  for(period in SVAR_AB_IRF.sim){
    k <- 0 # k表示把5*5的矩陣攤平到25個col的df時，要攤到第幾個columns上
    h <- h+1 # h表示第幾期的IRF
    for(j in 1:kk){
      for(i in 1:kk){
        k <- k+1 # k表示把5*5的矩陣攤平到25個col的df時，要攤到第幾個columns上
        df_IRF_plot.sim[h,k] <- period[i,j]
      }
    }
  }
  # 把這一次重抽得到的IRF append進`df_IRF.sim`中
  df_IRF.sim[,,counter] <- df_IRF_plot.sim
  ### SVAR.sim Ends ###
  if(counter>=N){
    break
  }
  counter <- counter+1
}
# Save
saveRDS(df_IRF.sim, file = "df_IRF.sim.rds")

df_IRF.sim <- read_rds("df_IRF.sim.rds")

# 看某一頁
head(df_IRF.sim[,,1000])
print(sum(is.na(df_IRF.sim)))

# 畫IRF & Bootstrap C.I.
df_IRF_plot.BS.L <- matrix(NA, nrow = hrz+1, ncol = kk^2)
df_IRF_plot.BS.U <- matrix(NA, nrow = hrz+1, ncol = kk^2)
df_IRF_plot.BS.Median <- matrix(NA, nrow = hrz+1, ncol = kk^2)
df_IRF_plot.BS.Mean <- matrix(NA, nrow = hrz+1, ncol = kk^2)
for(col in 1:(kk^2)){
  for(row in 1:(hrz+1) ){
    df_IRF_plot.BS.L[row,col] <- quantile(df_IRF.sim[row,col,], probs = 0.025)
    df_IRF_plot.BS.U[row,col] <- quantile(df_IRF.sim[row,col,], probs = 0.975)
    df_IRF_plot.BS.Median[row,col] <- quantile(df_IRF.sim[row,col,], probs = 0.5)
    df_IRF_plot.BS.Mean[row,col] <- mean(df_IRF.sim[row,col,])
  }
}

df_IRF_plot.BS.L <- df_IRF_plot.BS.L %>% as_tibble()
df_IRF_plot.BS.U <- df_IRF_plot.BS.U %>% as_tibble()
df_IRF_plot.BS.Median <- df_IRF_plot.BS.Median %>% as_tibble()
df_IRF_plot.BS.Mean <- df_IRF_plot.BS.Mean %>% as_tibble()

ind <- 0
for(i in 1:kk){
  for(j in 1:kk){
    ind <- ind+1
    nam <- paste("shock", j, "y", i, sep = '')
    assign(nam, bind_cols(df_IRF_plot.BS.L[ind], df_IRF_plot.BS.U[ind],
                          df_IRF_plot.BS.Median[ind], df_IRF_plot.BS.Mean[ind],
                          df_IRF_plot[ind]))
    # 改名
    evalStr <- paste0("colnames(", nam, ") <- c('Lower', 'Upper', 'Median', 'Mean', 'Actual')")
    eval(parse(text=evalStr))
    # 圖層
    evalStr <- paste0("p", ind, " <- ", "ggplot(",nam,") +geom_hline(yintercept=0, color = 'grey')+ geom_line(aes(x = 1:nrow(", nam, "), y = Lower), linetype = 'dashed', col='red')+geom_line(aes(x = 1:nrow(", nam, "), y = Upper), linetype = 'dashed', col='red')+geom_line(aes(x = 1:nrow(", nam, "), y = Median), col = 'Blue')")
    eval(parse(text=evalStr))
  }
}  

Text_Size_Theme = theme(
  axis.title.x = element_text(size = 12),
  axis.text.x = element_text(size = 12),
  axis.title.y = element_text(size = 12),
  axis.title = element_text(size = 12),
  plot.title = element_text(size=12))

add_layout <- function(p = p1,
                       response_of = "this variable", react_to = "the shock"){
  title_text <- paste("Response of", response_of, "to", react_to, sep = ' ')
  plot.return <- p + labs(x = "Period", y = "", title = title_text) + Text_Size_Theme
  return(plot.return)
}

vlb_list = c("Stock Price",
             "House Price",
             "Employment"
)
p1 <- add_layout(p1, vlb_list[1], "Stock shock")
p2 <- add_layout(p2, vlb_list[2], "Stock shock")
p3 <- add_layout(p3, vlb_list[3], "Stock shock")

p4 <- add_layout(p4, vlb_list[1], "Property Price shock")
p5 <- add_layout(p5, vlb_list[2], "Property Price shock")
p6 <- add_layout(p6, vlb_list[3], "Property Price shock")

p7 <- add_layout(p7, vlb_list[1], "Employment shock")
p8 <- add_layout(p8, vlb_list[2], "Employment shock")
p9 <- add_layout(p9, vlb_list[3], "Employment shock")

multiplot(p1,p2,p3,p4,p5,p6,
          p7,p8,p9,cols = 3)

multiplot(p1,p2,p3, cols = 1)
ggsave(filename = "/Users/caiyawei/Desktop/VAR/IRF_shock1.png", 
       plot = multiplot(p1,p2,p3, cols = 1),
       width = 30, height = 20, units = "cm",
       device = "png")

multiplot(p4,p5,p6, cols = 1)
ggsave(filename = "/Users/caiyawei/Desktop/VAR/IRF_shock2.png", 
       plot = multiplot(p4,p5,p6, cols = 1),
       width = 30, height = 20, units = "cm",
       device = "png")

multiplot(p7,p8,p9, cols = 1)
ggsave(filename = "/Users/caiyawei/Desktop/VAR/IRF_shock3.png", 
       plot = multiplot(p7,p8,p9, cols = 1),
       width = 30, height = 20, units = "cm",
       device = "png")

### Obtain the exogenous shock###
Coef.noc
Coef.OLS
const

# check the order of variables is correct (stock shock) (employment)
lm(stock_price ~ lag(stock_price)+lag(property_price)+lag(employment)+
     +lag(stock_price, n=2)+lag(property_price, n=2)+lag(employment, n=2)+lag(stock_price, n=3)+lag(property_price, n=3)+lag(employment, n=3)+lag(stock_price, n=4)+lag(property_price, n=4)+lag(employment, n=4),
   data = data)

get_residual <- function(vlb_order = 1){
  X.Mat <- By %>% as_tibble() %>%
    mutate(V1 = lag(stock_price),
           V2 = lag(property_price),
           V3 = lag(employment),
           V4 = lag(stock_price, n = 2),
           V5 = lag(property_price, n = 2),
           V6 = lag(employment, n = 2),
           V7 = lag(stock_price, n = 3),
           V8 = lag(property_price, n = 3),
           V9 = lag(employment, n = 3),
           V10 = lag(stock_price, n = 4),
           V11 = lag(property_price, n = 4),
           V12 = lag(employment, n = 4),
           V13 = rep(1, nrow(By))
    ) %>%
    select(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13
    ) %>%
    drop_na() %>%
    as.matrix
  Y.Mat <- By[(VAR.P+1):nrow(By),vlb_order] %>% as.matrix
  Y.Mat_hat <- X.Mat %*% (Coef.OLS[vlb_order, ] %>% as.matrix)
  residual <- Y.Mat - Y.Mat_hat
  return(residual)
}

# resid1 <- get_residual(1)
for(i in 1:3){
  assign(paste0("resid_",i), get_residual(i))
}

resid.Mat <- cbind(resid_1, resid_2, resid_3)
error.Mat <- (solve(B0) %*% A0 %*% t(resid.Mat)) %>% t()
error.Mat[,1] %>% hist

# Extract the stock shock

error.Mat %>% as_tibble() %>%
  select(e1 = 1, e2 = 2, e3 = 3) %>%
  mutate(no = row_number()) %>%
  ggplot()+
  geom_line(aes(no, e1))

data.err <- data %>% 
  bind_cols(e1 = c(rep(NA, VAR.P), error.Mat[,1])) %>%
  mutate(e1_pos = fifelse(e1 > 0, e1, 0),
         e1_neg = fifelse(e1 < 0, e1, 0))
data.err %>% head(10)

data.err %>% ggplot()+
  geom_line(aes(year, e1))

# check the order of variables is correct (house shock) (employment)
lm(property_price ~ lag(stock_price)+lag(property_price)+lag(employment)+
     +lag(stock_price, n=2)+lag(property_price, n=2)+lag(employment, n=2)+lag(stock_price, n=3)+lag(property_price, n=3)+lag(employment, n=3)+lag(stock_price, n=4)+lag(property_price, n=4)+lag(employment, n=4),
   data = data)

get_residual <- function(vlb_order = 1){
  X.Mat <- By %>% as_tibble() %>%
    mutate(V1 = lag(stock_price),
           V2 = lag(property_price),
           V3 = lag(employment),
           V4 = lag(stock_price, n = 2),
           V5 = lag(property_price, n = 2),
           V6 = lag(employment, n = 2),
           V7 = lag(stock_price, n = 3),
           V8 = lag(property_price, n = 3),
           V9 = lag(employment, n = 3),
           V10 = lag(stock_price, n = 4),
           V11 = lag(property_price, n = 4),
           V12 = lag(employment, n = 4),
           V13 = rep(1, nrow(By))
    ) %>%
    select(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13
    ) %>%
    drop_na() %>%
    as.matrix
  Y.Mat <- By[(VAR.P+1):nrow(By),vlb_order] %>% as.matrix
  Y.Mat_hat <- X.Mat %*% (Coef.OLS[vlb_order, ] %>% as.matrix)
  residual <- Y.Mat - Y.Mat_hat
  return(residual)
}

# resid1 <- get_residual(1)
for(i in 1:3){
  assign(paste0("resid_",i), get_residual(i))
}

resid.Mat <- cbind(resid_1, resid_2, resid_3)
error.Mat <-(solve(B0) %*% A0 %*% t(resid.Mat)) %>% t()
error.Mat[,2] %>% hist

# Extract the house shock

error.Mat %>% as_tibble() %>%
  select(e1 = 1, e2 = 2, e3 = 3) %>%
  mutate(no = row_number()) %>%
  ggplot()+
  geom_line(aes(no, e2))

data.err2 <- data %>% 
  bind_cols(e2 = c(rep(NA, VAR.P), error.Mat[,2])) %>%
  mutate(e2_pos = fifelse(e2 > 0, e2, 0),
         e2_neg = fifelse(e2 < 0, e2, 0))
data.err2 %>% head(10)

data.err2 %>% ggplot()+
  geom_line(aes(year, e2))

#data.err["Period"] = data$year - 7
data.err["e2"] = data.err2$e2
### Regression ###
#reg1 <- lm(em_d ~ Period + e1 + e2 + interest_rate_d, data = data.err)
#summary(reg1)

###衝擊也減去1991Q4###
data.err["e2"] = data.err["e2"] - data.err$e2[6]
data.err["e1"] = data.err["e1"] - data.err$e1[6]
reg1 <- lm(em_d ~ Period + e1 + e2 + interest_rate_d, data = data.err)
summary(reg1)

### VAR CPI in Order 1 ###

### DATA CPI ###
file = "/Users/caiyawei/Desktop/Japan Data/SVAR_forR.xlsx"
data_cpi = read_excel(file)
data_cpi = data_cpi[-c(107),]

By <- data_cpi %>% 
  select(stock_price, property_price, cpi) %>% 
  as.matrix

### RAW DATA ###
dim(By)
kk <- dim(By)[2]
draw_raw <- function(vlb_name = ""){
  vlb_name <- sym(vlb_name)
  plot.return <- data_cpi %>%
    ggplot()+
    geom_line(aes(year, !!vlb_name))
  return(plot.return)
}

p1 <- draw_raw("stock_price")
p2 <- draw_raw("property_price")
p3 <- draw_raw("cpi")

multiplot(p1, p2, p3)

### Setting Model ###
#----- 模型設定 -----#
VAR.P = 4                       # 最大的落後項數
CONST = TRUE                    # 是否有常數項
Y     = VAR.Y(By, VAR.P)        # 設定 Y
X     = VAR.X(By, VAR.P)        # 設定 X

hrz=39 # the length of response

### Reduced Form VAR ###
(Coef.OLS    = VAR.OLS(Y, X, CONST)                  )
(Sigma.OLS   = VAR.Sigma.OLS(Y, X, Coef.OLS, CONST)  )
(Sigma.MLE   = VAR.Sigma.MLE(Y, X, Coef.OLS, CONST))

VAR.P = 4 # By AIC

Amat = diag(3)
# Identification Conditions

Bmat = diag(3)


Amat;Bmat

### Estimate A^hat and B^hat ###
#C.Prime <- chol(Sigma.OLS)
#C <- t(C.Prime)
#C

### Solving system of linear equations ###
B0 <-Bmat


A0 <-Amat

SVAR_AB_est <- list("A0.svar" = A0, "B0.svar" = B0)

### IRF Without Bootstrap CI ###
### IRF
SVAR_AB_IRF <- VAR.svarirf.AB(By, VAR.P, Amat, Bmat, h = hrz, CONST, SVAR_AB_est = SVAR_AB_est)

# 6*6個圖的time series
df_IRF_plot <- matrix(NA, hrz+1, kk^2) #%>% as.tibble() ## hrz+1
#dim(df_IRF_plot)
h <- 0 # h表示第幾期的IRF
for(period in SVAR_AB_IRF){
  k <- 0 # k表示把5*5的矩陣攤平到25個col的df時，要攤到第幾個columns上
  h <- h+1 # h表示第幾期的IRF
  for(j in 1:kk){
    for(i in 1:kk){
      k <- k+1 # k表示把5*5的矩陣攤平到25個col的df時，要攤到第幾個columns上
      df_IRF_plot[h,k] <- period[i,j]
    }
  }
}
df_IRF_plot <- df_IRF_plot %>% as_tibble()

kk*1:kk

#output entire table
IRF_TABLE <- df_IRF_plot[,kk*1:kk] %>% select(stock_price = 1,
                                              property_price = 2,
                                              cpi = 3)
write.table(IRF_TABLE, file = "/Users/caiyawei/Desktop/IRF_TABLE3.csv", sep = ",", row.names = FALSE)

draw_IRF <- function(df = df_IRF_plot, V1 = 1){
  V1 <- paste0("V", V1) %>% sym()
  plot.return <- ggplot(df) + 
    geom_line(aes(x = 1:nrow(df), y = !!V1))
  return(plot.return)
}

for(i in 1:(kk^2)){
  assign(paste0("p",i), draw_IRF(df_IRF_plot, i))
}

multiplot(p1,p2,p3,p4,p5,p6,
          p7,p8,p9,cols = 3)

### IRF (Bootstrap C.I.) ###
lower = 0.025                                        # 控制成 95% CI
upper = 1-lower
kk = ncol(By)
ddY = VAR.ddY(By, VAR.P)
ddX = VAR.ddX(By, VAR.P)

# dim(ddY); dim(ddX)

T   = nrow(ddY)
T.total= nrow(By)
Ik  = diag(rep(1, kk))
# 16 coef if 4 variables; 55 coef if 5 variables
Coef = t(VAR.EbyE(ddY, ddX, CONST)$ddA)              # Step 1 估計模型
# residuals
U    = VAR.EbyE(ddY, ddX, CONST)$ddU
BSigma.u = VAR.ddSigma.OLS(ddY, ddX, CONST)
if(CONST == TRUE){
  const = Coef[, ncol(Coef)]
  Coef.noc= Coef[,-ncol(Coef)]                      # 刪掉 const
}else{
  const = matrix(0, kk, 1)
  Coef.noc = Coef
}

Theta.unit= VAR.Theta(Coef, h, BSigma.u, CONST)$unit # 估算 Theta.unit
Theta.std = VAR.Theta(Coef, h, BSigma.u, CONST)$std  # 估算 Theta.std

# dm.U <- U-mean(U)
dm.U <- U

N = 2000 #重抽次數
Theta.unit.sim = vector("list", N)
Theta.std.sim  = vector("list", N)

# check dimension
print("check dimensionality")
dim(ddX); dim(Coef.noc); dim(dm.U)


# 存N次重抽的IRF
df_IRF.sim <- array(NA, c(hrz+1,kk^2,N)) #dimensions are: Time Period, Number of shock interacts with variables, page (number of Bootstrap resamplings)
counter <- 1
while(TRUE){
  
  #cat("Now, there are ", counter-1, " sets of resamples.\n")
  Y.sim = matrix(0, nrow = T.total, ncol = kk)          # Y.sim = 0 #pseudo time series
  Y.sim[c(1:VAR.P),] = By[c(1:VAR.P), ] #initial values
  
  boot.number = sample(c(1:T), replace = TRUE)      # Step 3 取出放回
  U.sim = dm.U[boot.number,]
  
  # predicted values given the above initial values
  last.y= c(t(By[VAR.P:1,]))
  for(ii in 1:T){
    last.y = last.y[1:(kk*VAR.P)]
    Y.sim[ii+VAR.P, ] = Coef.noc %*% last.y + const + U.sim[ii,]      # Step 4 模擬資料
    last.y = c(Y.sim[ii+VAR.P,], last.y)
  }
  
  #   Y.sim[-c(1:VAR.P),] <- matrix(const, nrow = T.total-VAR.P, ncol = kk, byrow = T) + ddX %*% t(Coef.noc) + U.sim
  
  
  #`Y.sim` is the pseudo time series
  # Step 5 重新估算SVAR
  
  ### SVAR.sim Start ###
  
  Y_pseudo     = VAR.Y(Y.sim, VAR.P)        # 設定 Y
  X_pseudo     = VAR.X(Y.sim, VAR.P)        # 設定 X
  Coef.OLS_pseudo    = VAR.OLS(Y_pseudo, X_pseudo, CONST)
  Sigma.OLS_pseudo   = VAR.Sigma.OLS(Y_pseudo, X_pseudo, Coef.OLS_pseudo, CONST)
  C.Prime_pseudo <- chol(Sigma.OLS_pseudo)
  C_pseudo <- t(C.Prime_pseudo)
  B0_pseudo <- diag(diag(C_pseudo), ncol = kk, nrow = kk)
  A0_pseudo <- B0_pseudo %*% solve(C_pseudo)
  SVAR_AB_est.sim <- list("A0.svar" = A0_pseudo, "B0.svar" = B0_pseudo)
  SVAR_AB_IRF.sim <- VAR.svarirf.AB(Y.sim, VAR.P, Amat, Bmat, h = hrz, CONST, SVAR_AB_est = SVAR_AB_est.sim)
  
  # 5*5個圖的time series
  df_IRF_plot.sim <- matrix(NA, hrz+1, kk^2) #%>% as.tibble()
  # df_IRF.sim <- array(1:(120*25*N), c(120,25,N))
  # df_IRF.sim[2,1,1] # slicing
  
  h <- 0 # h表示第幾期的IRF
  for(period in SVAR_AB_IRF.sim){
    k <- 0 # k表示把5*5的矩陣攤平到25個col的df時，要攤到第幾個columns上
    h <- h+1 # h表示第幾期的IRF
    for(j in 1:kk){
      for(i in 1:kk){
        k <- k+1 # k表示把5*5的矩陣攤平到25個col的df時，要攤到第幾個columns上
        df_IRF_plot.sim[h,k] <- period[i,j]
      }
    }
  }
  # 把這一次重抽得到的IRF append進`df_IRF.sim`中
  df_IRF.sim[,,counter] <- df_IRF_plot.sim
  ### SVAR.sim Ends ###
  if(counter>=N){
    break
  }
  counter <- counter+1
}
# Save
saveRDS(df_IRF.sim, file = "df_IRF.sim.rds")

df_IRF.sim <- read_rds("df_IRF.sim.rds")

# 看某一頁
head(df_IRF.sim[,,1000])
print(sum(is.na(df_IRF.sim)))

# 畫IRF & Bootstrap C.I.
df_IRF_plot.BS.L <- matrix(NA, nrow = hrz+1, ncol = kk^2)
df_IRF_plot.BS.U <- matrix(NA, nrow = hrz+1, ncol = kk^2)
df_IRF_plot.BS.Median <- matrix(NA, nrow = hrz+1, ncol = kk^2)
df_IRF_plot.BS.Mean <- matrix(NA, nrow = hrz+1, ncol = kk^2)
for(col in 1:(kk^2)){
  for(row in 1:(hrz+1) ){
    df_IRF_plot.BS.L[row,col] <- quantile(df_IRF.sim[row,col,], probs = 0.025)
    df_IRF_plot.BS.U[row,col] <- quantile(df_IRF.sim[row,col,], probs = 0.975)
    df_IRF_plot.BS.Median[row,col] <- quantile(df_IRF.sim[row,col,], probs = 0.5)
    df_IRF_plot.BS.Mean[row,col] <- mean(df_IRF.sim[row,col,])
  }
}

df_IRF_plot.BS.L <- df_IRF_plot.BS.L %>% as_tibble()
df_IRF_plot.BS.U <- df_IRF_plot.BS.U %>% as_tibble()
df_IRF_plot.BS.Median <- df_IRF_plot.BS.Median %>% as_tibble()
df_IRF_plot.BS.Mean <- df_IRF_plot.BS.Mean %>% as_tibble()

ind <- 0
for(i in 1:kk){
  for(j in 1:kk){
    ind <- ind+1
    nam <- paste("shock", j, "y", i, sep = '')
    assign(nam, bind_cols(df_IRF_plot.BS.L[ind], df_IRF_plot.BS.U[ind],
                          df_IRF_plot.BS.Median[ind], df_IRF_plot.BS.Mean[ind],
                          df_IRF_plot[ind]))
    # 改名
    evalStr <- paste0("colnames(", nam, ") <- c('Lower', 'Upper', 'Median', 'Mean', 'Actual')")
    eval(parse(text=evalStr))
    # 圖層
    evalStr <- paste0("p", ind, " <- ", "ggplot(",nam,") +geom_hline(yintercept=0, color = 'grey')+ geom_line(aes(x = 1:nrow(", nam, "), y = Lower), linetype = 'dashed', col='red')+geom_line(aes(x = 1:nrow(", nam, "), y = Upper), linetype = 'dashed', col='red')+geom_line(aes(x = 1:nrow(", nam, "), y = Median), col = 'Blue')")
    eval(parse(text=evalStr))
  }
}  

Text_Size_Theme = theme(
  axis.title.x = element_text(size = 12),
  axis.text.x = element_text(size = 12),
  axis.title.y = element_text(size = 12),
  axis.title = element_text(size = 12),
  plot.title = element_text(size=12))

add_layout <- function(p = p1,
                       response_of = "this variable", react_to = "the shock"){
  title_text <- paste("Response of", response_of, "to", react_to, sep = ' ')
  plot.return <- p + labs(x = "Period", y = "", title = title_text) + Text_Size_Theme
  return(plot.return)
}

vlb_list = c("Stock Price",
             "House Price",
             "CPI"
)
p1 <- add_layout(p1, vlb_list[1], "Stock shock")
p2 <- add_layout(p2, vlb_list[2], "Stock shock")
p3 <- add_layout(p3, vlb_list[3], "Stock shock")

p4 <- add_layout(p4, vlb_list[1], "Property Price shock")
p5 <- add_layout(p5, vlb_list[2], "Property Price shock")
p6 <- add_layout(p6, vlb_list[3], "Property Price shock")

p7 <- add_layout(p7, vlb_list[1], "CPI shock")
p8 <- add_layout(p8, vlb_list[2], "CPI shock")
p9 <- add_layout(p9, vlb_list[3], "CPI shock")

multiplot(p1,p2,p3,p4,p5,p6,
          p7,p8,p9,cols = 3)

multiplot(p1,p2,p3, cols = 1)
ggsave(filename = "/Users/caiyawei/Desktop/VAR/IRF_shock_cpi1.png", 
       plot = multiplot(p1,p2,p3, cols = 1),
       width = 30, height = 20, units = "cm",
       device = "png")

multiplot(p4,p5,p6, cols = 1)
ggsave(filename = "/Users/caiyawei/Desktop/VAR/IRF_shock_cpi2.png", 
       plot = multiplot(p4,p5,p6, cols = 1),
       width = 30, height = 20, units = "cm",
       device = "png")

multiplot(p7,p8,p9, cols = 1)
ggsave(filename = "/Users/caiyawei/Desktop/VAR/IRF_shock_cpi3.png", 
       plot = multiplot(p7,p8,p9, cols = 1),
       width = 30, height = 20, units = "cm",
       device = "png")

### Obtain the exogenous shock###
Coef.noc
Coef.OLS
const

# check the order of variables is correct (stock shock) (cpi)
lm(stock_price ~ lag(stock_price)+lag(property_price)+lag(cpi)+
     +lag(stock_price, n=2)+lag(property_price, n=2)+lag(cpi, n=2)+lag(stock_price, n=3)+lag(property_price, n=3)+lag(cpi, n=3)+lag(stock_price, n=4)+lag(property_price, n=4)+lag(cpi, n=4),
   data = data_cpi)

get_residual <- function(vlb_order = 1){
  X.Mat <- By %>% as_tibble() %>%
    mutate(V1 = lag(stock_price),
           V2 = lag(property_price),
           V3 = lag(cpi),
           V4 = lag(stock_price, n = 2),
           V5 = lag(property_price, n = 2),
           V6 = lag(cpi, n = 2),
           V7 = lag(stock_price, n = 3),
           V8 = lag(property_price, n = 3),
           V9 = lag(cpi, n = 3),
           V10 = lag(stock_price, n = 4),
           V11 = lag(property_price, n = 4),
           V12 = lag(cpi, n = 4),
           V13 = rep(1, nrow(By))
    ) %>%
    select(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13
    ) %>%
    drop_na() %>%
    as.matrix
  Y.Mat <- By[(VAR.P+1):nrow(By),vlb_order] %>% as.matrix
  Y.Mat_hat <- X.Mat %*% (Coef.OLS[vlb_order, ] %>% as.matrix)
  residual <- Y.Mat - Y.Mat_hat
  return(residual)
}

# resid1 <- get_residual(1)
for(i in 1:3){
  assign(paste0("resid_",i), get_residual(i))
}

resid.Mat <- cbind(resid_1, resid_2, resid_3)
error.Mat <- (solve(B0) %*% A0 %*% t(resid.Mat)) %>% t()
error.Mat[,1] %>% hist

# Extract the stock shock

error.Mat %>% as_tibble() %>%
  select(e1 = 1, e2 = 2, e3 = 3) %>%
  mutate(no = row_number()) %>%
  ggplot()+
  geom_line(aes(no, e1))

data.err <- data_cpi %>% 
  bind_cols(e1 = c(rep(NA, VAR.P), error.Mat[,1])) %>%
  mutate(e1_pos = fifelse(e1 > 0, e1, 0),
         e1_neg = fifelse(e1 < 0, e1, 0))
data.err %>% head(10)

data.err %>% ggplot()+
  geom_line(aes(year, e1))

# check the order of variables is correct (house shock) (cpi)
lm(property_price ~ lag(stock_price)+lag(property_price)+lag(cpi)+
     +lag(stock_price, n=2)+lag(property_price, n=2)+lag(cpi, n=2)+lag(stock_price, n=3)+lag(property_price, n=3)+lag(cpi, n=3)+lag(stock_price, n=4)+lag(property_price, n=4)+lag(cpi, n=4),
   data = data_cpi)

get_residual <- function(vlb_order = 1){
  X.Mat <- By %>% as_tibble() %>%
    mutate(V1 = lag(stock_price),
           V2 = lag(property_price),
           V3 = lag(cpi),
           V4 = lag(stock_price, n = 2),
           V5 = lag(property_price, n = 2),
           V6 = lag(cpi, n = 2),
           V7 = lag(stock_price, n = 3),
           V8 = lag(property_price, n = 3),
           V9 = lag(cpi, n = 3),
           V10 = lag(stock_price, n = 4),
           V11 = lag(property_price, n = 4),
           V12 = lag(cpi, n = 4),
           V13 = rep(1, nrow(By))
    ) %>%
    select(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13
    ) %>%
    drop_na() %>%
    as.matrix
  Y.Mat <- By[(VAR.P+1):nrow(By),vlb_order] %>% as.matrix
  Y.Mat_hat <- X.Mat %*% (Coef.OLS[vlb_order, ] %>% as.matrix)
  residual <- Y.Mat - Y.Mat_hat
  return(residual)
}

# resid1 <- get_residual(1)
for(i in 1:3){
  assign(paste0("resid_",i), get_residual(i))
}

resid.Mat <- cbind(resid_1, resid_2, resid_3)
error.Mat <-  (solve(B0) %*% A0 %*% t(resid.Mat)) %>% t()
error.Mat[,2] %>% hist

# Extract the house shock

error.Mat %>% as_tibble() %>%
  select(e1 = 1, e2 = 2, e3 = 3) %>%
  mutate(no = row_number()) %>%
  ggplot()+
  geom_line(aes(no, e2))

data.err2 <- data_cpi %>% 
  bind_cols(e2 = c(rep(NA, VAR.P), error.Mat[,2])) %>%
  mutate(e2_pos = fifelse(e2 > 0, e2, 0),
         e2_neg = fifelse(e2 < 0, e2, 0))
data.err2 %>% head(10)

data.err2 %>% ggplot()+
  geom_line(aes(year, e2))


data.err["e2"] = data.err2$e2
### Regression ###
#reg3 <- lm(cpi_d ~ Period + e1 + e2 + interest_rate_d, data = data.err)
#summary(reg3)

###衝擊減掉1994Q4###
data.err["e2"] = data.err["e2"] - data.err$e2[6]
data.err["e1"] = data.err["e1"] - data.err$e1[6]
reg3 <- lm(cpi_d ~ Period + e1 + e2 + interest_rate_d, data = data.err)
summary(reg3)








### VAR DATA GDP ORDER 1 ###
file = "/Users/caiyawei/Desktop/Japan Data/SVAR_forR.xlsx"
data_gdp = read_excel(file)
data_gdp = data_gdp[-c(107),]



By <- data_gdp %>% 
  select(stock_price, property_price, gdp) %>% 
  as.matrix

### RAW DATA ###
dim(By)
kk <- dim(By)[2]
draw_raw <- function(vlb_name = ""){
  vlb_name <- sym(vlb_name)
  plot.return <- data_gdp %>%
    ggplot()+
    geom_line(aes(year, !!vlb_name))
  return(plot.return)
}

p1 <- draw_raw("stock_price")
p2 <- draw_raw("property_price")
p3 <- draw_raw("gdp")

multiplot(p1, p2, p3)

### Setting Model ###
#----- 模型設定 -----#
VAR.P = 4                       # 最大的落後項數
CONST = TRUE                    # 是否有常數項
Y     = VAR.Y(By, VAR.P)        # 設定 Y
X     = VAR.X(By, VAR.P)        # 設定 X

hrz=39 # the length of response

### Reduced Form VAR ###
(Coef.OLS    = VAR.OLS(Y, X, CONST)                  )
(Sigma.OLS   = VAR.Sigma.OLS(Y, X, Coef.OLS, CONST)  )
(Sigma.MLE   = VAR.Sigma.MLE(Y, X, Coef.OLS, CONST))

VAR.P = 4 # By AIC


Amat = diag(3)
# Identification Conditions

Bmat = diag(3)


Amat;Bmat

### Estimate A^hat and B^hat ###
#C.Prime <- chol(Sigma.OLS)
#C <- t(C.Prime)
#C

### Solving system of linear equations ###
B0 <-Bmat


A0 <-Amat



Amat;Bmat



SVAR_AB_est <- list("A0.svar" = A0, "B0.svar" = B0)

### IRF Without Bootstrap CI ###
### IRF
SVAR_AB_IRF <- VAR.svarirf.AB(By, VAR.P, Amat, Bmat, h = hrz, CONST, SVAR_AB_est = SVAR_AB_est)


# 6*6個圖的time series
df_IRF_plot <- matrix(NA, hrz+1, kk^2) #%>% as.tibble() ## hrz+1
#dim(df_IRF_plot)
h <- 0 # h表示第幾期的IRF
for(period in SVAR_AB_IRF){
  k <- 0 # k表示把5*5的矩陣攤平到25個col的df時，要攤到第幾個columns上
  h <- h+1 # h表示第幾期的IRF
  for(j in 1:kk){
    for(i in 1:kk){
      k <- k+1 # k表示把5*5的矩陣攤平到25個col的df時，要攤到第幾個columns上
      df_IRF_plot[h,k] <- period[i,j]
    }
  }
}
df_IRF_plot <- df_IRF_plot %>% as_tibble()

kk*1:kk





#output entire table
IRF_TABLE <- df_IRF_plot[,kk*1:kk] %>% select(stock_price = 1,
                                              property_price = 2,
                                              gdp = 3)
write.table(IRF_TABLE, file = "/Users/caiyawei/Desktop/IRF_TABLE2.csv", sep = ",", row.names = FALSE)

draw_IRF <- function(df = df_IRF_plot, V1 = 1){
  V1 <- paste0("V", V1) %>% sym()
  plot.return <- ggplot(df) + 
    geom_line(aes(x = 1:nrow(df), y = !!V1))
  return(plot.return)
}

for(i in 1:(kk^2)){
  assign(paste0("p",i), draw_IRF(df_IRF_plot, i))
}

multiplot(p1,p2,p3,p4,p5,p6,
          p7,p8,p9,cols = 3)

### IRF (Bootstrap C.I.) ###
lower = 0.025                                        # 控制成 95% CI
upper = 1-lower
kk = ncol(By)
ddY = VAR.ddY(By, VAR.P)
ddX = VAR.ddX(By, VAR.P)

# dim(ddY); dim(ddX)

T   = nrow(ddY)
T.total= nrow(By)
Ik  = diag(rep(1, kk))
# 16 coef if 4 variables; 55 coef if 5 variables
Coef = t(VAR.EbyE(ddY, ddX, CONST)$ddA)              # Step 1 估計模型
# residuals
U    = VAR.EbyE(ddY, ddX, CONST)$ddU
BSigma.u = VAR.ddSigma.OLS(ddY, ddX, CONST)
if(CONST == TRUE){
  const = Coef[, ncol(Coef)]
  Coef.noc= Coef[,-ncol(Coef)]                      # 刪掉 const
}else{
  const = matrix(0, kk, 1)
  Coef.noc = Coef
}

Theta.unit= VAR.Theta(Coef, h, BSigma.u, CONST)$unit # 估算 Theta.unit
Theta.std = VAR.Theta(Coef, h, BSigma.u, CONST)$std  # 估算 Theta.std

# dm.U <- U-mean(U)
dm.U <- U

N = 2000 #重抽次數
Theta.unit.sim = vector("list", N)
Theta.std.sim  = vector("list", N)

# check dimension
print("check dimensionality")
dim(ddX); dim(Coef.noc); dim(dm.U)


# 存N次重抽的IRF
df_IRF.sim <- array(NA, c(hrz+1,kk^2,N)) #dimensions are: Time Period, Number of shock interacts with variables, page (number of Bootstrap resamplings)
counter <- 1
while(TRUE){
  
  #cat("Now, there are ", counter-1, " sets of resamples.\n")
  Y.sim = matrix(0, nrow = T.total, ncol = kk)          # Y.sim = 0 #pseudo time series
  Y.sim[c(1:VAR.P),] = By[c(1:VAR.P), ] #initial values
  
  boot.number = sample(c(1:T), replace = TRUE)      # Step 3 取出放回
  U.sim = dm.U[boot.number,]
  
  # predicted values given the above initial values
  last.y= c(t(By[VAR.P:1,]))
  for(ii in 1:T){
    last.y = last.y[1:(kk*VAR.P)]
    Y.sim[ii+VAR.P, ] = Coef.noc %*% last.y + const + U.sim[ii,]      # Step 4 模擬資料
    last.y = c(Y.sim[ii+VAR.P,], last.y)
  }
  
  #   Y.sim[-c(1:VAR.P),] <- matrix(const, nrow = T.total-VAR.P, ncol = kk, byrow = T) + ddX %*% t(Coef.noc) + U.sim
  
  
  #`Y.sim` is the pseudo time series
  # Step 5 重新估算SVAR
  
  ### SVAR.sim Start ###
  
  Y_pseudo     = VAR.Y(Y.sim, VAR.P)        # 設定 Y
  X_pseudo     = VAR.X(Y.sim, VAR.P)        # 設定 X
  Coef.OLS_pseudo    = VAR.OLS(Y_pseudo, X_pseudo, CONST)
  Sigma.OLS_pseudo   = VAR.Sigma.OLS(Y_pseudo, X_pseudo, Coef.OLS_pseudo, CONST)
  C.Prime_pseudo <- chol(Sigma.OLS_pseudo)
  C_pseudo <- t(C.Prime_pseudo)
  B0_pseudo <- diag(diag(C_pseudo), ncol = kk, nrow = kk)
  A0_pseudo <- B0_pseudo %*% solve(C_pseudo)
  SVAR_AB_est.sim <- list("A0.svar" = A0_pseudo, "B0.svar" = B0_pseudo)
  SVAR_AB_IRF.sim <- VAR.svarirf.AB(Y.sim, VAR.P, Amat, Bmat, h = hrz, CONST, SVAR_AB_est = SVAR_AB_est.sim)
  
  # 5*5個圖的time series
  df_IRF_plot.sim <- matrix(NA, hrz+1, kk^2) #%>% as.tibble()
  # df_IRF.sim <- array(1:(120*25*N), c(120,25,N))
  # df_IRF.sim[2,1,1] # slicing
  
  h <- 0 # h表示第幾期的IRF
  for(period in SVAR_AB_IRF.sim){
    k <- 0 # k表示把5*5的矩陣攤平到25個col的df時，要攤到第幾個columns上
    h <- h+1 # h表示第幾期的IRF
    for(j in 1:kk){
      for(i in 1:kk){
        k <- k+1 # k表示把5*5的矩陣攤平到25個col的df時，要攤到第幾個columns上
        df_IRF_plot.sim[h,k] <- period[i,j]
      }
    }
  }
  # 把這一次重抽得到的IRF append進`df_IRF.sim`中
  df_IRF.sim[,,counter] <- df_IRF_plot.sim
  ### SVAR.sim Ends ###
  if(counter>=N){
    break
  }
  counter <- counter+1
}
# Save
saveRDS(df_IRF.sim, file = "df_IRF.sim.rds")

df_IRF.sim <- read_rds("df_IRF.sim.rds")

# 看某一頁
head(df_IRF.sim[,,1000])
print(sum(is.na(df_IRF.sim)))

# 畫IRF & Bootstrap C.I.
df_IRF_plot.BS.L <- matrix(NA, nrow = hrz+1, ncol = kk^2)
df_IRF_plot.BS.U <- matrix(NA, nrow = hrz+1, ncol = kk^2)
df_IRF_plot.BS.Median <- matrix(NA, nrow = hrz+1, ncol = kk^2)
df_IRF_plot.BS.Mean <- matrix(NA, nrow = hrz+1, ncol = kk^2)
for(col in 1:(kk^2)){
  for(row in 1:(hrz+1) ){
    df_IRF_plot.BS.L[row,col] <- quantile(df_IRF.sim[row,col,], probs = 0.025)
    df_IRF_plot.BS.U[row,col] <- quantile(df_IRF.sim[row,col,], probs = 0.975)
    df_IRF_plot.BS.Median[row,col] <- quantile(df_IRF.sim[row,col,], probs = 0.5)
    df_IRF_plot.BS.Mean[row,col] <- mean(df_IRF.sim[row,col,])
  }
}

df_IRF_plot.BS.L <- df_IRF_plot.BS.L %>% as_tibble()
df_IRF_plot.BS.U <- df_IRF_plot.BS.U %>% as_tibble()
df_IRF_plot.BS.Median <- df_IRF_plot.BS.Median %>% as_tibble()
df_IRF_plot.BS.Mean <- df_IRF_plot.BS.Mean %>% as_tibble()

ind <- 0
for(i in 1:kk){
  for(j in 1:kk){
    ind <- ind+1
    nam <- paste("shock", j, "y", i, sep = '')
    assign(nam, bind_cols(df_IRF_plot.BS.L[ind], df_IRF_plot.BS.U[ind],
                          df_IRF_plot.BS.Median[ind], df_IRF_plot.BS.Mean[ind],
                          df_IRF_plot[ind]))
    # 改名
    evalStr <- paste0("colnames(", nam, ") <- c('Lower', 'Upper', 'Median', 'Mean', 'Actual')")
    eval(parse(text=evalStr))
    # 圖層
    evalStr <- paste0("p", ind, " <- ", "ggplot(",nam,") +geom_hline(yintercept=0, color = 'grey')+ geom_line(aes(x = 1:nrow(", nam, "), y = Lower), linetype = 'dashed', col='red')+geom_line(aes(x = 1:nrow(", nam, "), y = Upper), linetype = 'dashed', col='red')+geom_line(aes(x = 1:nrow(", nam, "), y = Median), col = 'Blue')")
    eval(parse(text=evalStr))
  }
}  

Text_Size_Theme = theme(
  axis.title.x = element_text(size = 12),
  axis.text.x = element_text(size = 12),
  axis.title.y = element_text(size = 12),
  axis.title = element_text(size = 12),
  plot.title = element_text(size=12))

add_layout <- function(p = p1,
                       response_of = "this variable", react_to = "the shock"){
  title_text <- paste("Response of", response_of, "to", react_to, sep = ' ')
  plot.return <- p + labs(x = "Period", y = "", title = title_text) + Text_Size_Theme
  return(plot.return)
}

vlb_list = c("Stock Price",
             "House Price",
             "GDP"
)
p1 <- add_layout(p1, vlb_list[1], "Stock shock")
p2 <- add_layout(p2, vlb_list[2], "Stock shock")
p3 <- add_layout(p3, vlb_list[3], "Stock shock")

p4 <- add_layout(p4, vlb_list[1], "Property Price shock")
p5 <- add_layout(p5, vlb_list[2], "Property Price shock")
p6 <- add_layout(p6, vlb_list[3], "Property Price shock")

p7 <- add_layout(p7, vlb_list[1], "GDP shock")
p8 <- add_layout(p8, vlb_list[2], "GDP shock")
p9 <- add_layout(p9, vlb_list[3], "GDP shock")

multiplot(p1,p2,p3,p4,p5,p6,
          p7,p8,p9,cols = 3)

multiplot(p1,p2,p3, cols = 1)
ggsave(filename = "/Users/caiyawei/Desktop/VAR/IRF_shock_gdp1.png", 
       plot = multiplot(p1,p2,p3, cols = 1),
       width = 30, height = 20, units = "cm",
       device = "png")

multiplot(p4,p5,p6, cols = 1)
ggsave(filename = "/Users/caiyawei/Desktop/VAR/IRF_shock_gdp2.png", 
       plot = multiplot(p4,p5,p6, cols = 1),
       width = 30, height = 20, units = "cm",
       device = "png")

multiplot(p3,p6,p9, cols = 1)
ggsave(filename = "/Users/caiyawei/Desktop/VAR/IRF_shock_gdp33.png", 
       plot = multiplot(p7,p8,p9, cols = 1),
       width = 30, height = 20, units = "cm",
       device = "png")

### Obtain the exogenous shock###
Coef.noc
Coef.OLS
const

# check the order of variables is correct (stock shock) (gdp)
lm(stock_price ~ lag(stock_price)+lag(property_price)+lag(gdp)+
     +lag(stock_price, n=2)+lag(property_price, n=2)+lag(gdp, n=2)+lag(stock_price, n=3)+lag(property_price, n=3)+lag(gdp, n=3)+lag(stock_price, n=4)+lag(property_price, n=4)+lag(gdp, n=4),
   data = data_gdp)

get_residual <- function(vlb_order = 1){
  X.Mat <- By %>% as_tibble() %>%
    mutate(V1 = lag(stock_price),
           V2 = lag(property_price),
           V3 = lag(gdp),
           V4 = lag(stock_price, n = 2),
           V5 = lag(property_price, n = 2),
           V6 = lag(gdp, n = 2),
           V7 = lag(stock_price, n = 3),
           V8 = lag(property_price, n = 3),
           V9 = lag(gdp, n = 3),
           V10 = lag(stock_price, n = 4),
           V11 = lag(property_price, n = 4),
           V12 = lag(gdp, n = 4),
           V13 = rep(1, nrow(By))
    ) %>%
    select(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13
    ) %>%
    drop_na() %>%
    as.matrix
  Y.Mat <- By[(VAR.P+1):nrow(By),vlb_order] %>% as.matrix
  Y.Mat_hat <- X.Mat %*% (Coef.OLS[vlb_order, ] %>% as.matrix)
  residual <- Y.Mat - Y.Mat_hat
  return(residual)
}

# resid1 <- get_residual(1)
for(i in 1:3){
  assign(paste0("resid_",i), get_residual(i))
}

resid.Mat <- cbind(resid_1, resid_2, resid_3)
error.Mat <-  (solve(B0) %*% A0 %*% t(resid.Mat)) %>% t()
error.Mat[,1] %>% hist

# Extract the stock shock

error.Mat %>% as_tibble() %>%
  select(e1 = 1, e2 = 2, e3 = 3) %>%
  mutate(no = row_number()) %>%
  ggplot()+
  geom_line(aes(no, e1))

data.err <- data_gdp %>% 
  bind_cols(e1 = c(rep(NA, VAR.P), error.Mat[,1])) %>%
  mutate(e1_pos = fifelse(e1 > 0, e1, 0),
         e1_neg = fifelse(e1 < 0, e1, 0))
data.err %>% head(10)

data.err %>% ggplot()+
  geom_line(aes(year, e1))

# check the order of variables is correct (house shock) (gdp)
lm(property_price ~ lag(stock_price)+lag(property_price)+lag(gdp)+
     +lag(stock_price, n=2)+lag(property_price, n=2)+lag(gdp, n=2)+lag(stock_price, n=3)+lag(property_price, n=3)+lag(gdp, n=3)+lag(stock_price, n=4)+lag(property_price, n=4)+lag(gdp, n=4),
   data = data_gdp)

get_residual <- function(vlb_order = 1){
  X.Mat <- By %>% as_tibble() %>%
    mutate(V1 = lag(stock_price),
           V2 = lag(property_price),
           V3 = lag(gdp),
           V4 = lag(stock_price, n = 2),
           V5 = lag(property_price, n = 2),
           V6 = lag(gdp, n = 2),
           V7 = lag(stock_price, n = 3),
           V8 = lag(property_price, n = 3),
           V9 = lag(gdp, n = 3),
           V10 = lag(stock_price, n = 4),
           V11 = lag(property_price, n = 4),
           V12 = lag(gdp, n = 4),
           V13 = rep(1, nrow(By))
    ) %>%
    select(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13
    ) %>%
    drop_na() %>%
    as.matrix
  Y.Mat <- By[(VAR.P+1):nrow(By),vlb_order] %>% as.matrix
  Y.Mat_hat <- X.Mat %*% (Coef.OLS[vlb_order, ] %>% as.matrix)
  residual <- Y.Mat - Y.Mat_hat
  return(residual)
}

# resid1 <- get_residual(1)
for(i in 1:3){
  assign(paste0("resid_",i), get_residual(i))
}

resid.Mat <- cbind(resid_1, resid_2, resid_3)
error.Mat <- (solve(B0) %*% A0 %*% t(resid.Mat)) %>% t()
error.Mat[,2] %>% hist

# Extract the house shock

error.Mat %>% as_tibble() %>%
  select(e1 = 1, e2 = 2, e3 = 3) %>%
  mutate(no = row_number()) %>%
  ggplot()+
  geom_line(aes(no, e2))

data.err2 <- data_gdp %>% 
  bind_cols(e2 = c(rep(NA, VAR.P), error.Mat[,2])) %>%
  mutate(e2_pos = fifelse(e2 > 0, e2, 0),
         e2_neg = fifelse(e2 < 0, e2, 0))
data.err2 %>% head(10)

data.err2 %>% ggplot()+
  geom_line(aes(year, e2))


data.err["e2"] = data.err2$e2
### Regression ###
#reg2 <- lm(gdp_d~ Period + e1 + e2 + interest_rate_d, data = data.err)
#summary(reg2)

###衝擊減掉1994Q4###
data.err["e2"] = data.err["e2"] - data.err$e2[6]
data.err["e1"] = data.err["e1"] - data.err$e1[6]
reg2 <- lm(gdp_d~ Period + e1 + e2 + interest_rate_d, data = data.err)
summary(reg2)

texreg(list(reg1, reg2, reg3), caption = "Three VAR models in Order 1.", digits = 5)



###VAR Order2 House -> Stock -> Macro###
file = "/Users/caiyawei/Desktop/Japan Data/SVAR_forR.xlsx"
data = read_excel(file)
data = data[-c(107),]

By <- data %>% 
  select(property_price, stock_price, employment) %>% 
  as.matrix

### RAW DATA ###
dim(By)
kk <- dim(By)[2]
draw_raw <- function(vlb_name = ""){
  vlb_name <- sym(vlb_name)
  plot.return <- data %>%
    ggplot()+
    geom_line(aes(year, !!vlb_name))
  return(plot.return)
}

p1 <- draw_raw("property_price")
p2 <- draw_raw("stock_price")
p3 <- draw_raw("employment")

multiplot(p1, p2, p3)

### Setting Model ###
#----- 模型設定 -----#
VAR.P = 4                       # 最大的落後項數
CONST = TRUE                    # 是否有常數項
Y     = VAR.Y(By, VAR.P)        # 設定 Y
X     = VAR.X(By, VAR.P)        # 設定 X

hrz=39 # the length of response

### Reduced Form VAR ###
(Coef.OLS    = VAR.OLS(Y, X, CONST)                  )
(Sigma.OLS   = VAR.Sigma.OLS(Y, X, Coef.OLS, CONST)  )
(Sigma.MLE   = VAR.Sigma.MLE(Y, X, Coef.OLS, CONST))

VAR.P = 4 # By AIC

Amat = diag(3)
# Identification Conditions

Bmat = diag(3)


Amat;Bmat

### Estimate A^hat and B^hat ###
#C.Prime <- chol(Sigma.OLS)
#C <- t(C.Prime)
#C

### Solving system of linear equations ###
B0 <-Bmat


A0 <-Amat

SVAR_AB_est <- list("A0.svar" = A0, "B0.svar" = B0)

### IRF Without Bootstrap CI ###
### IRF
SVAR_AB_IRF <- VAR.svarirf.AB(By, VAR.P, Amat, Bmat, h = hrz, CONST, SVAR_AB_est = SVAR_AB_est)

kk*1:kk

#output entire table
IRF_TABLE <- df_IRF_plot[,kk*1:kk] %>% select(property_price = 1,
                                              stock_price = 2,
                                              employment = 3)
write.table(IRF_TABLE, file = "/Users/caiyawei/Desktop/IRF_TABLE_order2.csv", sep = ",", row.names = FALSE)

draw_IRF <- function(df = df_IRF_plot, V1 = 1){
  V1 <- paste0("V", V1) %>% sym()
  plot.return <- ggplot(df) + 
    geom_line(aes(x = 1:nrow(df), y = !!V1))
  return(plot.return)
}

for(i in 1:(kk^2)){
  assign(paste0("p",i), draw_IRF(df_IRF_plot, i))
}

multiplot(p1,p2,p3,p4,p5,p6,
          p7,p8,p9,cols = 3)

### IRF (Bootstrap C.I.) ###
lower = 0.025                                        # 控制成 95% CI
upper = 1-lower
kk = ncol(By)
ddY = VAR.ddY(By, VAR.P)
ddX = VAR.ddX(By, VAR.P)

# dim(ddY); dim(ddX)

T   = nrow(ddY)
T.total= nrow(By)
Ik  = diag(rep(1, kk))
# 16 coef if 4 variables; 55 coef if 5 variables
Coef = t(VAR.EbyE(ddY, ddX, CONST)$ddA)              # Step 1 估計模型
# residuals
U    = VAR.EbyE(ddY, ddX, CONST)$ddU
BSigma.u = VAR.ddSigma.OLS(ddY, ddX, CONST)
if(CONST == TRUE){
  const = Coef[, ncol(Coef)]
  Coef.noc= Coef[,-ncol(Coef)]                      # 刪掉 const
}else{
  const = matrix(0, kk, 1)
  Coef.noc = Coef
}

Theta.unit= VAR.Theta(Coef, h, BSigma.u, CONST)$unit # 估算 Theta.unit
Theta.std = VAR.Theta(Coef, h, BSigma.u, CONST)$std  # 估算 Theta.std

# dm.U <- U-mean(U)
dm.U <- U

N = 2000 #重抽次數
Theta.unit.sim = vector("list", N)
Theta.std.sim  = vector("list", N)

# check dimension
print("check dimensionality")
dim(ddX); dim(Coef.noc); dim(dm.U)


# 存N次重抽的IRF
df_IRF.sim <- array(NA, c(hrz+1,kk^2,N)) #dimensions are: Time Period, Number of shock interacts with variables, page (number of Bootstrap resamplings)
counter <- 1
while(TRUE){
  
  #cat("Now, there are ", counter-1, " sets of resamples.\n")
  Y.sim = matrix(0, nrow = T.total, ncol = kk)          # Y.sim = 0 #pseudo time series
  Y.sim[c(1:VAR.P),] = By[c(1:VAR.P), ] #initial values
  
  boot.number = sample(c(1:T), replace = TRUE)     
  U.sim = dm.U[boot.number,]
  
  # predicted values given the above initial values
  last.y= c(t(By[VAR.P:1,]))
  for(ii in 1:T){
    last.y = last.y[1:(kk*VAR.P)]
    Y.sim[ii+VAR.P, ] = Coef.noc %*% last.y + const + U.sim[ii,]      
    last.y = c(Y.sim[ii+VAR.P,], last.y)
  }
  
  #   Y.sim[-c(1:VAR.P),] <- matrix(const, nrow = T.total-VAR.P, ncol = kk, byrow = T) + ddX %*% t(Coef.noc) + U.sim
  
  
  #`Y.sim` is the pseudo time series
  
  ### SVAR.sim Start ###
  
  Y_pseudo     = VAR.Y(Y.sim, VAR.P)        
  X_pseudo     = VAR.X(Y.sim, VAR.P)        
  Coef.OLS_pseudo    = VAR.OLS(Y_pseudo, X_pseudo, CONST)
  Sigma.OLS_pseudo   = VAR.Sigma.OLS(Y_pseudo, X_pseudo, Coef.OLS_pseudo, CONST)
  C.Prime_pseudo <- chol(Sigma.OLS_pseudo)
  C_pseudo <- t(C.Prime_pseudo)
  B0_pseudo <- diag(diag(C_pseudo), ncol = kk, nrow = kk)
  A0_pseudo <- B0_pseudo %*% solve(C_pseudo)
  SVAR_AB_est.sim <- list("A0.svar" = A0_pseudo, "B0.svar" = B0_pseudo)
  SVAR_AB_IRF.sim <- VAR.svarirf.AB(Y.sim, VAR.P, Amat, Bmat, h = hrz, CONST, SVAR_AB_est = SVAR_AB_est.sim)
  
  df_IRF_plot.sim <- matrix(NA, hrz+1, kk^2) #%>% as.tibble()
  # df_IRF.sim <- array(1:(120*25*N), c(120,25,N))
  # df_IRF.sim[2,1,1] # slicing
  
  h <- 0
  for(period in SVAR_AB_IRF.sim){
    k <- 0 #
    h <- h+1 
    for(j in 1:kk){
      for(i in 1:kk){
        k <- k+1
        df_IRF_plot.sim[h,k] <- period[i,j]
      }
    }
  }
  df_IRF.sim[,,counter] <- df_IRF_plot.sim
  ### SVAR.sim Ends ###
  if(counter>=N){
    break
  }
  counter <- counter+1
}
# Save
saveRDS(df_IRF.sim, file = "df_IRF.sim.rds")

df_IRF.sim <- read_rds("df_IRF.sim.rds")

head(df_IRF.sim[,,1000])
print(sum(is.na(df_IRF.sim)))

# IRF & Bootstrap C.I.
df_IRF_plot.BS.L <- matrix(NA, nrow = hrz+1, ncol = kk^2)
df_IRF_plot.BS.U <- matrix(NA, nrow = hrz+1, ncol = kk^2)
df_IRF_plot.BS.Median <- matrix(NA, nrow = hrz+1, ncol = kk^2)
df_IRF_plot.BS.Mean <- matrix(NA, nrow = hrz+1, ncol = kk^2)
for(col in 1:(kk^2)){
  for(row in 1:(hrz+1) ){
    df_IRF_plot.BS.L[row,col] <- quantile(df_IRF.sim[row,col,], probs = 0.025)
    df_IRF_plot.BS.U[row,col] <- quantile(df_IRF.sim[row,col,], probs = 0.975)
    df_IRF_plot.BS.Median[row,col] <- quantile(df_IRF.sim[row,col,], probs = 0.5)
    df_IRF_plot.BS.Mean[row,col] <- mean(df_IRF.sim[row,col,])
  }
}

df_IRF_plot.BS.L <- df_IRF_plot.BS.L %>% as_tibble()
df_IRF_plot.BS.U <- df_IRF_plot.BS.U %>% as_tibble()
df_IRF_plot.BS.Median <- df_IRF_plot.BS.Median %>% as_tibble()
df_IRF_plot.BS.Mean <- df_IRF_plot.BS.Mean %>% as_tibble()

ind <- 0
for(i in 1:kk){
  for(j in 1:kk){
    ind <- ind+1
    nam <- paste("shock", j, "y", i, sep = '')
    assign(nam, bind_cols(df_IRF_plot.BS.L[ind], df_IRF_plot.BS.U[ind],
                          df_IRF_plot.BS.Median[ind], df_IRF_plot.BS.Mean[ind],
                          df_IRF_plot[ind]))
    # 改名
    evalStr <- paste0("colnames(", nam, ") <- c('Lower', 'Upper', 'Median', 'Mean', 'Actual')")
    eval(parse(text=evalStr))
    # 圖層
    evalStr <- paste0("p", ind, " <- ", "ggplot(",nam,") +geom_hline(yintercept=0, color = 'grey')+ geom_line(aes(x = 1:nrow(", nam, "), y = Lower), linetype = 'dashed', col='red')+geom_line(aes(x = 1:nrow(", nam, "), y = Upper), linetype = 'dashed', col='red')+geom_line(aes(x = 1:nrow(", nam, "), y = Median), col = 'Blue')")
    eval(parse(text=evalStr))
  }
}  

Text_Size_Theme = theme(
  axis.title.x = element_text(size = 12),
  axis.text.x = element_text(size = 12),
  axis.title.y = element_text(size = 12),
  axis.title = element_text(size = 12),
  plot.title = element_text(size=12))

add_layout <- function(p = p1,
                       response_of = "this variable", react_to = "the shock"){
  title_text <- paste("Response of", response_of, "to", react_to, sep = ' ')
  plot.return <- p + labs(x = "Period", y = "", title = title_text) + Text_Size_Theme
  return(plot.return)
}

vlb_list = c("House Price",
             "Stock Price",
             "Employment"
)
p1 <- add_layout(p1, vlb_list[1], "Property Price shock")
p2 <- add_layout(p2, vlb_list[2], "Property Price shock")
p3 <- add_layout(p3, vlb_list[3], "Property Price shock")

p4 <- add_layout(p4, vlb_list[1], "Stock shock")
p5 <- add_layout(p5, vlb_list[2], "Stock shock")
p6 <- add_layout(p6, vlb_list[3], "Stock shock")

p7 <- add_layout(p7, vlb_list[1], "Employment shock")
p8 <- add_layout(p8, vlb_list[2], "Employment shock")
p9 <- add_layout(p9, vlb_list[3], "Employment shock")

multiplot(p1,p2,p3,p4,p5,p6,
          p7,p8,p9,cols = 3)

multiplot(p1,p2,p3, cols = 1)
ggsave(filename = "/Users/caiyawei/Desktop/VAR/IRF_shock1.png", 
       plot = multiplot(p1,p2,p3, cols = 1),
       width = 30, height = 20, units = "cm",
       device = "png")

multiplot(p4,p5,p6, cols = 1)
ggsave(filename = "/Users/caiyawei/Desktop/VAR/IRF_shock2.png", 
       plot = multiplot(p4,p5,p6, cols = 1),
       width = 30, height = 20, units = "cm",
       device = "png")

multiplot(p7,p8,p9, cols = 1)
ggsave(filename = "/Users/caiyawei/Desktop/VAR/IRF_shock3.png", 
       plot = multiplot(p7,p8,p9, cols = 1),
       width = 30, height = 20, units = "cm",
       device = "png")

### Obtain the exogenous shock###
Coef.noc
Coef.OLS
const

# check the order of variables is correct (stock shock) (employment)
lm(stock_price ~ lag(property_price)+lag(stock_price)+lag(employment)+
     +lag(property_price, n=2)+lag(stock_price, n=2)+lag(employment, n=2)+lag(property_price, n=3)+lag(stock_price, n=3)+lag(employment, n=3)+lag(property_price, n=4)+lag(stock_price, n=4)+lag(employment, n=4),
   data = data)

get_residual <- function(vlb_order = 1){
  X.Mat <- By %>% as_tibble() %>%
    mutate(V1 = lag(property_price),
           V2 = lag(stock_price),
           V3 = lag(employment),
           V4 = lag(property_price, n = 2),
           V5 = lag(stock_price, n = 2),
           V6 = lag(employment, n = 2),
           V7 = lag(property_price, n = 3),
           V8 = lag(stock_price, n = 3),
           V9 = lag(employment, n = 3),
           V10 = lag(property_price, n = 4),
           V11 = lag(stock_price, n = 4),
           V12 = lag(employment, n = 4),
           V13 = rep(1, nrow(By))
    ) %>%
    select(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13
    ) %>%
    drop_na() %>%
    as.matrix
  Y.Mat <- By[(VAR.P+1):nrow(By),vlb_order] %>% as.matrix
  Y.Mat_hat <- X.Mat %*% (Coef.OLS[vlb_order, ] %>% as.matrix)
  residual <- Y.Mat - Y.Mat_hat
  return(residual)
}

# resid1 <- get_residual(1)
for(i in 1:3){
  assign(paste0("resid_",i), get_residual(i))
}

resid.Mat <- cbind(resid_1, resid_2, resid_3)
error.Mat <- (solve(B0) %*% A0 %*% t(resid.Mat)) %>% t()
error.Mat[,1] %>% hist

# Extract the stock shock

error.Mat %>% as_tibble() %>%
  select(e1 = 1, e2 = 2, e3 = 3) %>%
  mutate(no = row_number()) %>%
  ggplot()+
  geom_line(aes(no, e1))

data.err <- data %>% 
  bind_cols(e1 = c(rep(NA, VAR.P), error.Mat[,1])) %>%
  mutate(e1_pos = fifelse(e1 > 0, e1, 0),
         e1_neg = fifelse(e1 < 0, e1, 0))
data.err %>% head(10)

data.err %>% ggplot()+
  geom_line(aes(year, e1))

# check the order of variables is correct (house shock) (employment)
lm(property_price ~ lag(property_price)+lag(stock_price)+lag(employment)+
     +lag(property_price, n=2)+lag(stock_price, n=2)+lag(employment, n=2)+lag(property_price, n=3)+lag(stock_price, n=3)+lag(employment, n=3)+lag(property_price, n=4)+lag(stock_price, n=4)+lag(employment, n=4),
   data = data)

get_residual <- function(vlb_order = 1){
  X.Mat <- By %>% as_tibble() %>%
    mutate(V1 = lag(property_price),
           V2 = lag(stock_price),
           V3 = lag(employment),
           V4 = lag(property_price, n = 2),
           V5 = lag(stock_price, n = 2),
           V6 = lag(employment, n = 2),
           V7 = lag(property_price, n = 3),
           V8 = lag(stock_price, n = 3),
           V9 = lag(employment, n = 3),
           V10 = lag(property_price, n = 4),
           V11 = lag(stock_price, n = 4),
           V12 = lag(employment, n = 4),
           V13 = rep(1, nrow(By))
    ) %>%
    select(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13
    ) %>%
    drop_na() %>%
    as.matrix
  Y.Mat <- By[(VAR.P+1):nrow(By),vlb_order] %>% as.matrix
  Y.Mat_hat <- X.Mat %*% (Coef.OLS[vlb_order, ] %>% as.matrix)
  residual <- Y.Mat - Y.Mat_hat
  return(residual)
}

# resid1 <- get_residual(1)
for(i in 1:3){
  assign(paste0("resid_",i), get_residual(i))
}

resid.Mat <- cbind(resid_1, resid_2, resid_3)
error.Mat <- (solve(B0) %*% A0 %*% t(resid.Mat)) %>% t()
error.Mat[,2] %>% hist

# Extract the house shock

error.Mat %>% as_tibble() %>%
  select(e1 = 1, e2 = 2, e3 = 3) %>%
  mutate(no = row_number()) %>%
  ggplot()+
  geom_line(aes(no, e2))

data.err2 <- data %>% 
  bind_cols(e2 = c(rep(NA, VAR.P), error.Mat[,2])) %>%
  mutate(e2_pos = fifelse(e2 > 0, e2, 0),
         e2_neg = fifelse(e2 < 0, e2, 0))
data.err2 %>% head(10)

data.err2 %>% ggplot()+
  geom_line(aes(year, e2))


data.err["e2"] = data.err2$e2
### Regression ###
#reg1 <- lm(em_d ~ Period + e1 + e2 + interest_rate_d, data = data.err)
#summary(reg1)

### remove 1994Q4###
data.err["e2"] = data.err["e2"] - data.err$e2[6]
data.err["e1"] = data.err["e1"] - data.err$e1[6]
reg1 <- lm(em_d ~ Period + e1 + e2 + interest_rate_d, data = data.err)
summary(reg1)




### DATA CPI ###
#file = "/Users/caiyawei/Desktop/Japan Data/paper_data_quart.xlsx"
file = "/Users/caiyawei/Desktop/Japan Data/SVAR_forR.xlsx"
data = read_excel(file)
data = data[-c(107),]
data_cpi = data

By <- data_cpi %>% 
  select(property_price, stock_price, cpi) %>% 
  as.matrix

### RAW DATA ###
dim(By)
kk <- dim(By)[2]
draw_raw <- function(vlb_name = ""){
  vlb_name <- sym(vlb_name)
  plot.return <- data_cpi %>%
    ggplot()+
    geom_line(aes(year, !!vlb_name))
  return(plot.return)
}

p1 <- draw_raw("property_price")
p2 <- draw_raw("stock_price")
p3 <- draw_raw("cpi")

multiplot(p1, p2, p3)

### Setting Model ###
VAR.P = 4                       
CONST = TRUE                    
Y     = VAR.Y(By, VAR.P)       
X     = VAR.X(By, VAR.P)        

hrz=39 # the length of response

### Reduced Form VAR ###
(Coef.OLS    = VAR.OLS(Y, X, CONST)                  )
(Sigma.OLS   = VAR.Sigma.OLS(Y, X, Coef.OLS, CONST)  )
(Sigma.MLE   = VAR.Sigma.MLE(Y, X, Coef.OLS, CONST))

VAR.P = 4 # By AIC

Amat = diag(3)
# Identification Conditions

Bmat = diag(3)


Amat;Bmat

### Estimate A^hat and B^hat ###
#C.Prime <- chol(Sigma.OLS)
#C <- t(C.Prime)
#C

### Solving system of linear equations ###
B0 <-Bmat


A0 <-Amat
SVAR_AB_est <- list("A0.svar" = A0, "B0.svar" = B0)

### IRF Without Bootstrap CI ###
### IRF
SVAR_AB_IRF <- VAR.svarirf.AB(By, VAR.P, Amat, Bmat, h = hrz, CONST, SVAR_AB_est = SVAR_AB_est)
kk*1:kk

#output entire table
IRF_TABLE <- df_IRF_plot[,kk*1:kk] %>% select(property_price = 1,
                                              stock_price = 2,
                                              cpi = 3)
write.table(IRF_TABLE, file = "/Users/caiyawei/Desktop/IRF_TABLE3_order2.csv", sep = ",", row.names = FALSE)

draw_IRF <- function(df = df_IRF_plot, V1 = 1){
  V1 <- paste0("V", V1) %>% sym()
  plot.return <- ggplot(df) + 
    geom_line(aes(x = 1:nrow(df), y = !!V1))
  return(plot.return)
}

for(i in 1:(kk^2)){
  assign(paste0("p",i), draw_IRF(df_IRF_plot, i))
}

multiplot(p1,p2,p3,p4,p5,p6,
          p7,p8,p9,cols = 3)

### IRF (Bootstrap C.I.) ###
lower = 0.025                                        #  95% CI
upper = 1-lower
kk = ncol(By)
ddY = VAR.ddY(By, VAR.P)
ddX = VAR.ddX(By, VAR.P)

# dim(ddY); dim(ddX)

T   = nrow(ddY)
T.total= nrow(By)
Ik  = diag(rep(1, kk))
# 16 coef if 4 variables; 55 coef if 5 variables
Coef = t(VAR.EbyE(ddY, ddX, CONST)$ddA)              
# residuals
U    = VAR.EbyE(ddY, ddX, CONST)$ddU
BSigma.u = VAR.ddSigma.OLS(ddY, ddX, CONST)
if(CONST == TRUE){
  const = Coef[, ncol(Coef)]
  Coef.noc= Coef[,-ncol(Coef)]                      # const
}else{
  const = matrix(0, kk, 1)
  Coef.noc = Coef
}

Theta.unit= VAR.Theta(Coef, h, BSigma.u, CONST)$unit #  Theta.unit
Theta.std = VAR.Theta(Coef, h, BSigma.u, CONST)$std  #  Theta.std

# dm.U <- U-mean(U)
dm.U <- U

N = 2000 
Theta.unit.sim = vector("list", N)
Theta.std.sim  = vector("list", N)

# check dimension
print("check dimensionality")
dim(ddX); dim(Coef.noc); dim(dm.U)


df_IRF.sim <- array(NA, c(hrz+1,kk^2,N)) #dimensions are: Time Period, Number of shock interacts with variables, page (number of Bootstrap resamplings)
counter <- 1
while(TRUE){
  
  #cat("Now, there are ", counter-1, " sets of resamples.\n")
  Y.sim = matrix(0, nrow = T.total, ncol = kk)          # Y.sim = 0 #pseudo time series
  Y.sim[c(1:VAR.P),] = By[c(1:VAR.P), ] #initial values
  
  boot.number = sample(c(1:T), replace = TRUE)      
  U.sim = dm.U[boot.number,]
  
  # predicted values given the above initial values
  last.y= c(t(By[VAR.P:1,]))
  for(ii in 1:T){
    last.y = last.y[1:(kk*VAR.P)]
    Y.sim[ii+VAR.P, ] = Coef.noc %*% last.y + const + U.sim[ii,]      
    last.y = c(Y.sim[ii+VAR.P,], last.y)
  }

  
  ### SVAR.sim Start ###
  
  Y_pseudo     = VAR.Y(Y.sim, VAR.P)        
  X_pseudo     = VAR.X(Y.sim, VAR.P)        
  Coef.OLS_pseudo    = VAR.OLS(Y_pseudo, X_pseudo, CONST)
  Sigma.OLS_pseudo   = VAR.Sigma.OLS(Y_pseudo, X_pseudo, Coef.OLS_pseudo, CONST)
  C.Prime_pseudo <- chol(Sigma.OLS_pseudo)
  C_pseudo <- t(C.Prime_pseudo)
  B0_pseudo <- diag(diag(C_pseudo), ncol = kk, nrow = kk)
  A0_pseudo <- B0_pseudo %*% solve(C_pseudo)
  SVAR_AB_est.sim <- list("A0.svar" = A0_pseudo, "B0.svar" = B0_pseudo)
  SVAR_AB_IRF.sim <- VAR.svarirf.AB(Y.sim, VAR.P, Amat, Bmat, h = hrz, CONST, SVAR_AB_est = SVAR_AB_est.sim)
  
  df_IRF_plot.sim <- matrix(NA, hrz+1, kk^2) #%>% as.tibble()
  # df_IRF.sim <- array(1:(120*25*N), c(120,25,N))
  # df_IRF.sim[2,1,1] # slicing
  
  h <- 0 
  for(period in SVAR_AB_IRF.sim){
    k <- 0 
    h <- h+1 
    for(j in 1:kk){
      for(i in 1:kk){
        k <- k+1 
        df_IRF_plot.sim[h,k] <- period[i,j]
      }
    }
  }
  df_IRF.sim[,,counter] <- df_IRF_plot.sim
  ### SVAR.sim Ends ###
  if(counter>=N){
    break
  }
  counter <- counter+1
}
# Save
saveRDS(df_IRF.sim, file = "df_IRF.sim.rds")

df_IRF.sim <- read_rds("df_IRF.sim.rds")

head(df_IRF.sim[,,1000])
print(sum(is.na(df_IRF.sim)))

# IRF & Bootstrap C.I.
df_IRF_plot.BS.L <- matrix(NA, nrow = hrz+1, ncol = kk^2)
df_IRF_plot.BS.U <- matrix(NA, nrow = hrz+1, ncol = kk^2)
df_IRF_plot.BS.Median <- matrix(NA, nrow = hrz+1, ncol = kk^2)
df_IRF_plot.BS.Mean <- matrix(NA, nrow = hrz+1, ncol = kk^2)
for(col in 1:(kk^2)){
  for(row in 1:(hrz+1) ){
    df_IRF_plot.BS.L[row,col] <- quantile(df_IRF.sim[row,col,], probs = 0.025)
    df_IRF_plot.BS.U[row,col] <- quantile(df_IRF.sim[row,col,], probs = 0.975)
    df_IRF_plot.BS.Median[row,col] <- quantile(df_IRF.sim[row,col,], probs = 0.5)
    df_IRF_plot.BS.Mean[row,col] <- mean(df_IRF.sim[row,col,])
  }
}

df_IRF_plot.BS.L <- df_IRF_plot.BS.L %>% as_tibble()
df_IRF_plot.BS.U <- df_IRF_plot.BS.U %>% as_tibble()
df_IRF_plot.BS.Median <- df_IRF_plot.BS.Median %>% as_tibble()
df_IRF_plot.BS.Mean <- df_IRF_plot.BS.Mean %>% as_tibble()

ind <- 0
for(i in 1:kk){
  for(j in 1:kk){
    ind <- ind+1
    nam <- paste("shock", j, "y", i, sep = '')
    assign(nam, bind_cols(df_IRF_plot.BS.L[ind], df_IRF_plot.BS.U[ind],
                          df_IRF_plot.BS.Median[ind], df_IRF_plot.BS.Mean[ind],
                          df_IRF_plot[ind]))
    evalStr <- paste0("colnames(", nam, ") <- c('Lower', 'Upper', 'Median', 'Mean', 'Actual')")
    eval(parse(text=evalStr))
    evalStr <- paste0("p", ind, " <- ", "ggplot(",nam,") +geom_hline(yintercept=0, color = 'grey')+ geom_line(aes(x = 1:nrow(", nam, "), y = Lower), linetype = 'dashed', col='red')+geom_line(aes(x = 1:nrow(", nam, "), y = Upper), linetype = 'dashed', col='red')+geom_line(aes(x = 1:nrow(", nam, "), y = Median), col = 'Blue')")
    eval(parse(text=evalStr))
  }
}  

Text_Size_Theme = theme(
  axis.title.x = element_text(size = 12),
  axis.text.x = element_text(size = 12),
  axis.title.y = element_text(size = 12),
  axis.title = element_text(size = 12),
  plot.title = element_text(size=12))

add_layout <- function(p = p1,
                       response_of = "this variable", react_to = "the shock"){
  title_text <- paste("Response of", response_of, "to", react_to, sep = ' ')
  plot.return <- p + labs(x = "Period", y = "", title = title_text) + Text_Size_Theme
  return(plot.return)
}

vlb_list = c("House Price",
             "Stock Price",
             "CPI"
)
p1 <- add_layout(p1, vlb_list[1], "Property Price shock")
p2 <- add_layout(p2, vlb_list[2], "Property Price shock")
p3 <- add_layout(p3, vlb_list[3], "Property Price shock")

p4 <- add_layout(p4, vlb_list[1], "Stock shock")
p5 <- add_layout(p5, vlb_list[2], "Stock shock")
p6 <- add_layout(p6, vlb_list[3], "Stock shock")

p7 <- add_layout(p7, vlb_list[1], "CPI shock")
p8 <- add_layout(p8, vlb_list[2], "CPI shock")
p9 <- add_layout(p9, vlb_list[3], "CPI shock")

multiplot(p1,p2,p3,p4,p5,p6,
          p7,p8,p9,cols = 3)

multiplot(p1,p2,p3, cols = 1)
ggsave(filename = "/Users/caiyawei/Desktop/VAR/IRF_shock_cpi1.png", 
       plot = multiplot(p1,p2,p3, cols = 1),
       width = 30, height = 20, units = "cm",
       device = "png")

multiplot(p4,p5,p6, cols = 1)
ggsave(filename = "/Users/caiyawei/Desktop/VAR/IRF_shock_cpi2.png", 
       plot = multiplot(p4,p5,p6, cols = 1),
       width = 30, height = 20, units = "cm",
       device = "png")

multiplot(p7,p8,p9, cols = 1)
ggsave(filename = "/Users/caiyawei/Desktop/VAR/IRF_shock_cpi3.png", 
       plot = multiplot(p7,p8,p9, cols = 1),
       width = 30, height = 20, units = "cm",
       device = "png")

### Obtain the exogenous shock###
Coef.noc
Coef.OLS
const

# check the order of variables is correct (stock shock) (cpi)
lm(stock_price ~ lag(property_price)+lag(stock_price)+lag(cpi)+
     +lag(property_price, n=2)+lag(stock_price, n=2)+lag(cpi, n=2)+lag(property_price, n=3)+lag(stock_price, n=3)+lag(cpi, n=3)+lag(property_price, n=4)+lag(stock_price, n=4)+lag(cpi, n=4),
   data = data_cpi)

get_residual <- function(vlb_order = 1){
  X.Mat <- By %>% as_tibble() %>%
    mutate(V1 = lag(property_price),
           V2 = lag(stock_price),
           V3 = lag(cpi),
           V4 = lag(property_price, n = 2),
           V5 = lag(stock_price, n = 2),
           V6 = lag(cpi, n = 2),
           V7 = lag(property_price, n = 3),
           V8 = lag(stock_price, n = 3),
           V9 = lag(cpi, n = 3),
           V10 = lag(property_price, n = 4),
           V11 = lag(stock_price, n = 4),
           V12 = lag(cpi, n = 4),
           V13 = rep(1, nrow(By))
    ) %>%
    select(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13
    ) %>%
    drop_na() %>%
    as.matrix
  Y.Mat <- By[(VAR.P+1):nrow(By),vlb_order] %>% as.matrix
  Y.Mat_hat <- X.Mat %*% (Coef.OLS[vlb_order, ] %>% as.matrix)
  residual <- Y.Mat - Y.Mat_hat
  return(residual)
}

# resid1 <- get_residual(1)
for(i in 1:3){
  assign(paste0("resid_",i), get_residual(i))
}

resid.Mat <- cbind(resid_1, resid_2, resid_3)
error.Mat <-  (solve(B0) %*% A0 %*% t(resid.Mat)) %>% t()
error.Mat[,1] %>% hist

# Extract the stock shock

error.Mat %>% as_tibble() %>%
  select(e1 = 1, e2 = 2, e3 = 3) %>%
  mutate(no = row_number()) %>%
  ggplot()+
  geom_line(aes(no, e1))

data.err <- data_cpi %>% 
  bind_cols(e1 = c(rep(NA, VAR.P), error.Mat[,1])) %>%
  mutate(e1_pos = fifelse(e1 > 0, e1, 0),
         e1_neg = fifelse(e1 < 0, e1, 0))
data.err %>% head(10)

data.err %>% ggplot()+
  geom_line(aes(year, e1))

# check the order of variables is correct (house shock) (cpi)
lm(property_price ~ lag(property_price)+lag(stock_price)+lag(cpi)+
     +lag(property_price, n=2)+lag(stock_price, n=2)+lag(cpi, n=2)+lag(property_price, n=3)+lag(stock_price, n=3)+lag(cpi, n=3)+lag(property_price, n=4)+lag(stock_price, n=4)+lag(cpi, n=4),
   data = data_cpi)

get_residual <- function(vlb_order = 1){
  X.Mat <- By %>% as_tibble() %>%
    mutate(V1 = lag(property_price),
           V2 = lag(stock_price),
           V3 = lag(cpi),
           V4 = lag(property_price, n = 2),
           V5 = lag(stock_price, n = 2),
           V6 = lag(cpi, n = 2),
           V7 = lag(property_price, n = 3),
           V8 = lag(stock_price, n = 3),
           V9 = lag(cpi, n = 3),
           V10 = lag(property_price, n = 4),
           V11 = lag(stock_price, n = 4),
           V12 = lag(cpi, n = 4),
           V13 = rep(1, nrow(By))
    ) %>%
    select(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13
    ) %>%
    drop_na() %>%
    as.matrix
  Y.Mat <- By[(VAR.P+1):nrow(By),vlb_order] %>% as.matrix
  Y.Mat_hat <- X.Mat %*% (Coef.OLS[vlb_order, ] %>% as.matrix)
  residual <- Y.Mat - Y.Mat_hat
  return(residual)
}

# resid1 <- get_residual(1)
for(i in 1:3){
  assign(paste0("resid_",i), get_residual(i))
}

resid.Mat <- cbind(resid_1, resid_2, resid_3)
error.Mat <-  (solve(B0) %*% A0 %*% t(resid.Mat)) %>% t()
error.Mat[,2] %>% hist

# Extract the house shock

error.Mat %>% as_tibble() %>%
  select(e1 = 1, e2 = 2, e3 = 3) %>%
  mutate(no = row_number()) %>%
  ggplot()+
  geom_line(aes(no, e2))

data.err2 <- data_cpi %>% 
  bind_cols(e2 = c(rep(NA, VAR.P), error.Mat[,2])) %>%
  mutate(e2_pos = fifelse(e2 > 0, e2, 0),
         e2_neg = fifelse(e2 < 0, e2, 0))
data.err2 %>% head(10)

data.err2 %>% ggplot()+
  geom_line(aes(year, e2))


data.err["e2"] = data.err2$e2
### Regression ###
#reg3 <- lm(cpi_d ~ Period + e1 + e2 + interest_rate_d, data = data.err)
#summary(reg3)
###衝擊減掉1991Q4###
data.err["e2"] = data.err["e2"] - data.err$e2[6]
data.err["e1"] = data.err["e1"] - data.err$e1[6]
reg3 <- lm(cpi_d ~ Period + e1 + e2 + interest_rate_d, data = data.err)
summary(reg3)







### DATA GDP ###
data = read_excel(file)
data = data[-c(107),]
data_gdp = data

By <- data_gdp %>% 
  select(property_price, stock_price, gdp) %>% 
  as.matrix

### RAW DATA ###
dim(By)
kk <- dim(By)[2]
draw_raw <- function(vlb_name = ""){
  vlb_name <- sym(vlb_name)
  plot.return <- data_gdp %>%
    ggplot()+
    geom_line(aes(year, !!vlb_name))
  return(plot.return)
}

p1 <- draw_raw("property_price")
p2 <- draw_raw("stock_price")
p3 <- draw_raw("gdp")

multiplot(p1, p2, p3)

### Setting Model ###
VAR.P = 4                       
CONST = TRUE                    
Y     = VAR.Y(By, VAR.P)        
X     = VAR.X(By, VAR.P)        

hrz=39 # the length of response

### Reduced Form VAR ###
(Coef.OLS    = VAR.OLS(Y, X, CONST)                  )
(Sigma.OLS   = VAR.Sigma.OLS(Y, X, Coef.OLS, CONST)  )
(Sigma.MLE   = VAR.Sigma.MLE(Y, X, Coef.OLS, CONST))

VAR.P = 4 # By AIC

Amat = diag(3)
# Identification Conditions

Bmat = diag(3)


Amat;Bmat

### Estimate A^hat and B^hat ###
#C.Prime <- chol(Sigma.OLS)
#C <- t(C.Prime)
#C

### Solving system of linear equations ###
B0 <-Bmat


A0 <-Amat

SVAR_AB_est <- list("A0.svar" = A0, "B0.svar" = B0)

### IRF Without Bootstrap CI ###
### IRF
SVAR_AB_IRF <- VAR.svarirf.AB(By, VAR.P, Amat, Bmat, h = hrz, CONST, SVAR_AB_est = SVAR_AB_est)

kk*1:kk

#output entire table
IRF_TABLE <- df_IRF_plot[,kk*1:kk] %>% select(property_price = 1,
                                              stock_price = 2,
                                              gdp = 3)
write.table(IRF_TABLE, file = "/Users/caiyawei/Desktop/IRF_TABLE2.csv", sep = ",", row.names = FALSE)

draw_IRF <- function(df = df_IRF_plot, V1 = 1){
  V1 <- paste0("V", V1) %>% sym()
  plot.return <- ggplot(df) + 
    geom_line(aes(x = 1:nrow(df), y = !!V1))
  return(plot.return)
}

for(i in 1:(kk^2)){
  assign(paste0("p",i), draw_IRF(df_IRF_plot, i))
}

multiplot(p1,p2,p3,p4,p5,p6,
          p7,p8,p9,cols = 3)

### IRF (Bootstrap C.I.) ###
lower = 0.025                                        #  95% CI
upper = 1-lower
kk = ncol(By)
ddY = VAR.ddY(By, VAR.P)
ddX = VAR.ddX(By, VAR.P)


T   = nrow(ddY)
T.total= nrow(By)
Ik  = diag(rep(1, kk))
Coef = t(VAR.EbyE(ddY, ddX, CONST)$ddA)              
U    = VAR.EbyE(ddY, ddX, CONST)$ddU
BSigma.u = VAR.ddSigma.OLS(ddY, ddX, CONST)
if(CONST == TRUE){
  const = Coef[, ncol(Coef)]
  Coef.noc= Coef[,-ncol(Coef)]                    
}else{
  const = matrix(0, kk, 1)
  Coef.noc = Coef
}

Theta.unit= VAR.Theta(Coef, h, BSigma.u, CONST)$unit 
Theta.std = VAR.Theta(Coef, h, BSigma.u, CONST)$std  

dm.U <- U

N = 2000 
Theta.unit.sim = vector("list", N)
Theta.std.sim  = vector("list", N)

print("check dimensionality")
dim(ddX); dim(Coef.noc); dim(dm.U)


df_IRF.sim <- array(NA, c(hrz+1,kk^2,N)) 
counter <- 1
while(TRUE){
  
  Y.sim = matrix(0, nrow = T.total, ncol = kk)          
  Y.sim[c(1:VAR.P),] = By[c(1:VAR.P), ] 
  
  boot.number = sample(c(1:T), replace = TRUE)      
  U.sim = dm.U[boot.number,]
  
  last.y= c(t(By[VAR.P:1,]))
  for(ii in 1:T){
    last.y = last.y[1:(kk*VAR.P)]
    Y.sim[ii+VAR.P, ] = Coef.noc %*% last.y + const + U.sim[ii,]      
    last.y = c(Y.sim[ii+VAR.P,], last.y)
  }
  

  
  Y_pseudo     = VAR.Y(Y.sim, VAR.P)        
  X_pseudo     = VAR.X(Y.sim, VAR.P)        
  Coef.OLS_pseudo    = VAR.OLS(Y_pseudo, X_pseudo, CONST)
  Sigma.OLS_pseudo   = VAR.Sigma.OLS(Y_pseudo, X_pseudo, Coef.OLS_pseudo, CONST)
  C.Prime_pseudo <- chol(Sigma.OLS_pseudo)
  C_pseudo <- t(C.Prime_pseudo)
  B0_pseudo <- diag(diag(C_pseudo), ncol = kk, nrow = kk)
  A0_pseudo <- B0_pseudo %*% solve(C_pseudo)
  SVAR_AB_est.sim <- list("A0.svar" = A0_pseudo, "B0.svar" = B0_pseudo)
  SVAR_AB_IRF.sim <- VAR.svarirf.AB(Y.sim, VAR.P, Amat, Bmat, h = hrz, CONST, SVAR_AB_est = SVAR_AB_est.sim)
  
  df_IRF_plot.sim <- matrix(NA, hrz+1, kk^2) 

  
  h <- 0 
  for(period in SVAR_AB_IRF.sim){
    k <- 0
    h <- h+1 
    for(j in 1:kk){
      for(i in 1:kk){
        k <- k+1
        df_IRF_plot.sim[h,k] <- period[i,j]
      }
    }
  }
  df_IRF.sim[,,counter] <- df_IRF_plot.sim
  ### SVAR.sim Ends ###
  if(counter>=N){
    break
  }
  counter <- counter+1
}
saveRDS(df_IRF.sim, file = "df_IRF.sim.rds")

df_IRF.sim <- read_rds("df_IRF.sim.rds")

head(df_IRF.sim[,,1000])
print(sum(is.na(df_IRF.sim)))

# draw IRF & Bootstrap C.I.
df_IRF_plot.BS.L <- matrix(NA, nrow = hrz+1, ncol = kk^2)
df_IRF_plot.BS.U <- matrix(NA, nrow = hrz+1, ncol = kk^2)
df_IRF_plot.BS.Median <- matrix(NA, nrow = hrz+1, ncol = kk^2)
df_IRF_plot.BS.Mean <- matrix(NA, nrow = hrz+1, ncol = kk^2)
for(col in 1:(kk^2)){
  for(row in 1:(hrz+1) ){
    df_IRF_plot.BS.L[row,col] <- quantile(df_IRF.sim[row,col,], probs = 0.025)
    df_IRF_plot.BS.U[row,col] <- quantile(df_IRF.sim[row,col,], probs = 0.975)
    df_IRF_plot.BS.Median[row,col] <- quantile(df_IRF.sim[row,col,], probs = 0.5)
    df_IRF_plot.BS.Mean[row,col] <- mean(df_IRF.sim[row,col,])
  }
}

df_IRF_plot.BS.L <- df_IRF_plot.BS.L %>% as_tibble()
df_IRF_plot.BS.U <- df_IRF_plot.BS.U %>% as_tibble()
df_IRF_plot.BS.Median <- df_IRF_plot.BS.Median %>% as_tibble()
df_IRF_plot.BS.Mean <- df_IRF_plot.BS.Mean %>% as_tibble()

ind <- 0
for(i in 1:kk){
  for(j in 1:kk){
    ind <- ind+1
    nam <- paste("shock", j, "y", i, sep = '')
    assign(nam, bind_cols(df_IRF_plot.BS.L[ind], df_IRF_plot.BS.U[ind],
                          df_IRF_plot.BS.Median[ind], df_IRF_plot.BS.Mean[ind],
                          df_IRF_plot[ind]))
    evalStr <- paste0("colnames(", nam, ") <- c('Lower', 'Upper', 'Median', 'Mean', 'Actual')")
    eval(parse(text=evalStr))
    evalStr <- paste0("p", ind, " <- ", "ggplot(",nam,") +geom_hline(yintercept=0, color = 'grey')+ geom_line(aes(x = 1:nrow(", nam, "), y = Lower), linetype = 'dashed', col='red')+geom_line(aes(x = 1:nrow(", nam, "), y = Upper), linetype = 'dashed', col='red')+geom_line(aes(x = 1:nrow(", nam, "), y = Median), col = 'Blue')")
    eval(parse(text=evalStr))
  }
}  

Text_Size_Theme = theme(
  axis.title.x = element_text(size = 12),
  axis.text.x = element_text(size = 12),
  axis.title.y = element_text(size = 12),
  axis.title = element_text(size = 12),
  plot.title = element_text(size=12))

add_layout <- function(p = p1,
                       response_of = "this variable", react_to = "the shock"){
  title_text <- paste("Response of", response_of, "to", react_to, sep = ' ')
  plot.return <- p + labs(x = "Period", y = "", title = title_text) + Text_Size_Theme
  return(plot.return)
}

vlb_list = c("House Price",
             "Stock Price",
             "GDP"
)
p1 <- add_layout(p1, vlb_list[1], "Property Price shock")
p2 <- add_layout(p2, vlb_list[2], "Property Price shock")
p3 <- add_layout(p3, vlb_list[3], "Property Price shock")

p4 <- add_layout(p4, vlb_list[1], "Stock shock")
p5 <- add_layout(p5, vlb_list[2], "Stock shock")
p6 <- add_layout(p6, vlb_list[3], "Stock shock")

p7 <- add_layout(p7, vlb_list[1], "GDP shock")
p8 <- add_layout(p8, vlb_list[2], "GDP shock")
p9 <- add_layout(p9, vlb_list[3], "GDP shock")

multiplot(p1,p2,p3,p4,p5,p6,
          p7,p8,p9,cols = 3)

multiplot(p1,p2,p3, cols = 1)
ggsave(filename = "/Users/caiyawei/Desktop/VAR/IRF_shock_gdp1.png", 
       plot = multiplot(p1,p2,p3, cols = 1),
       width = 30, height = 20, units = "cm",
       device = "png")

multiplot(p4,p5,p6, cols = 1)
ggsave(filename = "/Users/caiyawei/Desktop/VAR/IRF_shock_gdp2.png", 
       plot = multiplot(p4,p5,p6, cols = 1),
       width = 30, height = 20, units = "cm",
       device = "png")

multiplot(p7,p8,p9, cols = 1)
ggsave(filename = "/Users/caiyawei/Desktop/VAR/IRF_shock_gdp33.png", 
       plot = multiplot(p7,p8,p9, cols = 1),
       width = 30, height = 20, units = "cm",
       device = "png")

### Obtain the exogenous shock###
Coef.noc
Coef.OLS
const

# check the order of variables is correct (stock shock) (gdp)
lm(stock_price ~ lag(property_price)+lag(stock_price)+lag(gdp)+
     +lag(property_price, n=2)+lag(stock_price, n=2)+lag(gdp, n=2)+lag(property_price, n=3)+lag(stock_price, n=3)+lag(gdp, n=3)+lag(property_price, n=4)+lag(stock_price, n=4)+lag(gdp, n=4),
   data = data_gdp)

get_residual <- function(vlb_order = 1){
  X.Mat <- By %>% as_tibble() %>%
    mutate(V1 = lag(property_price),
           V2 = lag(stock_price),
           V3 = lag(gdp),
           V4 = lag(property_price, n = 2),
           V5 = lag(stock_price, n = 2),
           V6 = lag(gdp, n = 2),
           V7 = lag(property_price, n = 3),
           V8 = lag(stock_price, n = 3),
           V9 = lag(gdp, n = 3),
           V10 = lag(property_price, n = 4),
           V11 = lag(stock_price, n = 4),
           V12 = lag(gdp, n = 4),
           V13 = rep(1, nrow(By))
    ) %>%
    select(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13
    ) %>%
    drop_na() %>%
    as.matrix
  Y.Mat <- By[(VAR.P+1):nrow(By),vlb_order] %>% as.matrix
  Y.Mat_hat <- X.Mat %*% (Coef.OLS[vlb_order, ] %>% as.matrix)
  residual <- Y.Mat - Y.Mat_hat
  return(residual)
}

# resid1 <- get_residual(1)
for(i in 1:3){
  assign(paste0("resid_",i), get_residual(i))
}

resid.Mat <- cbind(resid_1, resid_2, resid_3)
error.Mat <-  (solve(B0) %*% A0 %*% t(resid.Mat)) %>% t()
error.Mat[,1] %>% hist

# Extract the stock shock

error.Mat %>% as_tibble() %>%
  select(e1 = 1, e2 = 2, e3 = 3) %>%
  mutate(no = row_number()) %>%
  ggplot()+
  geom_line(aes(no, e1))

data.err <- data_gdp %>% 
  bind_cols(e1 = c(rep(NA, VAR.P), error.Mat[,1])) %>%
  mutate(e1_pos = fifelse(e1 > 0, e1, 0),
         e1_neg = fifelse(e1 < 0, e1, 0))
data.err %>% head(10)

data.err %>% ggplot()+
  geom_line(aes(year, e1))

# check the order of variables is correct (house shock) (gdp)
lm(property_price ~ lag(property_price)+lag(stock_price)+lag(gdp)+
     +lag(property_price, n=2)+lag(stock_price, n=2)+lag(gdp, n=2)+lag(property_price, n=3)+lag(stock_price, n=3)+lag(gdp, n=3)+lag(property_price, n=4)+lag(stock_price, n=4)+lag(gdp, n=4),
   data = data_gdp)

get_residual <- function(vlb_order = 1){
  X.Mat <- By %>% as_tibble() %>%
    mutate(V1 = lag(property_price),
           V2 = lag(stock_price),
           V3 = lag(gdp),
           V4 = lag(property_price, n = 2),
           V5 = lag(stock_price, n = 2),
           V6 = lag(gdp, n = 2),
           V7 = lag(property_price, n = 3),
           V8 = lag(stock_price, n = 3),
           V9 = lag(gdp, n = 3),
           V10 = lag(property_price, n = 4),
           V11 = lag(stock_price, n = 4),
           V12 = lag(gdp, n = 4),
           V13 = rep(1, nrow(By))
    ) %>%
    select(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13
    ) %>%
    drop_na() %>%
    as.matrix
  Y.Mat <- By[(VAR.P+1):nrow(By),vlb_order] %>% as.matrix
  Y.Mat_hat <- X.Mat %*% (Coef.OLS[vlb_order, ] %>% as.matrix)
  residual <- Y.Mat - Y.Mat_hat
  return(residual)
}

# resid1 <- get_residual(1)
for(i in 1:3){
  assign(paste0("resid_",i), get_residual(i))
}

resid.Mat <- cbind(resid_1, resid_2, resid_3)
error.Mat <-  (solve(B0) %*% A0 %*% t(resid.Mat)) %>% t()
error.Mat[,2] %>% hist

# Extract the house shock

error.Mat %>% as_tibble() %>%
  select(e1 = 1, e2 = 2, e3 = 3) %>%
  mutate(no = row_number()) %>%
  ggplot()+
  geom_line(aes(no, e2))

data.err2 <- data_gdp %>% 
  bind_cols(e2 = c(rep(NA, VAR.P), error.Mat[,2])) %>%
  mutate(e2_pos = fifelse(e2 > 0, e2, 0),
         e2_neg = fifelse(e2 < 0, e2, 0))
data.err2 %>% head(10)

data.err2 %>% ggplot()+
  geom_line(aes(year...1, e2))


data.err["e2"] = data.err2$e2

### Regression ###
#reg2 <- lm(gdp_d ~ Period + e1 + e2 + interest_rate_d, data = data.err)
#summary(reg2)

data.err["e2"] = data.err["e2"] - data.err$e2[6]
data.err["e1"] = data.err["e1"] - data.err$e1[6]
reg2 <- lm(gdp_d ~ Period + e1 + e2 + interest_rate_d, data = data.err)
summary(reg2)

texreg(list(reg1, reg2, reg3), caption = "Three VAR models in Order2.", digits = 5)
















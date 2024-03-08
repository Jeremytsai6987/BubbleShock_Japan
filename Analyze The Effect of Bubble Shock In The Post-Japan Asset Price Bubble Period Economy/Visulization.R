library(ggplot2)
library(dplyr)
library(ggthemes)
library(e1071) 
library(lubridate)
library(tidyr)
theme_economist(
  base_size = 10,
  base_family = "sans",
  horizontal = TRUE,
  dkpanel = FALSE
)

theme_economist_white(
  base_size = 11,
  base_family = "sans",
  gray_bg = TRUE,
  horizontal = TRUE
)

#1.Real residential property price
real_res_property_price <- read.csv("/Users/caiyawei/Desktop/new data/data_for_last_part.csv")
real_res_property_price <- as.data.frame(real_res_property_price)
real_res_property_price$DATE <- as.Date(real_res_property_price$Year)

p <- ggplot(real_res_property_price, aes(x = DATE, y = Property_price))+
  ggtitle("Plot of Real Residential Property Price") +
  xlab("Date") +
  ylab("Real Residential Property Price") +
  theme_economist_white()+
  geom_line()+
  theme_economist()
p

mean(real_res_property_price$Property_price) 
sd(real_res_property_price$Property_price)
skewness(real_res_property_price$Property_price)
kurtosis(real_res_property_price$Property_price)

#2.employment Rate
umemployment <- read.csv("/Users/caiyawei/Desktop/new data/data_for_last_part.csv")
umemployment <- as.data.frame(umemployment)
umemployment$DATE <- as.Date(umemployment$Year)

p1 <- ggplot(umemployment, aes(x = DATE, y = Employment))+
  ggtitle("Plot of employment Rate") +
  xlab("Date") +
  ylab("employment Rate") +
  theme_economist_white()+
  geom_line()+
  theme_economist()
p1

mean(umemployment$Employment) 
sd(umemployment$Employment)
skewness(umemployment$Employment)
kurtosis(umemployment$Employment)

#3.Employment
employment <- read.csv("/Users/caiyawei/Desktop/Japan Data/Employment Rate All Persons for Japan /LREM64TTJPM156S.csv")
employment <- as.data.frame(employment)
employment$DATE <- as.Date(employment$DATE)

p2 <- ggplot(employment, aes(x = DATE, y = LREM64TTJPM156S))+
  ggtitle("Plot of Employment Rate All Person For Japan") +
  xlab("Date") +
  ylab("Employment Rate") +
  theme_economist_white()+
  geom_rect(data = employment,
            aes(xmin = as.Date("1986-12-01"), xmax = as.Date("1991-03-31"), ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "#99ccff", alpha = 0.2)+ 
  geom_line()+
  theme_economist()
p2

mean(employment$LREM64TTJPM156S) 
sd(employment$LREM64TTJPM156S)
skewness(employment$LREM64TTJPM156S)
kurtosis(employment$LREM64TTJPM156S)

#4.Comsumer Price of All Items for Japan
cpi <- read.csv("/Users/caiyawei/Desktop/new data/data_for_last_part.csv")
cpi <- as.data.frame(cpi)
cpi$DATE <- as.Date(cpi$Year)

p3 <- ggplot(cpi, aes(x = DATE, y = CPI))+
  ggtitle("Plot of Consumer Price Index") +
  xlab("Date") +
  ylab("Consumer Price Index ") +
  theme_economist_white()+
  geom_line()+
  theme_economist()
p3

mean(cpi$CPI) 
sd(cpi$CPI)
skewness(cpi$CPI)
kurtosis(cpi$CPI)

#5.Long-Term Government Bond Yields
gov_bond <- read.csv("/Users/caiyawei/Desktop/Japan Data/Long-Term Government Bond Yields/IRLTLT01JPM156N.csv")
gov_bond <- as.data.frame(gov_bond)
gov_bond$DATE <- as.Date(gov_bond$DATE)

p4 <- ggplot(gov_bond, aes(x = DATE, y = IRLTLT01JPM156N))+
  ggtitle("Plot of Long-Term Government Bond Yields In Japan ") +
  xlab("Date") +
  ylab("Long-Term Government Bond Yields ") +
  theme_economist_white()+
  geom_rect(data = gov_bond,
            aes(xmin = as.Date("1986-12-01"), xmax = as.Date("1991-03-31"), ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "#99ccff", alpha = 0.2)+ 
  geom_line()+
  theme_economist()
p4

mean(gov_bond$IRLTLT01JPM156N) 
sd(gov_bond$IRLTLT01JPM156N)
skewness(gov_bond$IRLTLT01JPM156N)
kurtosis(gov_bond$IRLTLT01JPM156N)

#6.Interest Rates, Discount Rate for Japan
interest <- read.csv("/Users/caiyawei/Desktop/new data/data_for_last_part.csv")
interest <- as.data.frame(interest)
interest$DATE <- as.Date(interest$Year)

p5 <- ggplot(interest, aes(x = DATE, y = interest))+
  ggtitle("Plot of Interest Rates") +
  xlab("Date") +
  ylab("Interest Rates")+
  theme_economist_white()+
  geom_line()+
  theme_economist()
p5

mean(interest$interest) 
sd(interest$interest)
skewness(interest$interest)
kurtosis(interest$interest)

#7./Gross National Income for Japan
gni <- read.csv("/Users/caiyawei/Desktop/Japan Data/Gross National Income for Japan/MKTGNIJPA646NWDB.csv")
gni <- as.data.frame(gni)
gni$DATE <- as.Date(gni$DATE)

p6 <- ggplot(gni, aes(x = DATE, y = MKTGNIJPA646NWDB))+
  ggtitle("Plot of Gross National Income for Japan ") +
  xlab("Date") +
  ylab("Gross National Income for Japan") +
  theme_economist_white()+
  geom_rect(data = gni,
            aes(xmin = as.Date("1986-12-01"), xmax = as.Date("1991-03-31"), ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "#99ccff", alpha = 0.2)+ 
  geom_line()+
  theme_economist()
p6

mean(gni$MKTGNIJPA646NWDB) 
sd(gni$MKTGNIJPA646NWDB)
skewness(gni$MKTGNIJPA646NWDB)
kurtosis(gni$MKTGNIJPA646NWDB)

#8.Credit to Households and NPISHs, Adjusted for Breaks, for Japan
credit2hous <- read.csv("//Users/caiyawei/Desktop/Japan Data/Total Credit to Households and NPISHs, Adjusted for Breaks, for Japan/QJPHAMUSDA.csv")
credit2hous <- as.data.frame(credit2hous)
credit2hous$DATE <- as.Date(credit2hous$DATE)

p7 <- ggplot(credit2hous, aes(x = DATE, y = QJPHAMUSDA))+
  ggtitle("Plot of Credit to Households and NPISHs ") +
  xlab("Date") +
  ylab("Credit to Households and NPISHs, Adjusted for Breaks, for Japan") +
  theme_economist_white()+
  geom_rect(data = credit2hous,
            aes(xmin = as.Date("1986-12-01"), xmax = as.Date("1991-03-31"), ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "#99ccff", alpha = 0.2)+ 
  geom_line()+
  theme_economist()
p7

mean(credit2hous$QJPHAMUSDA) 
sd(credit2hous$QJPHAMUSDA)
skewness(credit2hous$QJPHAMUSDA)
kurtosis(credit2hous$QJPHAMUSDA)

#9.Producer Prices Index Total Durable Consumer Goods for Japan 
ppi_durable <- read.csv("/Users/caiyawei/Desktop/Japan Data/Producer Prices Index Total Durable Consumer Goods for Japan /PITGCD01JPM661N.csv")
ppi_durable <- as.data.frame(ppi_durable)
ppi_durable$DATE <- as.Date(ppi_durable$DATE)

p8 <- ggplot(ppi_durable, aes(x = DATE, y = PITGCD01JPM661N))+
  ggtitle("Plot of PPI Durable Consumer Goods") +
  xlab("Date") +
  ylab("Producer Prices Index Total Durable Consumer Goods for Japan") +
  theme_economist_white()+
  geom_rect(data = ppi_durable,
            aes(xmin = as.Date("1986-12-01"), xmax = as.Date("1991-03-31"), ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "#99ccff", alpha = 0.2)+ 
  geom_line()+
  theme_economist()
p8

mean(ppi_durable$PITGCD01JPM661N) 
sd(ppi_durable$PITGCD01JPM661N)
skewness(ppi_durable$PITGCD01JPM661N)
kurtosis(ppi_durable$PITGCD01JPM661N)

#10.Producer Prices Index Type of goods Non durable consumer goods Total for Japan
ppi_nondurable <- read.csv("/Users/caiyawei/Desktop/Japan Data/Producer Prices Index Type of goods Non durable consumer goods Total for Japan/JPNPITGND01GYM.csv")
ppi_nondurable <- as.data.frame(ppi_nondurable)
ppi_nondurable$DATE <- as.Date(ppi_nondurable$DATE)

p9 <- ggplot(ppi_nondurable, aes(x = DATE, y = JPNPITGND01GYM))+
  ggtitle("Plot of Producer Prices Index Type of goods Non durable") +
  xlab("Date") +
  ylab("Producer Prices Index Type of goods Non durable consumer goods Total for Japan") +
  theme_economist_white()+
  geom_rect(data = ppi_nondurable,
            aes(xmin = as.Date("1986-12-01"), xmax = as.Date("1991-03-31"), ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "#99ccff", alpha = 0.2)+ 
  geom_line()+
  theme_economist()
p9

mean(ppi_nondurable$JPNPITGND01GYM) 
sd(ppi_nondurable$JPNPITGND01GYM)
skewness(ppi_nondurable$JPNPITGND01GYM)
kurtosis(ppi_nondurable$JPNPITGND01GYM)

#11.Gross Domestic Product Per Capita for Japan
gdp_per_capita <- read.csv("/Users/caiyawei/Desktop/new data/data_for_last_part.csv")
gdp_per_capita <- as.data.frame(gdp_per_capita)
gdp_per_capita$DATE <- as.Date(gdp_per_capita$Year)

p10 <- ggplot(gdp_per_capita, aes(x = DATE, y = GDP))+
  ggtitle("Plot of Gross Domestic Product") +
  xlab("Date") +
  ylab("Gross Domestic Product") +
  theme_economist_white()+
  geom_line()+
  theme_economist()
p10

mean(gdp_per_capita$GDP) 
sd(gdp_per_capita$GDP)
skewness(gdp_per_capita$GDP)
kurtosis(gdp_per_capita$GDP)

#12.Household Debt to GDP
household_debt2gdp<- read.csv("/Users/caiyawei/Desktop/Japan Data/HouseholdDebt to GDP/imf-dm-export-20210324.csv")
household_debt2gdp <- as.data.frame(household_debt2gdp)
household_debt2gdp$Date <- strptime(household_debt2gdp$Date, "%Y")
household_debt2gdp$Date <- as.Date(household_debt2gdp$Date)
class(household_debt2gdp$Date)

p11 <- ggplot(household_debt2gdp, aes(x = Date, y = Japan))+
  ggtitle("Plot of Household Debt to GDP") +
  xlab("Date") +
  ylab("Household Debt to GDP") +
  theme_economist_white()+
  geom_rect(data = household_debt2gdp,
            aes(xmin = as.Date("1986-12-01"), xmax = as.Date("1991-03-31"), ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "#99ccff", alpha = 0.2)+ 
  geom_line()+
  theme_economist()
p11

mean(household_debt2gdp$Japan) 
sd(household_debt2gdp$Japan)
skewness(household_debt2gdp$Japan)
kurtosis(household_debt2gdp$Japan)

#13.Assets and Liabilities of Domestically Licensed Banks (Banking Accounts)(End of Month)
boj_financial_inst <- read.csv("/Users/caiyawei/Desktop/Japan Data/BOJ Data/Financial Istisution/Annually Financial Instisution/nme_R031.32718.20210326220326.01.csv")
boj_financial_inst <- as.data.frame(boj_financial_inst)
boj_financial_inst$Series.code <- strptime(boj_financial_inst$Series.code, "%Y")
boj_financial_inst$Series.code <- as.Date(boj_financial_inst$Series.code)
class(boj_financial_inst$Series.code)
#boj_financial_inst <- drop_na(boj_financial_inst)

p12 <- ggplot(boj_financial_inst, aes(x = Series.code, y = BS02.FAABK_FAAB2DBEA02))+
  ggtitle("Plot of Assets and Liabilities") +
  xlab("Year") +
  ylab("Assets and Liabilities") +
  theme_economist_white()+
  geom_rect(data = boj_financial_inst,
            aes(xmin = as.Date("1986-12-01"), xmax = as.Date("1991-03-31"), ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "#99ccff", alpha = 0.2)+ 
  geom_line()+
  theme_economist()
p12
#boj_financial_inst <- drop_na(boj_financial_inst)
mean(na.omit(boj_financial_inst$BS02.FAABK_FAAB2DBEA02))
sd(na.omit(boj_financial_inst$BS02.FAABK_FAAB2DBEA02))
skewness(na.omit(boj_financial_inst$BS02.FAABK_FAAB2DBEA02))
kurtosis(na.omit(boj_financial_inst$BS02.FAABK_FAAB2DBEA02))


#14.Commercial Paper
boj_financial_inst <- read.csv("/Users/caiyawei/Desktop/Japan Data/BOJ Data/Financial Istisution/Annually Financial Instisution/nme_R031.32718.20210326220326.01.csv")
boj_financial_inst <- as.data.frame(boj_financial_inst)
boj_financial_inst$Series.code <- strptime(boj_financial_inst$Series.code, "%Y")
boj_financial_inst$Series.code <- as.Date(boj_financial_inst$Series.code)
class(boj_financial_inst$Series.code)
#boj_financial_inst <- drop_na(boj_financial_inst)

p13 <- ggplot(boj_financial_inst, aes(x = Series.code, y = BS02.FAABK_FAAB2DBEA13))+
  ggtitle("Plot of Commercial Paper of Assets and Liabilities") +
  xlab("Year") +
  ylab("Commercial Paper)") +
  theme_economist_white()+
  geom_rect(data = boj_financial_inst,
            aes(xmin = as.Date("1986-12-01"), xmax = as.Date("1991-03-31"), ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "#99ccff", alpha = 0.2)+ 
  geom_line()+
  theme_economist()
p13
#boj_financial_inst <- drop_na(boj_financial_inst)
mean(na.omit(boj_financial_inst$BS02.FAABK_FAAB2DBEA13))
sd(na.omit(boj_financial_inst$BS02.FAABK_FAAB2DBEA13))
skewness(na.omit(boj_financial_inst$BS02.FAABK_FAAB2DBEA13))
kurtosis(na.omit(boj_financial_inst$BS02.FAABK_FAAB2DBEA13))

#15.Trading-related Financial Derivatives/Assets

p14 <- ggplot(boj_financial_inst, aes(x = Series.code, y = BS02.FAABK_FAAB2DBEA145))+
  ggtitle("Plot of Trading-related Financial Derivatives/Assets of Assets and Liabilities of Domestically Licensed Banks (Banking Accounts)") +
  xlab("Year") +
  ylab("Trading-related Financial Derivatives/Assets") +
  theme_economist_white()+
  geom_rect(data = boj_financial_inst,
            aes(xmin = as.Date("1986-12-01"), xmax = as.Date("1991-03-31"), ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "#99ccff", alpha = 0.2)+ 
  geom_line()+
  theme_economist()
p14
#boj_financial_inst <- drop_na(boj_financial_inst)
mean(na.omit(boj_financial_inst$BS02.FAABK_FAAB2DBEA145))
sd(na.omit(boj_financial_inst$BS02.FAABK_FAAB2DBEA145))
skewness(na.omit(boj_financial_inst$BS02.FAABK_FAAB2DBEA145))
kurtosis(na.omit(boj_financial_inst$BS02.FAABK_FAAB2DBEA145))

#16.Government Bonds/Assets

p15 <- ggplot(boj_financial_inst, aes(x = Series.code, y = BS02.FAABK_FAAB2DBEA16))+
  ggtitle("Plot of Government Bonds/Assets") +
  xlab("Year") +
  ylab("Government Bonds/Assets") +
  theme_economist_white()+
  geom_rect(data = boj_financial_inst,
            aes(xmin = as.Date("1986-12-01"), xmax = as.Date("1991-03-31"), ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "#99ccff", alpha = 0.2)+ 
  geom_line()+
  theme_economist()
p15
#boj_financial_inst <- drop_na(boj_financial_inst)
mean(na.omit(boj_financial_inst$BS02.FAABK_FAAB2DBEA16))
sd(na.omit(boj_financial_inst$BS02.FAABK_FAAB2DBEA16))
skewness(na.omit(boj_financial_inst$BS02.FAABK_FAAB2DBEA16))
kurtosis(na.omit(boj_financial_inst$BS02.FAABK_FAAB2DBEA16))

#17.Local Government Bonds/Assets

p16 <- ggplot(boj_financial_inst, aes(x = Series.code, y = BS02.FAABK_FAAB2DBEA17))+
  ggtitle("Plot of Local Government Bonds/Assets of Assets and Liabilities of Domestically Licensed Banks (Banking Accounts)") +
  xlab("Year") +
  ylab("Local Government Bonds/Assets") +
  theme_economist_white()+
  geom_rect(data = boj_financial_inst,
            aes(xmin = as.Date("1986-12-01"), xmax = as.Date("1991-03-31"), ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "#99ccff", alpha = 0.2)+ 
  geom_line()+
  theme_economist()
p16
#boj_financial_inst <- drop_na(boj_financial_inst)
mean(na.omit(boj_financial_inst$BS02.FAABK_FAAB2DBEA17))
sd(na.omit(boj_financial_inst$BS02.FAABK_FAAB2DBEA17))
skewness(na.omit(boj_financial_inst$BS02.FAABK_FAAB2DBEA17))
kurtosis(na.omit(boj_financial_inst$BS02.FAABK_FAAB2DBEA17))

#18.Investment Securities/Assets

p17 <- ggplot(boj_financial_inst, aes(x = Series.code, y = BS02.FAABK_FAAB2DBEA21))+
  ggtitle("Plot of Investment Securities/Assets of Assets and Liabilities") +
  xlab("Year") +
  ylab("Investment Securities/Assets") +
  theme_economist_white()+
  geom_rect(data = boj_financial_inst,
            aes(xmin = as.Date("1986-12-01"), xmax = as.Date("1991-03-31"), ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "#99ccff", alpha = 0.2)+ 
  geom_line()+
  theme_economist()
p17
#boj_financial_inst <- drop_na(boj_financial_inst)
mean(na.omit(boj_financial_inst$BS02.FAABK_FAAB2DBEA21))
sd(na.omit(boj_financial_inst$BS02.FAABK_FAAB2DBEA21))
skewness(na.omit(boj_financial_inst$BS02.FAABK_FAAB2DBEA21))
kurtosis(na.omit(boj_financial_inst$BS02.FAABK_FAAB2DBEA21))

#19.Corporate Bonds/Assets

p18 <- ggplot(boj_financial_inst, aes(x = Series.code, y = BS02.FAABK_FAAB2DBEA25))+
  ggtitle("Plot of Corporate Bonds/Assets of Assets and Liabilities") +
  xlab("Year") +
  ylab("Corporate Bonds/Assets") +
  theme_economist_white()+
  geom_rect(data = boj_financial_inst,
            aes(xmin = as.Date("1986-12-01"), xmax = as.Date("1991-03-31"), ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "#99ccff", alpha = 0.2)+ 
  geom_line()+
  theme_economist()
p18
#boj_financial_inst <- drop_na(boj_financial_inst)
mean(na.omit(boj_financial_inst$BS02.FAABK_FAAB2DBEA25))
sd(na.omit(boj_financial_inst$BS02.FAABK_FAAB2DBEA25))
skewness(na.omit(boj_financial_inst$BS02.FAABK_FAAB2DBEA25))
kurtosis(na.omit(boj_financial_inst$BS02.FAABK_FAAB2DBEA25))


#20.Stocks/Assets

p21 <- ggplot(boj_financial_inst, aes(x = Series.code, y = BS02.FAABK_FAAB2DBEA29))+
  ggtitle("Plot of Stocks/Assets of Assets and Liabilities") +
  xlab("Year") +
  ylab("Stocks/Assets") +
  theme_economist_white()+
  geom_rect(data = boj_financial_inst,
            aes(xmin = as.Date("1986-12-01"), xmax = as.Date("1991-03-31"), ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "#99ccff", alpha = 0.2)+ 
  geom_line()+
  theme_economist()
p21
#boj_financial_inst <- drop_na(boj_financial_inst)
mean(na.omit(boj_financial_inst$BS02.FAABK_FAAB2DBEA29))
sd(na.omit(boj_financial_inst$BS02.FAABK_FAAB2DBEA29))
skewness(na.omit(boj_financial_inst$BS02.FAABK_FAAB2DBEA29))
kurtosis(na.omit(boj_financial_inst$BS02.FAABK_FAAB2DBEA29))


#21. Nikkei Industry Research Institute
nikkei <- read.csv("/Users/caiyawei/Desktop/new data/data_for_last_part.csv")
nikkei <- as.data.frame(nikkei)
nikkei$DATE <- as.Date(nikkei$Year)



p20 <- ggplot(nikkei, aes(x = DATE, y = nikkei))+
  ggtitle("Plot of NIKKEI225") +
  xlab("Date") +
  ylab("NIKKEI225") +
  theme_economist_white()+
  geom_line()+
  theme_economist()
p20
#boj_financial_inst <- drop_na(boj_financial_inst)
mean(na.omit(nikkei$nikkei))
sd(na.omit(nikkei$nikkei))
skewness(na.omit(nikkei$nikkei))
kurtosis(na.omit(nikkei$nikkei))


#22.Consumer Price Index: OECD Groups: Housing: Total for Japan
house_price <- read.csv("/Users/caiyawei/Desktop/Japan Data/Housing/JPNCPGRHO01GYM.csv")
house_price <- as.data.frame(house_price)
house_price$DATE <- as.Date(house_price$DATE)

p18 <- ggplot(house_price, aes(x = DATE, y = JPNCPGRHO01GYM))+
  ggtitle("Plot of CPI: Housing") +
  xlab("Year") +
  ylab("CPI: Housing") +
  theme_economist_white()+
  geom_rect(data = house_price,
            aes(xmin = as.Date("1986-12-01"), xmax = as.Date("1991-03-31"), ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "#99ccff", alpha = 0.2)+ 
  geom_line()+
  theme_economist()
p18
#boj_financial_inst <- drop_na(boj_financial_inst)
mean(na.omit(house_price$JPNCPGRHO01GYM))
sd(na.omit(house_price$JPNCPGRHO01GYM))
skewness(na.omit(house_price$JPNCPGRHO01GYM))
kurtosis(na.omit(house_price$JPNCPGRHO01GYM))


#Quartly Data
NIKKEI225 <- read.csv("/Users/caiyawei/Desktop/Quartly Data/NIKKEI/NIKKEI225 (1).csv")
NIKKEI225 <- as.data.frame(NIKKEI225)
NIKKEI225$DATE <- as.Date(NIKKEI225$DATE)

p18 <- ggplot(NIKKEI225, aes(x = DATE, y = NIKKEI225))+
  ggtitle("Plot of NIKKEI225") +
  xlab("Year") +
  ylab("NIKKEI225") +
  theme_economist_white()+
  geom_rect(data = NIKKEI225,
            aes(xmin = as.Date("1986-12-01"), xmax = as.Date("1991-03-31"), ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "#99ccff", alpha = 0.2)+ 
  geom_line()+
  theme_economist()
p18
#boj_financial_inst <- drop_na(boj_financial_inst)
mean(na.omit(NIKKEI225$NIKKEI225))
sd(na.omit(NIKKEI225$NIKKEI225))
skewness(na.omit(NIKKEI225$NIKKEI225))
kurtosis(na.omit(NIKKEI225$NIKKEI225))

#Quartly Data
house_price <- read.csv("/Users/caiyawei/Desktop/Quartly Data/Residential Price/QJPN628BIS (1).csv")
house_price <- as.data.frame(house_price)
house_price$DATE <- as.Date(house_price$DATE)

p18 <- ggplot(house_price, aes(x = DATE, y = QJPN628BIS))+
  ggtitle("Plot of House Price") +
  xlab("Year") +
  ylab("House Price") +
  theme_economist_white()+
  geom_rect(data = house_price,
            aes(xmin = as.Date("1986-12-01"), xmax = as.Date("1991-03-31"), ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "#99ccff", alpha = 0.2)+ 
  geom_line()+
  theme_economist()
p18
#boj_financial_inst <- drop_na(boj_financial_inst)
mean(na.omit(house_price$QJPN628BIS))
sd(na.omit(house_price$QJPN628BIS))
skewness(na.omit(house_price$QJPN628BIS))
kurtosis(na.omit(house_price$QJPN628BIS))

#Quartly Data
interest_rate <- read.csv("/Users/caiyawei/Desktop/Quartly Data/Interest Rate/INTDSRJPM193N (1).csv")
interest_rate <- as.data.frame(interest_rate)
interest_rate$DATE <- as.Date(interest_rate$DATE)
specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))
interest_rate$INTDSRJPM193N <- as.numeric(as.character(interest_rate$INTDSRJPM193N))
interest_rate$INTDSRJPM193N <- specify_decimal(interest_rate$INTDSRJPM193N, 4)
interest_rate$INTDSRJPM193N <- as.numeric(interest_rate$INTDSRJPM193N)
drop_na(interest_rate)
p18 <- ggplot(interest_rate, aes(x = DATE, y = INTDSRJPM193N))+
  ggtitle("Plot of Interest Rate") +
  xlab("Year") +
  ylab("Interest Rate") +
  theme_economist_white()+
  geom_rect(data = interest_rate,
            aes(xmin = as.Date("1986-12-01"), xmax = as.Date("1991-03-31"), ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "#99ccff", alpha = 0.2)+ 
  geom_line()+
  theme_economist()
p18
#boj_financial_inst <- drop_na(boj_financial_inst)
mean(na.omit(interest_rate$INTDSRJPM193N))
sd(na.omit(interest_rate$INTDSRJPM193N))
skewness(na.omit(interest_rate$INTDSRJPM193N))
kurtosis(na.omit(interest_rate$INTDSRJPM193N))

#Quartly Data
employment <- read.csv("/Users/caiyawei/Desktop/Quartly Data/Employment/LREM64TTJPM156S (1).csv")
employment <- as.data.frame(employment)
employment$DATE <- as.Date(employment$DATE)
#specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))
employment$LREM64TTJPM156S <- as.numeric(as.character(employment$LREM64TTJPM156S))
employment$LREM64TTJPM156S <- specify_decimal(employment$LREM64TTJPM156S, 5)
drop_na(employment)
p18 <- ggplot(employment, aes(x = DATE, y = LREM64TTJPM156S))+
  ggtitle("Plot of Employment") +
  xlab("Year") +
  ylab("Employment") +
  theme_economist_white()+
  geom_rect(data = employment,
            aes(xmin = as.Date("1986-12-01"), xmax = as.Date("1991-03-31"), ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "#99ccff", alpha = 0.2)+ 
  geom_line()+
  theme_economist()
p18
#boj_financial_inst <- drop_na(boj_financial_inst)
mean(na.omit(employment$LREM64TTJPM156S))
sd(na.omit(employment$LREM64TTJPM156S))
skewness(na.omit(employment$LREM64TTJPM156S))
kurtosis(na.omit(employment$LREM64TTJPM156S))


#Quartly Data
cpi <- read.csv("/Users/caiyawei/Desktop/Quartly Data/CPI/CPALTT01JPM659N.csv")
cpi <- as.data.frame(cpi)
cpi$DATE <- as.Date(cpi$DATE)
#specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))
cpi$CPALTT01JPM659N <- as.numeric(as.character(cpi$CPALTT01JPM659N))
drop_na(cpi)
#cpi$JPNCPIALLMINMEI <- specify_decimal(cpi$CPALTT01JPM659N, 3)
cpi$CPALTT01JPM659N <- as.numeric(cpi$CPALTT01JPM659N)
drop_na(cpi)
p18 <- ggplot(cpi, aes(x = DATE, y = CPALTT01JPM659N))+
  ggtitle("Plot of CPI") +
  xlab("Year") +
  ylab("CPI") +
  theme_economist_white()+
  geom_rect(data = cpi,
            aes(xmin = as.Date("1986-12-01"), xmax = as.Date("1991-03-31"), ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "#99ccff", alpha = 0.2)+ 
  geom_line()+
  theme_economist()
p18
#boj_financial_inst <- drop_na(boj_financial_inst)
mean(na.omit(cpi$CPALTT01JPM659N))
sd(na.omit(cpi$CPALTT01JPM659N))
skewness(na.omit(cpi$CPALTT01JPM659N))
kurtosis(na.omit(cpi$CPALTT01JPM659N))
######
ppi_nondurable <- read.csv("/Users/caiyawei/Desktop/Quartly Data/producer price/PITGCD01JPM661N (1).csv")
ppi_nondurable <- as.data.frame(ppi_nondurable)
ppi_nondurable$DATE <- as.Date(ppi_nondurable$DATE)
ppi_nondurable$PITGCD01JPM661N <- as.numeric(as.character(ppi_nondurable$PITGCD01JPM661N))
drop_na(ppi_nondurable)

p9 <- ggplot(ppi_nondurable, aes(x = DATE, y = PITGCD01JPM661N))+
  ggtitle("Plot of Producer Prices Index Type of goods Non durable") +
  xlab("Date") +
  ylab("Producer Prices Index Type of goods Non durable consumer goods Total for Japan") +
  theme_economist_white()+
  geom_rect(data = ppi_nondurable,
            aes(xmin = as.Date("1986-12-01"), xmax = as.Date("1991-03-31"), ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "#99ccff", alpha = 0.2)+ 
  geom_line()+
  theme_economist()
p9

mean(na.omit(ppi_nondurable$PITGCD01JPM661N)) 
sd(na.omit(ppi_nondurable$PITGCD01JPM661N))
skewness(na.omit(ppi_nondurable$PITGCD01JPM661N))
kurtosis(na.omit(ppi_nondurable$PITGCD01JPM661N))

#####Seasonal Adjustment#####
###House Price###
library(readxl)
House_price <- read_excel("/Users/caiyawei/Desktop/Quartly Data/Residential Price/QJPN368BIS.xls")
House_price$D1 <- ifelse(House_price$Season == '1', 1, 0)
House_price$D2 <- ifelse(House_price$Season == '2', 1, 0)
House_price$D3 <- ifelse(House_price$Season == '3', 1, 0)
House_price$D4 <- ifelse(House_price$Season == '4', 1, 0)
House_price$DATE <- as.Date(House_price$observation_date)
model <- lm(formula = House_price$QJPN368BIS_PCH ~ House_price$D1 + House_price$D2 + House_price$D3 + House_price$D4 -1)
House_price$res = resid(model)
p10 <- ggplot(House_price, aes(x = DATE, y = res))+
  ggtitle("Plot of Real Residential Property Price") +
  xlab("Date") +
  ylab("Real Residential Property Price") +
  theme_economist_white()+
  geom_rect(data = House_price,
            aes(xmin = as.Date("1986-12-01"), xmax = as.Date("1991-03-31"), ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "#99ccff", alpha = 0.2)+ 
  geom_line()+
  theme_economist()
p10
write.csv(House_price,file="/Users/caiyawei/Desktop/Quartly Data/Seasonal Adjustment/Residential Price/House_price.csv",row.names = FALSE)

mean(na.omit(House_price$res))
sd(na.omit(House_price$res))
skewness(na.omit(House_price$res))
kurtosis(na.omit(House_price$res))
###GDP###

gdp <- read.csv('/Users/caiyawei/Desktop/Quartly Data/GDP Original Series for Japan/LORSGPORJPQ659S.csv')
gdp$D1 <- ifelse(gdp$Season == '1', 1, 0)
gdp$D2 <- ifelse(gdp$Season == '2', 1, 0)
gdp$D3 <- ifelse(gdp$Season == '3', 1, 0)
gdp$D4 <- ifelse(gdp$Season == '4', 1, 0)
#gdp$observation_date<-as.character(gdp$observation_date)
gdp$DATE <- as.Date(gdp$observation_date)
model <- lm(formula = gdp$LORSGPORJPQ659S ~ gdp$D1 + gdp$D2 + gdp$D3 + gdp$D4 -1)
gdp$res = resid(model)
p11 <- ggplot(gdp, aes(x = DATE, y = res))+
  ggtitle("Plot of GDP") +
  xlab("Date") +
  ylab("GDP") +
  theme_economist_white()+
  geom_rect(data = gdp,
            aes(xmin = as.Date("1986-12-01"), xmax = as.Date("1991-03-31"), ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "#99ccff", alpha = 0.2)+ 
  geom_line()+
  theme_economist()
p11
write.csv(gdp,file="/Users/caiyawei/Desktop/Quartly Data/Seasonal Adjustment/GDP/GDP.csv",row.names = FALSE)

mean(na.omit(gdp$res))
sd(na.omit(gdp$res))
skewness(na.omit(gdp$res))
kurtosis(na.omit(gdp$res))

###Interest###
interest <- read.csv('/Users/caiyawei/Desktop/Quartly Data/Interest Rate/INTDSRJPM193N .csv')
interest$D1 <- ifelse(interest$Season == '1', 1, 0)
interest$D2 <- ifelse(interest$Season == '2', 1, 0)
interest$D3 <- ifelse(interest$Season == '3', 1, 0)
interest$D4 <- ifelse(interest$Season == '4', 1, 0)
#gdp$observation_date<-as.character(gdp$observation_date)
interest$DATE <- as.Date(interest$DATE)
model <- lm(formula = interest$INTDSRJPM193N ~ interest$D1 + interest$D2 + interest$D3 + interest$D4 -1)
interest$res <- resid(model)
p11 <- ggplot(interest, aes(x = DATE, y = res))+
  ggtitle("Plot of Interest Rate") +
  xlab("Date") +
  ylab("Interest Rate") +
  theme_economist_white()+
  geom_rect(data = interest,
            aes(xmin = as.Date("1986-12-01"), xmax = as.Date("1991-03-31"), ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "#99ccff", alpha = 0.2)+ 
  geom_line()+
  theme_economist()
p11
write.csv(interest,file="/Users/caiyawei/Desktop/Quartly Data/Seasonal Adjustment/interest/interest.csv",row.names = FALSE)

mean(na.omit(interest$res))
sd(na.omit(interest$res))
skewness(na.omit(interest$res))
kurtosis(na.omit(interest$res))

###CPI###
cpi <- read.csv('/Users/caiyawei/Desktop/Quartly Data/CPI/CPALTT01JPM659N.csv')
cpi$D1 <- ifelse(cpi$Season == '1', 1, 0)
cpi$D2 <- ifelse(cpi$Season == '2', 1, 0)
cpi$D3 <- ifelse(cpi$Season == '3', 1, 0)
cpi$D4 <- ifelse(cpi$Season == '4', 1, 0)
#gdp$observation_date<-as.character(gdp$observation_date)
cpi$DATE <- as.Date(cpi$DATE)
model <- lm(formula = cpi$CPALTT01JPM659N ~ cpi$D1 + cpi$D2 + cpi$D3 +  cpi$D4 -1)
cpi$res <- resid(model)
p11 <- ggplot(cpi, aes(x = DATE, y = res))+
  ggtitle("Plot of CPI") +
  xlab("Date") +
  ylab("CPI") +
  theme_economist_white()+
  geom_rect(data = cpi,
            aes(xmin = as.Date("1986-12-01"), xmax = as.Date("1991-03-31"), ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "#99ccff", alpha = 0.2)+ 
  geom_line()+
  theme_economist()
p11
write.csv(cpi,file="/Users/caiyawei/Desktop/Quartly Data/Seasonal Adjustment/cpi/cpi.csv",row.names = FALSE)
mean(na.omit(cpi$res))
sd(na.omit(cpi$res))
skewness(na.omit(cpi$res))
kurtosis(na.omit(cpi$res))
###NIKKEI###
NIKKEI <- read.csv('/Users/caiyawei/Desktop/Quartly Data/NIKKEI/NIKKEI225 (1).csv')
NIKKEI$D1 <- ifelse(NIKKEI$Season == '1', 1, 0)
NIKKEI$D2 <- ifelse(NIKKEI$Season == '2', 1, 0)
NIKKEI$D3 <- ifelse(NIKKEI$Season == '3', 1, 0)
NIKKEI$D4 <- ifelse(NIKKEI$Season == '4', 1, 0)
#gdp$observation_date<-as.character(gdp$observation_date)
NIKKEI$DATE <- as.Date(NIKKEI$DATE)
model <- lm(formula = NIKKEI$NIKKEI225 ~ NIKKEI$D1 + NIKKEI$D2 + NIKKEI$D3 +  NIKKEI$D4 -1)
NIKKEI$res <- resid(model)
p11 <- ggplot(NIKKEI, aes(x = DATE, y = res))+
  ggtitle("Plot of NIKKEI") +
  xlab("Date") +
  ylab("NIKKEI") +
  theme_economist_white()+
  geom_rect(data = NIKKEI,
            aes(xmin = as.Date("1986-12-01"), xmax = as.Date("1991-03-31"), ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "#99ccff", alpha = 0.2)+ 
  geom_line()+
  theme_economist()
p11
write.csv(NIKKEI,file="/Users/caiyawei/Desktop/Quartly Data/Seasonal Adjustment/NIKKEI/NIKKEI.csv",row.names = FALSE)

mean(na.omit(NIKKEI$res))
sd(na.omit(NIKKEI$res))
skewness(na.omit(NIKKEI$res))
kurtosis(na.omit(NIKKEI$res))
###Employment###
employment <- read.csv('/Users/caiyawei/Desktop/Quartly Data/Employment/LREM64TTJPM156S (1).csv')
employment$D1 <- ifelse(employment$Season == '1', 1, 0)
employment$D2 <- ifelse(employment$Season == '2', 1, 0)
employment$D3 <- ifelse(employment$Season == '3', 1, 0)
employment$D4 <- ifelse(employment$Season == '4', 1, 0)
#gdp$observation_date<-as.character(gdp$observation_date)
employment$DATE <- as.Date(employment$DATE)
model <- lm(formula = employment$LREM64TTJPM156S ~ employment$D1 + employment$D2 + employment$D3 +  employment$D4 -1)
employment$res <- resid(model)
p11 <- ggplot(employment, aes(x = DATE, y = res))+
  ggtitle("Plot of Employment") +
  xlab("Date") +
  ylab("Employment") +
  theme_economist_white()+
  geom_rect(data = employment,
            aes(xmin = as.Date("1986-12-01"), xmax = as.Date("1991-03-31"), ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "#99ccff", alpha = 0.2)+ 
  geom_line()+
  theme_economist()
p11
write.csv(employment,file="/Users/caiyawei/Desktop/Quartly Data/Seasonal Adjustment/employment/employment.csv",row.names = FALSE)
mean(na.omit(employment$res))
sd(na.omit(employment$res))
skewness(na.omit(employment$res))
kurtosis(na.omit(employment$res))


















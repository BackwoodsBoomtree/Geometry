library(mcr)


x <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
y <- c(0, 2, 1, 3, 4, 6, 5, 7, 8, 10)

dem.reg <- mcreg(x, y, method.reg = "Deming")

str(dem.reg)

dem.reg@para

plot(x,y, main = "Regression Comparison", xlab = "Current Method", ylab = "New Method")
abline(dem.reg@para[1:2], col = "blue")

MCResult.plot(dem.reg, equal.axis = TRUE)

i <- dem.reg@para[1]
s <- dem.reg@para[2]
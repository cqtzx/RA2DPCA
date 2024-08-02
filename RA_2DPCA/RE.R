A <- readMat("Yale.mat")
x <- A$x
label <- A$label
pcs <- 25
ss <- seq(5,25,5)
ratio <- 0.4
size_x <- dim(x)
height <- size_x[1]; width <- size_x[2]; n <- size_x[3]
n1 <- n
ix1 <- randperm(n1)
ix_noise <- ix1[1:floor(n1*ratio)]
ix_clean <- seq(1:n1)[-ix_noise]
label_noise <- label[ix_noise]
label_clean <- label[ix_clean]
x_train <- x
x_n <- x_train[,,ix_noise]
x_clean <- x_train[,,ix_clean]
for(i in 1:length(ix_noise)){
  k1 <- sample(1:(height-floor(ratio*height)),1)
  k11 <- k1+floor(ratio*height)-1
  k2 <- sample(1:(width-floor(ratio*width)),1)
  k22 <- k2+floor(ratio*width)-1
  noise <- matrix(sample(c(0,255),floor(ratio*height)*floor(ratio*width),replace = T),floor(ratio*height),floor(ratio*width))
  x_n[k1:k11,k2:k22,i] <- noise
}
x_train[,,ix_noise] <- x_n
#2DPCA
pca2d <- RE_2DPCA(x_train,ix_clean,pcs,ss)
pra <- quantile(x_train, probs = 0.9)
re1 <- RE_RA_2DPCA(x_train,ix_clean,pra,1,pcs,ss)
re2 <- RE_OMRA_2DPCA(x_train,ix_clean,pra,1,pcs,ss)

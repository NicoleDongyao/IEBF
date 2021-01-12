#读取原图片
library(jpeg)
img = readJPEG("..input/test_nlm/lena.jpg")

#图片加噪声
A = rnorm(n = 426 * 498, mean=0, sd=1)
A = matrix(A, nrow=426, ncol=498)
noise = A * (50^2 / 255^2)
B = img + noise
m = nrow(B)
n = ncol(B)
for (i in 1:m) {
  for (j in 1:n){
    if (B[i, j] < 0){
      B[i, j] = 0
    }
    if (B[i, j] > 1){
      B[i, j] = 1
    }
  }
}
img_noise = B

#图片处理
image = nlm(img,21,5,0.001)

#对比结果（原图片，噪声图片，处理后的图片）
par(mar = c(0,0,0,0))
plot(0,0,type = "n")
rasterImage(img_noise, -1, -1, 0, 1)
rasterImage(image, 0, -1, 1, 1)
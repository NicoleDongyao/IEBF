library(jpeg)
img0 <- readJPEG("fptest.jpg")
img0 <- img0[, , 1]

# library(rtiff)
# img0 <- readTiff("101_2.tif")
# img0 = img0@red

# dispose
img <- normal(img0)
pic_con(img0, img)

imgg <- gabor(img)
pic_con(img0, imgg)

imgl <- laplace(img)
pic_con(img0, imgl)

imgv = volterra(img)
pic_con(img0, imgv)

## evaluate
#distribution
library(ggplot2)
img_ <- as.vector(img0)
df <- data.frame(pixel = img_)
p <- ggplot(df, aes(x = pixel))
p + geom_density(color = "black", fill = "gray")
mean(img_)

# ordinary
img0t <- binthin(img0, 0.87)

# gabor
imgt <- binthin(imgg, 0.87)
simi_score(img0t, imgt)

# laplacian
imgt <- binthin(imgl, 0.87)
simi_score(img0t, imgt)

# volterra
imgt <- binthin(imgv, 0.87)
simi_score(img0t, imgt)

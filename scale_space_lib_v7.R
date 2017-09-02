#======================================================================================
# Functions: scale space analysis and more
#======================================================================================

# Number of points on the circle to draw
Nd <- 200
phid <- seq(0,2*pi,length.out=Nd)

histstretch<-function(data) {
  cur.lim<-quantile(data,c(0.025,0.975),na.rm=TRUE)
  data<-pmax(cur.lim[1],pmin(cur.lim[2],data))
  data<-floor(255*(data-cur.lim[1])/(cur.lim[2]-cur.lim[1]))
  data[is.na(data)]<-0
  data
}

# Function to convert geographic coordinates into image row,column
# input: range in geographic coordinates
# requires xr,yr, ps, ysign to be set prior to calling this function
# Note: the code has to be properly tested
xy_to_rowcol_old <- function(xyr){
	xrl <- xyr[,1]	# x range
	yrl <- xyr[,2]	# y range

	i1 <- floor((xrl[1] - xr[1])/ps[1])
	i2 <- ceiling((xrl[2] - xr[1])/ps[1])

	if(ysign<0)
	{
		j1 <- floor((yr[2]-yrl[2])/ps[2])
		j2 <- ceiling((yr[2]-yrl[1])/ps[2])
	}else
	{
		j1 <- floor((yrl[1] - yr[1])/ps[2])
		j2 <- ceiling((yrl[2] - yr[1])/ps[2])
	}
	
	return(c(i1,i2,j1,j2))
}

xy_to_rowcol <- function(xyr,imageinfo){

	N0 <- imageinfo[["rows"]]
	M0 <- imageinfo[["columns"]]
	Nb <- imageinfo[["bands"]]

	xyLL <- c(imageinfo[["ll.x"]],imageinfo[["ll.y"]])
	ps <- c(imageinfo[["res.x"]],imageinfo[["res.x"]])

	ysign <- attributes(imageinfo)$ysign
	proj_raster <- attributes(imageinfo)$projection

	max_val_ms <- max(attributes(imageinfo)$df$Bmax)

	# Entire coordinate range of the image
	xr <- c(xyLL[1],xyLL[1]+M0*ps[1])

	if(ysign<0)
	{
		yr <- c(xyLL[2],xyLL[2]+N0*ps[2])
	}else
	{
		yr <- c(xyLL[2]-N0*ps[2],xyLL[2])
	}
	
	xrl <- xyr[,1]	# x range
	yrl <- xyr[,2]	# y range

	i1 <- floor((xrl[1] - xr[1])/ps[1])
	i2 <- ceiling((xrl[2] - xr[1])/ps[1])

	if(ysign<0)
	{
		j1 <- floor((yr[2]-yrl[2])/ps[2])
		j2 <- ceiling((yr[2]-yrl[1])/ps[2])
	}else
	{
		j1 <- floor((yrl[1] - yr[1])/ps[2])
		j2 <- ceiling((yrl[2] - yr[1])/ps[2])
	}

	#i1 <- max(i1,1)
	#j1 <- max(j1,1)
	
	#i2 <- min(M0,i2)
	#j2 <- min(N0,j2)
	
	return(c(i1,i2,j1,j2))
}


# Function to read image subset
read_subset <- function(imagefn,i1,i2,j1,j2){

	work.region.tif <- new("GDALReadOnlyDataset", paste(Path_in,"/",imagefn,sep=""))

	val <- work.region.tif[j1:j2,i1:i2]
	# --- close image
	GDAL.close(work.region.tif)
	
	return(val)
}


spat_average_v3 <- function(P,W)
{
	d <- dim(W)
	nx <- d[1]
	ny <- d[2]

	d <- dim(P)
	M <- d[1]
	N <- d[2]

	kx <- (nx-1)/2
	ky <- (ny-1)/2

	
	Pcum <- P
	Pcum[,] <- 0
	
	for(dx in -kx:kx)
	for(dy in -ky:ky){
	
		Pshift <- P
		
		if(dx!=0)if(dx>0){
			Pshift[1:(M-dx),] <- Pshift[(1+dx):M,]
			Pshift[(M-dx+1):M,] <- 0    # zero padding
		}else{
			Pshift[(1-dx):M,] <- Pshift[1:(M+dx),]
			Pshift[(M+dx+1):M,] <- 0    # zero padding
		}
		
		if(dy!=0)if(dy>0){
			Pshift[,1:(N-dy)] <- Pshift[,(1+dy):N]
			Pshift[,(N-dy+1):N] <- 0    # zero padding
		}else{
			Pshift[,(1-dy):N] <- Pshift[,1:(N+dy)]
			Pshift[,(N+dy+1):N] <- 0    # zero padding
		}

		Pcum <- Pcum + Pshift*W[dx+kx+1,dy+ky+1]
	}

	#Pcum[1:kx,] <- 0
	#Pcum[(M-kx+1):M,] <- 0

	#Pcum[,1:ky] <- 0
	#Pcum[,(N-ky+1):N] <- 0
	
	return(Pcum)
}


# Computes determinant of Hessian matrix (or Laplacian) for a fixed value of scale t
det_Hessian <- function(P,t){
	
	val <- array(0,dim(P))
	
	sigma <- sqrt(c(t,t))

	#kern_size <- pmax(6*sigma,c(3,3))
	kern_size <- 6*sigma
	if(min(kern_size)<3)return(val)

	W <- gaussianKernel(sigma, dim = length(sigma), size = kern_size, normalised = TRUE)
	
	nx <- dim(W)[1]
	ny <- dim(W)[2]

	i0 <- (nx-1)/2
	j0 <- (ny-1)/2

	dy <- array(rep(1:nx,times=ny),c(nx,ny))
	dx <- array(rep(1:ny,each=nx),c(nx,ny))

	dx <- dx-(i0+1)
	dy <- dy-(j0+1)

	Wxx <- W * ((dx^2)/t - 1) / t
	Wxy <- W * dx*dy/(t^2)
	Wyy <- W * ((dy^2)/t - 1) / t

	#Gxx <- morph(P,Wxx,operator="*")
	#Gxy <- morph(P,Wxy,operator="*")
	#Gyy <- morph(P,Wyy,operator="*")

	Gxx <- spat_average_v3(P,Wxx)
	Gxy <- spat_average_v3(P,Wxy)
	Gyy <- spat_average_v3(P,Wyy)

	val <- (t^2)*(Gxx*Gyy - (Gxy)^2)

	return(val)
}

# Computes determinant of Hessian matrix (or Laplacian) for a fixed value of scale t
det_Hessian_v1 <- function(P,t){
	
	val <- array(0,dim(P))
	
	sigma <- sqrt(c(t,t))

	#kern_size <- pmax(6*sigma,c(3,3))
	kern_size <- 6*sigma
	if(min(kern_size)<3)return(val)

	W <- gaussianKernel(sigma, dim = length(sigma), size = kern_size, normalised = TRUE)
	
	nx <- dim(W)[1]
	ny <- dim(W)[2]

	i0 <- (nx-1)/2
	j0 <- (ny-1)/2

	dy <- array(rep(1:nx,times=ny),c(nx,ny))
	dx <- array(rep(1:ny,each=nx),c(nx,ny))

	dx <- dx-(i0+1)
	dy <- dy-(j0+1)

	Wxx <- W * ((dx^2)/t - 1) / t
	Wxy <- W * dx*dy/(t^2)
	Wyy <- W * ((dy^2)/t - 1) / t

	Gxx <- convolve_2d(P,Wxx)
	Gxy <- convolve_2d(P,Wxy)
	Gyy <- convolve_2d(P,Wyy)
	
	val <- (t^2)*(Gxx*Gyy - (Gxy)^2)

	d <- dim(P)
	M<-d[1]
	N<-d[2]
	
	val <- val[(i0+1):(i0+M),(j0+1):(j0+N)]
	
	return(val)
}

# Computes determinant of Hessian matrix (or Laplacian) for a fixed value of scale t
det_Hessian_v2 <- function(P,t){
	
	val <- array(0,dim(P))
	
	sigma <- sqrt(c(t,t))

	#kern_size <- pmax(6*sigma,c(3,3))
	kern_size <- 6*sigma
	if(min(kern_size)<3)return(val)

	W <- gaussianKernel(sigma, dim = length(sigma), size = kern_size, normalised = TRUE)
	
	nx <- dim(W)[1]
	ny <- dim(W)[2]

	i0 <- (nx-1)/2
	j0 <- (ny-1)/2

	dy <- array(rep(1:nx,times=ny),c(nx,ny))
	dx <- array(rep(1:ny,each=nx),c(nx,ny))

	dx <- dx-(i0+1)
	dy <- dy-(j0+1)

	Wxx <- W * ((dx^2)/t - 1) / t
	Wxy <- W * dx*dy/(t^2)
	Wyy <- W * ((dy^2)/t - 1) / t

	Gxx <- convolve_2d_Eigen(P,Wxx)
	Gxy <- convolve_2d_Eigen(P,Wxy)
	Gyy <- convolve_2d_Eigen(P,Wyy)
	
	val <- (t^2)*(Gxx*Gyy - (Gxy)^2)
	
	d <- dim(P)
	M<-d[1]
	N<-d[2]
	
	val <- val[(i0+1):(i0+M),(j0+1):(j0+N)]

	return(val)
}

# Computes determinant of Hessian matrix (or Laplacian) for a fixed value of scale t
det_Hessian_v3 <- function(P,t){
	
	val <- array(0,dim(P))
	
	sigma <- sqrt(c(t,t))

	#kern_size <- pmax(6*sigma,c(3,3))
	kern_size <- 6*sigma
	if(min(kern_size)<3)return(val)

	W <- gaussianKernel(sigma, dim = length(sigma), size = kern_size, normalised = TRUE)
	
	nx <- dim(W)[1]
	ny <- dim(W)[2]

	i0 <- (nx-1)/2
	j0 <- (ny-1)/2

	dy <- array(rep(1:nx,times=ny),c(nx,ny))
	dx <- array(rep(1:ny,each=nx),c(nx,ny))

	dx <- dx-(i0+1)
	dy <- dy-(j0+1)

	Wxx <- W * ((dx^2)/t - 1) / t
	Wxy <- W * dx*dy/(t^2)
	Wyy <- W * ((dy^2)/t - 1) / t


	Gxx <- convolve_2d_tricks(P,Wxx)
	Gxy <- convolve_2d_tricks(P,Wxy)
	Gyy <- convolve_2d_tricks(P,Wyy)
	
	val <- (t^2)*(Gxx*Gyy - (Gxy)^2)
	
	d <- dim(P)
	M<-d[1]
	N<-d[2]
	
	val <- val[(i0+1):(i0+M),(j0+1):(j0+N)]

	return(val)
}


# unsigned Hessian feature strength measure I; for a fixed value of scale t
Hessian_feature_I <- function(P,k,t){
	
	sigma <- sqrt(c(t,t))

	W <- gaussianKernel(sigma, dim = length(sigma), size = 6*sigma, normalised = TRUE)

	nx <- dim(W)[1]
	ny <- dim(W)[2]

	i0 <- (nx-1)/2
	j0 <- (ny-1)/2

	dy <- array(rep(1:nx,times=ny),c(nx,ny))
	dx <- array(rep(1:ny,each=nx),c(nx,ny))

	dx <- dx-(i0+1)
	dy <- dy-(j0+1)

	Wxx <- W * ((dx^2)/t - 1) / t
	Wxy <- W * dx*dy/(t^2)
	Wyy <- W * ((dy^2)/t - 1) / t
			
	Gxx <- morph(P,Wxx,operator="*")
	Gxy <- morph(P,Wxy,operator="*")
	Gyy <- morph(P,Wyy,operator="*")
			
	DH <- (t^2)*(Gxx*Gyy - (Gxy)^2)
	#DH <- (t)*(Gxx+Gyy)

	val <- DH - (t^2)*k*(Gxx+Gyy)^2
	val[val<0] <- 0
	
	return(val)
}


# Computes DOH in scale space tarr
DOH_space <- function(P,tarr){

	Ns <- length(tarr)

	detHes <- array(0,c(dim(P),Ns))

	for(k in 1:Ns){
		detHes[,,k] <- det_Hessian_v2(P,tarr[k])
		#cat(k,"out of ",Ns,"\n")
	}
	return(detHes)
}

# Computes scale space tarr; differential operator is unsigned Hessian feature strength measure I
Scale_space <- function(P,tarr){

	k_par <- 0.04

	Ns <- length(tarr)

	val <- array(0,c(dim(P),Ns))

	for(k in 1:Ns){
		val[,,k] <- Hessian_feature_I(P,k_par,tarr[k])
		
		#cat(k,"out of ",Ns,"\n")
	}
	return(val)
}


# Find local maxima in 1D
# returns only the first of multiple adjacent maxima
localMaxima <- function(x)
{
	y <- diff(c(-Inf, x)) > 0L
	y <- cumsum(rle(y)$lengths)
	y <- y[seq.int(1L, length(y), 2L)]

	if (x[[1]] == x[[2]]) {
		y <- y[-1]
	}
	
	y
}

# get all adjacent maxima
expandlocalmax <- function(x,y){
	z <- array(0,0)
	for(i in 1:length(y))
	{
		k <- y[i]
		z <- c(z,k)
		
		j<-k+1
		if(j<=length(x)) while(x[j]==x[k])
		{
			z <- c(z,j)
			j<- j+1
			if(j>=length(x)) break
		}
	}
	z
}


# Detect local maxima in 2D array
local_max_2D <- function(F){
	
	d <- dim(F)
	
	M <- d[1]
	N <- d[2]
	
	cand <- array(0,d)
	cand1 <- cand
	cand2 <- cand

	for(i in 2:(M-1))
	{
		val <- F[i,]
		
		y <- localMaxima(val)

		if(length(y)>0)
		{
			z <- expandlocalmax(val,y)
			
			cand1[i,z] <- 1
		}
		
	}

	for(j in 2:(N-1))
	{
		val <- F[,j]
		
		y <- localMaxima(val)

		if(length(y)>0)
		{
			z <- expandlocalmax(val,y)
			
			cand2[z,j] <- 1
		}

	}
	
	cand[,] <- cand1 & cand2

	ij <- which(cand==1,arr.ind=TRUE)

	return(ij)
}

# Detect local maxima in 2D array
# Vectorized version of v1
local_max_2D_v2 <- function(F,M,N){
	
	n <- length(F)
	
	cand <- array(0,n)
	cand1 <- cand
	cand2 <- cand

	for(i in 2:(M-1))
	{
		val <- F[i,]
		
		y <- localMaxima(val)

		if(length(y)>0)
		{
			z <- expandlocalmax(val,y)
			
			cand1[i,z] <- 1
		}
		
	}

	for(j in 2:(N-1))
	{
		val <- F[,j]
		
		y <- localMaxima(val)

		if(length(y)>0)
		{
			z <- expandlocalmax(val,y)
			
			cand2[z,j] <- 1
		}

	}
	
	cand[,] <- cand1 & cand2

	ij <- which(cand==1,arr.ind=TRUE)

	return(ij)
}


# Detect local maxima in 3D array
local_max_3D <- function(F){
	pts <- array(0,c(0,5))
	
	for(k in (Ns-1):2)
	{
		ij <- local_max_2D(F[,,k])
		
		if(nrow(ij)>0)
		for(m in 1:nrow(ij))
		{
			i <- ij[m,1]
			j <- ij[m,2]
			j <- ij[m,2]
			
			tmp <- F[i,j,(k-1):(k+1)]
			v0 <- F[i,j,k]

			if(v0>0)
			{
				tmp <- tmp - v0
				if((tmp[1]<=0)&(tmp[3]<=0))
				{
					# Point is accepted
					pts <- rbind(pts,c(i,j,k,v0,1))
				}
			}
		}
	}	
	return(pts)
}



# Detect local maxima in 3D array
# Vectorized version of v1
local_max_3D_v2 <- function(F){
	pts <- array(0,c(0,5))
	
	for(k in (Ns-1):2)
	{
		ij <- local_max_2D(F[,k])
		
		if(nrow(ij)>0)
		for(m in 1:nrow(ij))
		{
			i <- ij[m,1]
			j <- ij[m,2]
			
			tmp <- F[,(k-1):(k+1)]
			v0 <- F[i,j,k]

			if(v0>0)
			{
				tmp <- tmp - v0
				if((tmp[1]<=0)&(tmp[3]<=0))
				{
					# Point is accepted
					pts <- rbind(pts,c(i,j,k,v0,1))
				}
			}
		}
	}	
	return(pts)
}


# 2nd order polynomial interpolation of the maximum based on 3 points
parabol_interp <- function(x,y)
{
	a <- array(0,3)
	if((length(x)==3)&(length(y)==3))
	{
		a[3] <- (2*(y[2]-y[1])*(x[3]+x[1])-(y[3]-y[1])*(x[2]+x[1]))/(2*(2*y[2]-y[1]-y[3]))
		a[2] <- (y[2]-y[1])/((x[2]-x[1])*(x[1]+x[2]-2*a[3]))
		a[1] <- y[1] - a[2]*(x[1]-a[3])^2
		val <- a
	}else(val <- -1)
	return(val)
}

# 2nd order polynomial interpolation of the maximum based on 3 points
# Vectorized version: x and y are 2D matrices with 3 columns; rows are the extrema points
parabol_interp_vect <- function(x,y)
{
	if(!all.equal(dim(y),dim(x))) return(NA)

	n <- nrow(x)
	
	if(n==1)
	a <- array(0,c(n,3))
	a[,3] <- (2*(y[,2]-y[,1])*(x[,3]+x[,1])-(y[,3]-y[,1])*(x[,2]+x[,1]))/(2*(2*y[,2]-y[,1]-y[,3]))
	a[,2] <- (y[,2]-y[,1])/((x[,2]-x[,1])*(x[,1]+x[,2]-2*a[,3]))
	a[,1] <- y[,1] - a[,2]*(x[,1]-a[,3])^2

	return(a)
}


# 2nd order polynomial interpolation of the maximum based on n points
parabol_interp_v2 <- function(x,y)
{
	a <- array(0,3)
	if((length(x)<3)&(length(y)<3)) return(NA)

	xav   <- mean(x) 
	x2av  <- mean(x^2)
	x3av  <- mean(x^3)
	x4av  <- mean(x^4)
	yav   <- mean(y)
	xyav  <- mean(x*y)
	x2yav <- mean((x^2)*y)
	
	D <- rbind(c(x2av,xav,1),c(x3av,x2av,xav),c(x4av,x3av,x2av))
	vec <- c(yav,xyav,x2yav)
	
	Da <- D
	Da[,1] <- vec
	Db <- D
	Db[,2] <- vec
	Dc <- D
	Dc[,3] <- vec
	
	detD <- det(D)
	a <- det(Da)/detD
	b <- det(Db)/detD
	c <- det(Dc)/detD

	if(Debug){
		plot(x,y)
		xdisp <- seq(min(x),max(x),length.out=250)
		lines(xdisp,a*xdisp^2+b*xdisp+c,col="blue")
	}

	# Is that a maximum?
	if(abs(a)<1e-08){
		val <- x[which.max(y)]
		return(val)
	}

	if(a>0){
		val <- x[which.max(y)]
		return(val)
	}

	val <- -b/(2*a)

	return(val)
}

determine_max <- function(x,y){
	
	# Interpolate the max
	ind <- which(y==max(y))
	if(length(ind)>1){
		# Problem, requires handling. Especially if the function has multiple peaks!
		if(Debug)plot(x,y)
		return(mean(x[ind]))
	}
	
	j1 <- max(c(ind-3,1))
	j2 <- min(c(ind+3,length(x)))

	x1 <- x[j1:j2]
	y1 <- y[j1:j2]
	#plot(x1,y1)

	val <- interpolate_max(x1,y1)
	
	return(val)
}

determine_max_v2 <- function(x,y){

	#if(Debug)plot(x,y)
	
	# Interpolate the max
	ind <- which(y==max(y))
	if(length(ind)>1){
		# Problem, requires handling. Especially if the function has multiple peaks!
		return(mean(x[ind]))
	}
	
	if(ind==1 | ind==length(y)) return(x[ind])
	
	j1 <- max(c(ind-1,1))
	j2 <- min(c(ind+1,length(x)))

	x1 <- x[j1:j2]
	y1 <- y[j1:j2]

	a <- parabol_interp(x1,y1)
	
	val <- a[3]
	
	return(val)
}


# Display detected blobs over the image 
# Note: requires the image to be open in active window!
display_blobs <- function(ind,color){
	for(k in ind)
	{
		x0 <- Special_points[k,1]
		y0 <- Special_points[k,2]
		rd <- rep(0.5*Special_points[k,3],Nd)
		
		xd <- x0 + rd*cos(phid)
		yd <- y0 + rd*sin(phid)

		points(x0,y0,col=color,pch=16)
		lines(xd,yd,col=color,lwd=1)
	}
}

# Display best blobs over the image 
# Note: requires the image to be open in active window!
display_best_blobs <- function(ind,color){
	for(k in ind)
	{
		x0 <- Best_points[k,1]
		y0 <- Best_points[k,2]
		rd <- rep(0.5*Best_points[k,3],Nd)
		
		xd <- x0 + rd*cos(phid)
		yd <- y0 + rd*sin(phid)

		points(x0,y0,col=color,pch=16)
		lines(xd,yd,col=color,lwd=1)
	}
}

# Display all blobs over the entire image 
# Note: requires the image to be open in active window!

display_sel_blobs <- function(Blobs,ind,color){
	for(k in ind)
	{
		x0 <- Blobs[k,1]
		y0 <- Blobs[k,2]
		rd <- rep(0.5*Blobs[k,3],Nd)
		
		xd <- x0 + rd*cos(phid)
		yd <- y0 + rd*sin(phid)

		points(x0,y0,col=color,pch=16)
		lines(xd,yd,col=color,lwd=1)
	}
}

display_all_blobs <- function(Blobs,color){
	for(k in 1:nrow(Blobs))
	{
		x0 <- Blobs$x[k]
		y0 <- Blobs$y[k]
		rd <- rep(0.5*Blobs$d[k],Nd)
		
		xd <- x0 + rd*cos(phid)
		yd <- y0 + rd*sin(phid)

		points(x0,y0,col=color,pch=16)
		lines(xd,yd,col=color,lwd=1)
	}
}

display_blobs_v2 <- function(Blobs,color,lwd){
	for(k in 1:nrow(Blobs))
	{
		x0 <- Blobs$x[k]
		y0 <- Blobs$y[k]
		rd <- rep(0.5*Blobs$d[k],Nd)
		
		xd <- x0 + rd*cos(phid)
		yd <- y0 + rd*sin(phid)

		points(x0,y0,col=color,pch=16)
		lines(xd,yd,col=color,lwd=lwd)
	}
}


display_image <- function(A,nR,Ng,nB,title_str){
	
	d <- A@grid@cells.dim

	M <- d[1]
	N <- d[2]
	
	ps <- A@grid@cellsize

	fx <- M/800
	fy <- N/600
	
	#f <- 1/ceiling(max(fx,fy))
	f <- 1/4

	tmp <- data.matrix(A@data[,c(nR,nG,nB)], rownames.force = NA)
	yR <- array(tmp[,1],c(M,N))
	yG <- array(tmp[,2],c(M,N))
	yB <- array(tmp[,3],c(M,N))

	tmpR <- rescale(yR,f,boxKernel())	
	tmpG <- rescale(yG,f,boxKernel())	
	tmpB <- rescale(yB,f,boxKernel())	
	mydata <- data.frame(cbind(as.vector(tmpR[,]),as.vector(tmpG[,]),as.vector(tmpB[,])))
	
	for(k in 1:3)mydata[,k] <- histstretch(mydata[,k])

	cellcentre.offset <- A@grid@cellcentre.offset - 0.5*ps + 0.5*ps/f
	mygrid <- GridTopology(cellcentre.offset, ps/f, ceiling(c(M,N)*f))
	

	Adisp <- SpatialGridDataFrame(mygrid, mydata, proj4string = A@proj4string)
	
	windows()
	image(Adisp,red=nR,green=nG,blue=nB,axes=TRUE)
	title(main=title_str)
}

# Remove overlapping points (due to tile overlap effects)
# Clean spatially overlapping points
clean_overlap_old <- function(Blobs){

	pt <- 1
	val <- Blobs

	while(pt<=nrow(val))
	{
		x0 <- val[pt,1]
		y0 <- val[pt,2]
		r0 <- 0.5*val[pt,3]

		xall <- val[,1]
		yall <- val[,2]
		rall <- 0.5*val[,3]
		
		dall <- sqrt((xall-x0)^2+(yall-y0)^2)

		ind <- which(dall+rall<r0)
		
		if(length(ind)>1)
		{
			tmp <- max(val[ind,]$magn)
			ind0 <- which(val[ind,]$magn==tmp)
			
			ind0 <- ind[ind0]
			val[pt,] <- val[ind0,]
			
			ind <- ind[ind!=ind0]
			val <- val[-ind,]
		}else pt <- pt+1
	}

	return(val)
}


clean_overlap <- function(Blobs,r_thr){
	ind <- order(Blobs$magn,decreasing=TRUE)
	Blobs <- Blobs[ind,]

	val <- Blobs
	weak <- array(0,nrow(val))

	pt <- 0
	while(pt<nrow(val))
	{
		pt <- pt +1
		if(weak[pt]==1) next

		x0 <- val[pt,1]
		y0 <- val[pt,2]
		r0 <- 0.5*val[pt,3]

		xall <- val[,1]
		yall <- val[,2]
		rall <- 0.5*val[,3]

		dall <- sqrt((xall-x0)^2+(yall-y0)^2)

		ind <- which((dall<=rall+r0)&(rall>r0)&(weak==0))

		subind <- which(ind<=pt)

		if(length(subind)>0) ind <- ind[-subind]

		if(length(ind)==0) next

		ind2 <- which(val$magn[ind]<val$magn[pt])
		if(length(ind2)>0){
			weak[ind[ind2]] <- 1
		}
	}

	val <- val[weak==0,]

	weak <- array(0,nrow(val))

	pt <- 0
	while(pt<nrow(val))
	{
		pt <- pt +1
		if(weak[pt]==1) next

		x0 <- val[pt,1]
		y0 <- val[pt,2]
		r0 <- 0.5*val[pt,3]

		xall <- val[,1]
		yall <- val[,2]
		rall <- 0.5*val[,3]

		dall <- sqrt((xall-x0)^2+(yall-y0)^2)

		#ind <- which((dall<=pmin(rall,r0))&(weak==0))
		ind <- which((dall<=0.5*(rall+r0))&(weak==0))

		subind <- which(ind<=pt)

		if(length(subind)>0) ind <- ind[-subind]

		if(length(ind)==0) next
		
		ind2 <- which(val$magn[ind]<val$magn[pt])
		if(length(ind2)>0){
			weak[ind[ind2]] <- 1
		}
	}

	val <- val[weak==0,]

	return(val)
}


clean_overlap_new <- function(Blobs){

	image(MSdisp,red=nR,green=nG,blue=nB,axes=TRUE)
	title(main="phase1")

	ind <- order(Blobs$d,decreasing=TRUE)
	Blobs <- Blobs[ind,]

	
	val <- Blobs
	
	weak <- array(0,nrow(val))

	pt <- 0
	while(pt<nrow(val))
	{
		pt <- pt +1
		if(weak[pt]==1) next
		
		display_all_blobs(val[pt,],"white")
		
		x0 <- val[pt,1]
		y0 <- val[pt,2]
		r0 <- 0.5*val[pt,3]

		xall <- val[,1]
		yall <- val[,2]
		rall <- 0.5*val[,3]
		
		dall <- sqrt((xall-x0)^2+(yall-y0)^2)

#		ind <- which((dall<=r0)&(weak==0))
		ind <- which((dall<r0)&(weak==0))
		
		subind <- which(ind==pt)
		if(length(subind)>0) ind <- ind[-subind]
		
		if(length(ind)==0) next

		pt2 <- ind[1]
		#display_all_blobs(val[pt2,],"green")
		
		if(val$magn[pt2]>val$magn[pt]){
			weak[pt] <- 1
		} else {
			weak[pt2] <-1
			pt <- pt-1
		}
				
	}


	display_all_blobs(val[weak==0,],"blue") 
	display_all_blobs(val[weak==1,],"white")
	#text(Blobs$x,Blobs$y,labels=1:nrow(Blobs),pos=4)

	return(val[weak==0,])
}


clean_overlap_v3 <- function(Blobs){

	image(MSdisp,red=nR,green=nG,blue=nB,axes=TRUE)
	title(main="phase1")

	ind <- order(Blobs$d,decreasing=TRUE)
	Blobs <- Blobs[ind,]

	
	val <- Blobs
	weak <- array(0,nrow(val))

	pt <- 0
	while(pt<nrow(val))
	{
		pt <- pt +1
		
		#display_all_blobs(val[pt,],"white")
		
		x0 <- val[pt,1]
		y0 <- val[pt,2]
		r0 <- 0.5*val[pt,3]

		xall <- val[,1]
		yall <- val[,2]
		rall <- 0.5*val[,3]
		
		dall <- sqrt((xall-x0)^2+(yall-y0)^2)

#		ind <- which((dall+rall<=r0)&(weak==0))
#		ind <- which((dall<=rall+r0)&(weak==0))
		ind <- which((dall<=r0)&(weak==0))
		
		subind <- which(ind<=pt)
		
		if(length(subind)>0) ind <- ind[-subind]
		
		if(length(ind)==0) next

		#display_all_blobs(val[ind,],"green")
		
		if(max(val$magn[ind])>val$magn[pt]){
			weak[pt] <- 1
			next
		}
				
	}


	display_all_blobs(val[weak==0,],"blue") 
	display_all_blobs(val[weak==1,],"white")

	val <- val[weak==0,]

	#windows()
	image(MSdisp,red=nR,green=nG,blue=nB,axes=TRUE)
	title(main="phase2")

#	display_all_blobs(val,"blue") 

	weak <- array(0,nrow(val))
	
	pt <- 0
	while(pt<nrow(val))
	{
		pt <- pt+1
		if(weak[pt]==1) next
		
		display_all_blobs(val[pt,],"white")
		
		x0 <- val[pt,1]
		y0 <- val[pt,2]
		r0 <- 0.5*val[pt,3]

		xall <- val[,1]
		yall <- val[,2]
		rall <- 0.5*val[,3]
		
		dall <- sqrt((xall-x0)^2+(yall-y0)^2)

#		ind <- which((dall<=r0)&(weak==0)) 
#		ind <- which((dall<=r0)&(weak==0))
		ind <- which((dall<=r0)&(weak==0)) 

		subind <- which(ind<=pt)

		if(length(subind)>0) ind <- ind[-subind]

		if(length(ind)==0) next

		display_all_blobs(val[ind,],"green")
		text(val$x[ind],val$y[ind],labels=ind,pos=4,col="white")
		
		ind2 <- which(val$magn[ind]<val$magn[pt])
		if(length(ind2)>0){
			weak[ind[ind2]] <- 1
			#display_all_blobs(val[ind[ind2],],"yellow")
		}

	}


	
	
	display_all_blobs(val[weak==0,],"blue") 
	display_all_blobs(val[weak==1,],"white")
	#text(Blobs$x,Blobs$y,labels=1:nrow(Blobs),pos=4)

	
	return(val[weak==0,])
}


clean_overlap_v4 <- function(Blobs){

	ind <- order(Blobs$d,decreasing=TRUE)
	Blobs <- Blobs[ind,]

	pt <- 0
	val <- Blobs

	while(pt<nrow(val))
	{
		pt <- pt+1
	
		display_all_blobs(val[pt,],"blue") 
	
		x0 <- val[pt,1]
		y0 <- val[pt,2]
		r0 <- 0.5*val[pt,3]

		rest <- val[(pt+1):nrow(val),]
		
		xall <- rest[,1]
		yall <- rest[,2]
		rall <- 0.5*rest[,3]
		
		dall <- sqrt((xall-x0)^2+(yall-y0)^2)

		# Fully contained smaller objects
		ind <- which(rall+dall<=r0) 
#		ind <- which(dall<=r0) 
		
		if(length(ind)==0) next
		
		ind <- ind+pt
		
		if(max(val$magn[ind])<=val$magn[pt]){

			next
		}

		display_all_blobs(val[pt,],"white") 
		val <- val[-pt,]
		
		pt <- pt - 1
		
	}

	
	
	text(Blobs$x,Blobs$y,labels=1:nrow(Blobs),pos=4)

	
	return(val[weak==0,])
}

# r1 is scalar, whereas r2 and d are vectors of same lengths
circle_overlap <- function(r1,r2,d){

	x <- (d*d + r1*r1 - r2*r2)/(2*d)


	a <- (d+r2-r1)*(d+r1-r2)*(r1+r2-d)*(r1+r2+d)
	a[a <0] <- NA
	a <- sqrt(a)/d
	
	ind <- which(d>=r1+r2)
	
	if(length(ind)>0)a[ind] <- 0
	
	ind <- which(is.na(a))
	if(length(ind)>0) a[ind] <- 0
	
	val <- r1*r1* acos((d*d+r1*r1-r2*r2)/(2*d*r1)) + r2*r2* acos((d*d+r2*r2-r1*r1)/(2*d*r2)) - 0.5 * a
	
	
	val <- val/sqrt((pi*r1*r1)*(pi*r2*r2))
	
	return(val)
}

clean_overlap_tree_mask <- function(Blobs){

	ind <- order(Blobs$d,decreasing=TRUE)
	Blobs <- Blobs[ind,]

	val <- Blobs
	
	weak <- array(0,nrow(val))
	
	#windows()
	#image(MSdisp,red=nR,green=nG,blue=nB,axes=TRUE)

	for(pt in 1:nrow(val)){
	
		if(weak[pt]==1) next
	
		#display_all_blobs(val[pt,],"blue") 
	
		x1 <- val$x[pt]
		y1 <- val$y[pt]
		r1 <- 0.5*val$d[pt]

		rest <- val
		#rest[pt,] <- NA

		x2 <- rest$x
		y2 <- rest$y
		r2 <- 0.5*rest$d

		d <- sqrt((x2-x1)^2+(y2-y1)^2)
		
		# measure area overlap
		#ar_ov <- circle_overlap(r1,r2,d)

		# Fully contained objects
		ind <- which((r2+d<=r1)&(weak==0))

		if(length(ind)>0){

			#ind3 <- ind[which(Blobs$d[ind]>2*mean(ps.ms))]
			ind3 <- ind
			if(length(ind3)>0){
				#display_all_blobs(Blobs[ind3,],"green")
				
				#if((max(val$magn[ind3])>val$magn[pt])&(max(val$ndvi[ind3])>val$ndvi[pt])){
				#if((max(val$magn[ind])>val$magn[pt])&(max(val$ndvi[ind])>val$ndvi[pt])){
				if((max(val$magn[ind])>val$magn[pt])){
					#display_all_blobs(val[pt,],"white") 
					weak[pt] <- 1
				} else{
					small_blobs <- Blobs[ind,]
					
					ind2 <- ind[which(small_blobs$magn<Blobs$magn[pt])]
					#display_all_blobs(Blobs[ind2,],"white")
					weak[ind2] <- 1
				}
			}
		}
		
		# objects crossing outlines 
		ind <- which((d>r1-r2)&(d<=r1+r2)&(weak==0))

		if(length(ind)>0){

			ind3 <- ind[which(Blobs$d[ind]>2*mean(ps.ms))]
			if(length(ind3)>0){
				#display_all_blobs(Blobs[ind3,],"green")
				
				if((max(val$magn[ind3])>val$magn[pt])&(max(val$ndvi[ind3])>val$ndvi[pt])){
				#if((max(val$magn[ind])>val$magn[pt])&(max(val$ndvi[ind])>val$ndvi[pt])){
				#if((max(val$magn[ind])>val$magn[pt])){
					#display_all_blobs(val[pt,],"white") 
					weak[pt] <- 1
				} else{
					small_blobs <- Blobs[ind,]
					
					ind2 <- ind[which(small_blobs$magn<Blobs$magn[pt])]
					#display_all_blobs(Blobs[ind2,],"white")
					weak[ind2] <- 1
				}
			}
		}
		
	}

	#image(MSdisp,red=nR,green=nG,blue=nB,axes=TRUE)
	#display_all_blobs(val[weak==1,],"white") 
	#display_all_blobs(val[weak==0,],"blue") 
	#text(Blobs$x,Blobs$y,labels=1:nrow(Blobs),pos=4)
	
	return(val[weak==0,])
}

blobs_contained <- function(Blobs, Large_blobs){

	val <- Blobs
	
	weak <- array(0,nrow(val))
	
	#windows()
	#image(MSdisp,red=nR,green=nG,blue=nB,axes=TRUE)

	for(pt in 1:nrow(val)){
	
		if(weak[pt]==1) next
	
		#display_all_blobs(val[pt,],"blue") 
	
		x1 <- val$x[pt]
		y1 <- val$y[pt]
		r1 <- 0.5*val$d[pt]

		rest <- Large_blobs

		x2 <- rest$x
		y2 <- rest$y
		r2 <- 0.5*rest$d

		d <- sqrt((x2-x1)^2+(y2-y1)^2)
		
		# measure area overlap
		#ar_ov <- circle_overlap(r1,r2,d)

		# Fully contained objects
		ind <- which(d<=r1+r2)

		if(length(ind)>0){
			weak[pt] <- 1
			#display_all_blobs(Blobs[pt,],"white")
		}
	}

	return(val[weak==0,])
}


clean_cocentric_blobs <- function(Blobs){

	val <- Blobs
	
	weak <- array(0,nrow(val))
	
	#windows()
	#image(MSdisp,red=nR,green=nG,blue=nB,axes=TRUE)

	for(pt in 1:nrow(val)){
	
		if(weak[pt]==1) next
	
		#display_all_blobs(val[pt,],"blue") 
	
		x1 <- val$x[pt]
		y1 <- val$y[pt]
		r1 <- 0.5*val$d[pt]

		rest <- Blobs

		x2 <- rest$x
		y2 <- rest$y
		r2 <- 0.5*rest$d

		d <- sqrt((x2-x1)^2+(y2-y1)^2)
		
		# measure area overlap
		#ar_ov <- circle_overlap(r1,r2,d)

		# Fully contained objects
		ind <- which(d<=0.125*(r1+r2))

		if(length(ind)>0){
			if(max(Blobs$magn[ind])>Blobs$magn[pt]){
				weak[pt] <- 1
				#display_all_blobs(Blobs[pt,],"white")
			}
		}
	}

	return(val[weak==0,])
}


geo_transform <- function(Blobs,dx,dy){

	Blobs$x <- Blobs$x + dx
	Blobs$y <- Blobs$y + dy

	return(Blobs)
}

valild_match <- function(ind_a,ind_b){

	a <- Blobs_A[ind_a,] 
	b <- Blobs_B[ind_b,] 
	
	#val <- (abs(a$d-b$d)<0.25*(a$d+b$d)/2) & (abs(a$ndvi-b$ndvi)<0.25) & (abs(a$magn-b$magn)<0.5)
	val_d <- (abs(a$d-b$d)<0.25*max(a$d,b$d)) 
	val_ndvi <- (abs(a$ndvi-b$ndvi)<0.25) 
	val_magn <- (abs(a$magn-b$magn)<0.5)

	val <- val_d 
	
	val
}

qual_match <- function(ind_a,ind_b){

	a <- Blobs_A[ind_a,] 
	b <- Blobs_B[ind_b,] 

	b <- geo_transform(b,dx,dy)
	
	val_r <- sqrt((a$x - b$x)^2 + (a$y - b$y)^2)
	val_d <- min(a$d,b$d) / max(a$d,b$d)
		
	val_ndvi <- min(a$ndvi,b$ndvi)/max(a$ndvi,b$ndvi)
	val_magn <- min(a$magn,b$magn)/max(a$magn,b$magn)

	
	c(val_r,val_d,val_ndvi,val_magn)
}



interpolate_point_v3 <- function(i,j,k,detHes,ndvi,ps,xy){

	d <- dim(detHes)
	M <- d[1]
	N <- d[2]

	val <- array(0,5)

	ind1 <- (j-1)*M + i

	nn <- 1
	k1 <- max(k - nn,1)
	k2 <- min(k + nn,Ns)
	x <- tarr[k1:k2]
	#x <- sqrt(tarr[k1:k2])
	y <- detHes[i,j,k1:k2]

	a <- parabol_interp_vect(t(x),t(y))

	#t0 <- a[3]^2
	t0 <- a[3]
	
	val[4] <- a[1]
	val[3] <- 2*sqrt(2*t0)*mean(ps)

	# Position
	i1 <- max(i - nn,1)
	i2 <- min(i + nn,M)

	i_arr <- i1:i2
	j_arr <- rep(j,length(i_arr))

	pn_arr <- (j_arr-1)*M + i_arr

	x <- xy[pn_arr,1]

	y <- detHes[i1:i2,j,k]

	a <- parabol_interp(x,y)

	val[1] <- a[3]

	j1 <- max(j - nn,1)
	j2 <- min(j + nn,N)

	j_arr <- j1:j2
	i_arr <- rep(i,length(j_arr))
	
	pn_arr <- (j_arr-1)*M + i_arr

	x <- xy[pn_arr,2]

	y <- detHes[i,j1:j2,k]

	a <- parabol_interp(x,y)

	val[2] <- a[3]

	x0 <- val[1]
	y0 <- val[2]
	
	#ndvi test
	xyp <- xy
	xyp[,1] <- xyp[,1] - x0
	xyp[,2] <- xyp[,2] - y0
	
	dp <- sqrt(rowSums(xyp^2))
	
	indxy <- which(dp <=val[3]/2)
	val[5] <- mean(ndvi[indxy])

	return(val)
}


interpolate_point_v2 <- function(i,j,k,detHes,ndvi,ps,xy){
	d <- dim(detHes)
	M <- d[1]
	N <- d[2]

	val <- array(0,5)

	ind1 <- (j-1)*M + i

	#nn <- 1
	#k1 <- max(k - nn,1)
	#k2 <- min(k + nn,Ns)
	#x <- tarr[k1:k2]
	#y <- detHes[i,j,k1:k2]
	#a <- parabol_interp(x,y)

	x <- sqrt(tarr[2:Ns])
	y <- detHes[i,j,2:Ns]
	
	t0 <- determine_max(x,y)
	t0 <- t0^2
	val[4] <- detHes[i,j,k]
	val[3] <- 2*sqrt(2*t0)*mean(ps)

	# Interpolation for x coordinate

	# Check if function is decreasing
	y <- detHes[,j,k]

	if(Debug){
		z <- c(0,diff(y))
		z2 <- c(0,diff(z))

		plot(y)
		lines(y)
		lines(z,col="blue")
		lines(z2,col="green")
		abline(h=0,lty=2)
	}
	
	nn_max <- 4
	nn1 <- min(c(nn_max,i-1))

	s1 <- sign(diff(y[i:(i-nn1)]))

	ind <- which(s1>0)
	if(length(ind)==0) ind <- nn1

	i1 <- i-ind[1]+1

	nn2 <- min(c(nn_max,M-i))

	s2 <- sign(diff(y[i:(i+nn2)]))

	ind <- which(s2>0)
	if(length(ind)==0) ind <- nn2+1

	i2 <- i+ind[1] - 1

	if(Debug){
		abline(v=i,lty=2,col="blue")
		abline(v=i1,lty=2,col="blue")
		abline(v=i2,lty=2,col="blue")
	}
	
	pn <- (j-1)*M + i

	x <- i1:i2 - i
	y <- detHes[i1:i2,j,k]
	
	a <- determine_max(x,y)
	
	val[1] <- a*ps[1]+xy[pn,1]

	# Interpolation for y coordinate
	
	# Check if function is decreasing
	y <- detHes[i,,k]
	
	if(Debug){
		z <- c(0,diff(y))
		z2 <- c(0,diff(z))

		plot(y)
		lines(y)
		lines(z,col="blue")
		lines(z2,col="green")
		abline(h=0,lty=2)
	}
	nn_max <- 4
	nn1 <- min(c(nn_max,j-1))

	s1 <- sign(diff(y[j:(j-nn1)]))

	ind <- which(s1>0)
	if(length(ind)==0) ind <- nn1

	j1 <- j-ind[1]+1

	nn2 <- min(c(nn_max,N-j))

	s2 <- sign(diff(y[j:(j+nn2)]))

	ind <- which(s2>0)
	if(length(ind)==0) ind <- nn2+1

	j2 <- j+ind[1] - 1

	if(Debug){
		abline(v=j,lty=2,col="blue")
		abline(v=j1,lty=2,col="blue")
		abline(v=j2,lty=2,col="blue")
	}

	
	x <- j1:j2 - j
	y <- detHes[i,j1:j2,k]
	
	a <- determine_max(x,y)
	
	val[2] <- a*ps[2]+xy[pn,2]

	x0 <- val[1]
	y0 <- val[2]
	
	# Compute average ndvi value inside
	xy[,1] <- xy[,1] - x0
	xy[,2] <- xy[,2] - y0
	
	dp <- sqrt(rowSums(xy^2))
	
	indxy <- which(dp <=val[3]/2)
	val[5] <- mean(ndvi[indxy])

	return(val)
}

detect_blobs_all <- function(A, tarr){
	d <- A@grid@cells.dim
	
	P <- array(A$ndvi,d)

	detHes <- DOH_space(P,tarr)

	# Delete negative values
	detHes[detHes<0] <- 0

	# Remove NAs and NaNs
	ind <- which(is.na(detHes))
	detHes[ind] <- 0
	
	# Find local extrema in scale space
	pts <- local_max_3D(detHes)	
		
	ptsdf <- data.frame(pts)
	names(ptsdf) <- c("col","row","s","magn","sort")

	Special_points <- array(0,c(0,5))

	for(k in Ns:2)
	{
		ind <- which(ptsdf$s==k)

		Npts <- length(ind)

		if(Npts>0)
		{
			tmpdf <- ptsdf[ind,]
						
			tmp <- tmpdf$magn
			inds <- order(tmp,decreasing=TRUE)
			
			tmpdf <- tmpdf[inds,]

			#interpolate the point
			for(pt in 1:Npts)
			{
				i <- tmpdf$col[pt]
				j <- tmpdf$row[pt]

				Spec_pt <- interpolate_point(i,j,k, detHes,A)

				Special_points <- rbind(Special_points,Spec_pt,deparse.level = 0)
			}
		}
	}

	Special_points  <- data.frame(Special_points)
	names(Special_points) <- c("x","y","d","magn","ndvi")

	return(Special_points)
}


detect_blobs <- function(A, tarr, ndvi_thresh, magn_thresh){
	
	d <- A@grid@cells.dim
	
	P <- array(A$ndvi,d)
	
	detHes <- DOH_space(P,tarr)

	# Delete negative values
	detHes[detHes<0] <- 0

	# Remove NAs and NaNs
	ind <- which(is.na(detHes))
	detHes[ind] <- 0
	
	# Find local extrema in scale space
	pts <- local_max_3D(detHes)	
		
	ptsdf <- data.frame(pts)
	names(ptsdf) <- c("col","row","s","magn","sort")

	Special_points <- array(0,c(0,5))

	for(k in Ns:2)
	{
		ind <- which(ptsdf$s==k)

		Npts <- length(ind)

		if(Npts==0) next
		tmpdf <- ptsdf[ind,]
					
		tmp <- tmpdf$magn
		inds <- order(tmp,decreasing=TRUE)
		
		tmpdf <- tmpdf[inds,]

		#interpolate the point
		for(pt in 1:Npts)
		{
			i <- tmpdf$col[pt]
			j <- tmpdf$row[pt]

			Spec_pt <- interpolate_point(i,j,k, detHes,A)

			Special_points <- rbind(Special_points,Spec_pt,deparse.level = 0)
		}
	
	}

	Special_points  <- data.frame(Special_points)
	names(Special_points) <- c("x","y","d","magn","ndvi")


	# Apply magnitude	& ndvi thresholds
	ind <- which((Special_points$ndvi>=ndvi_thresh)&(Special_points$magn>=magn_thresh))
	
	Special_points <- Special_points[ind,]

	# Clean spatially overlapping points
	Best_points <- Special_points
	#Best_points <- clean_overlap_v3(Special_points)

	return(Best_points)
}

# Using C++
detect_blobs_v2 <- function(P, ps, xTL, yTL, tarr, magn_thresh){

	detHes <- DOH_space(P,tarr)

	# Remove NAs and NaNs
	detHes[which(is.na(detHes))] <- 0

	# Find local extrema in scale space
	#pts <- local_max_3D(detHes)
	d <- dim(detHes)
	pts <- local_max_3d(detHes,d[1],d[2],d[3])

	magn_arr <- detHes[pts[,1] + (pts[,2]-1)*d[1] + (pts[,3]-1)*d[1]*d[2]]
	ind <- which(magn_arr>=magn_thresh)

	pts <- pts[ind, , drop=F]

	z <- interpolate_blob(detHes,d[1],d[2],d[3], pts, ps, xTL, yTL, tarr)
	Special_points  <- data.frame(z)
	names(Special_points) <- c("x","y","d","magn")

	
	# Apply magnitude	& ndvi thresholds
	ind <- which(Special_points$magn>=magn_thresh)
	
	Special_points <- Special_points[ind,]

	# Clean spatially overlapping points
	Best_points <- Special_points
	#Best_points <- clean_overlap_v3(Special_points)

	return(Best_points)
}


# Using C++
detect_blobs_v3 <- function(P, n_col, n_row, ps, xTL, yTL, tarr, magn_thresh){

	n_s <- length(tarr)

	detHes <- DOH_spaceC(P,n_col,n_row,tarr)

	# Remove NAs and NaNs
	detHes[which(is.na(detHes))] <- 0

	# Find local extrema in scale space
	#pts <- local_max_3D(detHes)
	#d <- dim(detHes)
	d <- c(n_col,n_row,length(tarr))
	pts <- local_max_3d(detHes,n_col,n_row,n_s)

	# Apply magnitude threshold
	magn_arr <- detHes[pts[,1] + (pts[,2]-1)*n_col + (pts[,3]-1)*n_col*n_row]
	ind <- which(magn_arr>=magn_thresh)
	pts <- pts[ind, , drop=F]

	z <- interpolate_blob(detHes,n_col,n_row,n_s, pts, ps, xTL, yTL, tarr)
	Special_points  <- data.frame(z)
	names(Special_points) <- c("x","y","d","magn")

	return(Special_points)
}





detect_gaps <- function(A, tarr, ndvi_thresh, magn_thresh){
	
	d <- A@grid@cells.dim
	
	P <- array(A$ndvi,d)
	
	detHes <- DOH_space(P,tarr)

	# Delete negative values
	#detHes[detHes<0] <- 0

	# Remove NAs and NaNs
	ind <- which(is.na(detHes))
	detHes[ind] <- min(detHes,na.rm=TRUE)
	
	# Find local extrema in scale space
	pts <- local_max_3D(detHes)	
		
	ptsdf <- data.frame(pts)
	names(ptsdf) <- c("col","row","s","magn","sort")

	Special_points <- array(0,c(0,5))

	for(k in Ns:2)
	{
		ind <- which(ptsdf$s==k)

		Npts <- length(ind)

		if(Npts==0) next
		tmpdf <- ptsdf[ind,]
					
		tmp <- tmpdf$magn
		inds <- order(tmp,decreasing=TRUE)
		
		tmpdf <- tmpdf[inds,]

		#interpolate the point
		for(pt in 1:Npts)
		{
			i <- tmpdf$col[pt]
			j <- tmpdf$row[pt]

			Spec_pt <- interpolate_point(i,j,k, detHes,A)

			Special_points <- rbind(Special_points,Spec_pt,deparse.level = 0)
		}
	
	}

	Special_points  <- data.frame(Special_points)
	names(Special_points) <- c("x","y","d","magn","ndvi")


	# Apply magnitude	& ndvi thresholds
	#ind <- which((Special_points$ndvi<ndvi_thresh)&(Special_points$magn>=magn_thresh))
	ind <- which(Special_points$magn>=magn_thresh)
	
	Special_points <- Special_points[ind,]

	# Clean spatially overlapping points
	Best_points <- Special_points
	#Best_points <- clean_overlap_v3(Special_points)

	return(Best_points)
}


	
# Custom GLCM functions

global_glcm <- function(xyz,pix,lag_x=lag_x,lag_y=lag_y)
{
	npix <- nrow(xyz)

	zr <- range(xyz[,3])
	zq <- xyz[,3]

	GM <- array(0,c(ql,ql))

	count <- 0

	xy <- xyz[,1:2]

	len_arr <- array(0,npix)
	
	for(i in 1:npix)
	{
		tmp <- xy[,]
		tmp[,1] <- tmp[,1]- xy[i,1]
		tmp[,2] <- tmp[,2]- xy[i,2]
		
		ind <- which((tmp[,1]>=(lag_x-0.5)*pix) & (tmp[,1]<(lag_x+0.5)*pix) & (tmp[,2]>=(lag_y-0.5)*pix) & (tmp[,2]<(lag_y+0.5)*pix))
		nn <- length(ind)
		len_arr[i] <- nn
		
		count <- count + nn
		
		GM[zq[i],zq[ind]] <- GM[zq[i],zq[ind]]+1
	}

	return(GM/count)
}

global_glcm_isotropic <- function(xyz,pix,lag=1,width=0.5)
{
	npix <- nrow(xyz)

	zr <- range(xyz[,3])
	zq <- xyz[,3]

	GM <- array(0,c(ql,ql))

	count <- 0

	xy <- xyz[,1:2]

	len_arr <- array(0,npix)
	
	for(i in 1:npix)
	{
		tmp <- xy[,]
		tmp[,1] <- tmp[,1]- xy[i,1]
		tmp[,2] <- tmp[,2]- xy[i,2]
		d <- sqrt(rowSums(tmp^2))
		
		ind <- which((d>=(lag-0.5*width)*pix) & (d<(lag+0.5*width)*pix))
		nn <- length(ind)
		if(nn>0)
		{
			len_arr[i] <- nn
			count <- count + nn
			for(j in 1:ql) GM[zq[i],j] <- GM[zq[i],j] + sum(zq[ind]==j)
		}
	}

	return(GM/count)
}

global_glcm_polar <- function(xyz,pix,lag_r=lag_r,dr=dr,lag_phi=lag_phi,dphi=dphi)
{
	npix <- nrow(xyz)

	zr <- range(xyz[,3])
	zq <- xyz[,3]

	GM <- array(0,c(ql,ql))

	count <- 0

	xy <- xyz[,1:2]

	len_arr <- array(0,npix)
	
	for(i in 1:npix)
	{
		x <- xy[,1]- xy[i,1]
		y <- xy[,2]- xy[i,2]
		
		r <- sqrt(x^2+y^2)
		phi <- atan2(y,x)
		ind <- phi<0
		phi[ind] <- phi[ind]+pi
		
		ind <- which((r>=(lag_r-dr)*pix) & (r<(lag_r+dr)*pix) & (phi>=(lag_phi-dphi)) & (phi<(lag_phi+dphi)))

		nn <- length(ind)
		len_arr[i] <- nn
		
		count <- count + nn
		
		for(j in 1:nn)
		{
			a <- min(zq[i],zq[ind[j]])
			b <- max(zq[i],zq[ind[j]])
			
			GM[a,b] <- GM[a,b]+1
		}
	}
	
	# symmetrize the matrix
	for(a in 1:(ql-1))
	for(b in (a+1):ql)
	{
		GM[b,a]<-GM[a,b]
	}
	
	return(GM/sum(GM))
}


global_glcm_fast <- function(xyz,pix,lag_r=lag_r,lag_dir)
{
	npix <- nrow(xyz)

	zr <- range(xyz[,3])
	zq <- xyz[,3]

	GM <- array(0,c(ql,ql))

	count <- 0

	xy <- xyz[,1:2]

	len_arr <- array(0,npix)
	

	xy_sp <- SpatialPixelsDataFrame(xy, data.frame(xyz[,3]), proj4string = CRS(as.character(NA)), round = NULL, grid = NULL)
	
	xy_shift <- xy
	
	if(lag_dir==1)
	{
		xy_shift[,1] <- xy[,1] + lag_r*pix
	}

	xy2_sp <- SpatialPixelsDataFrame(xy_shift, data.frame(xyz[,3]), proj4string = CRS(as.character(NA)), round = NULL, grid = NULL)

	over(xy_sp,xy2_sp)
	
	for(i in 1:nrow(xy)){
		xy_shift[i,]
	}
	
	for(i in 1:npix)
	{
		x0 <- xy[i,1]
		y0 <- xy[i,1]
		
		switch()
	
		x <- xy[,1]- xy[i,1]
		y <- xy[,2]- xy[i,2]
		
		r <- sqrt(x^2+y^2)
		phi <- atan2(y,x)
		ind <- phi<0
		phi[ind] <- phi[ind]+pi
		
		ind <- which((r>=(lag_r-dr)*pix) & (r<(lag_r+dr)*pix) & (phi>=(lag_phi-dphi)) & (phi<(lag_phi+dphi)))

		nn <- length(ind)
		len_arr[i] <- nn
		
		count <- count + nn
		
		for(j in 1:nn)
		{
			a <- min(zq[i],zq[ind[j]])
			b <- max(zq[i],zq[ind[j]])
			
			GM[a,b] <- GM[a,b]+1
		}
	}
	
	# symmetrize the matrix
	for(a in 1:(ql-1))
	for(b in (a+1):ql)
	{
		GM[b,a]<-GM[a,b]
	}
	
	return(GM/sum(GM))
}


glcm_mean <- function(A)
{
	return(sum(row(A)*A))
}

glcm_variance <- function(A)
{
	m <- glcm_mean(A)
	
	return(sum(A*(row(A)-m)^2))
}

glcm_energy <- function(A)
{
	return(sqrt(sum(A^2)))
}

glcm_entropy <- function(A)
{
	val <- -A*log(A)

	val[A==0] <- 0
	
	return(sum(val))
}
	

pan_profile <- function(phic,rc,Smax,r,ps,psi){
	dpsi_corr <- array(dpsi,length(rc))
	dpsi_corr <- dpsi_corr * mean(ps)/rc
	indc <- which((abs(phic-psi)<=dpsi_corr)&(rc <= Smax)&(rc>r))
	rs   <- rc[indc]
	indr <- order(rs)
	indc <- indc[indr]
	return(indc)
}	

height_shadow <- function(Blobs,A,psi){

	tree_height <- array(0,nrow(Blobs))

	for(id in 1:nrow(Blobs)){
		pan <- A$band1
		xy <- coordinates(A)
		ps <- A@grid@cellsize

		val <- 0
		xy0 <- as.numeric(Blobs[id,1:2])
		r <- as.numeric(0.5*Blobs$d[id])
		
		xyc <- xy[,]
		xyc[,1] <- xyc[,1] - xy0[1]
		xyc[,2] <- xyc[,2] - xy0[2]
		rc <- sqrt(rowSums(xyc^2))
		phic <- atan2(xyc[,2],xyc[,1])

		hmax <- r*4
		Smax <- hmax*tan(theta_sun)

		indc <- pan_profile(phic,rc,Smax,r,ps,psi)

		if(length(indc)<4){
			next
		}
		
		rs <- rc[indc]
		pans <- pan[indc] - shad_thr
		phis <- phic[indc]

		if(all(pans>0)) next
		
		# Delete positive points at the small r
		if(pans[1]>0){
			tmp <- rle(pans>0)
			pt <- tmp$lengths[1]
			pans[1:pt] <- pans[pt+1]
			pan[indc[1:pt]] <- pans[pt+1]
		}
		
		#Delete short spikes in data
		tmp <- rle(pans>0)

		ind_spike <- which(tmp$lengths<3)
		while(length(ind_spike)>0){

			t2 <- sum(tmp$lengths[1:ind_spike[1]])
			t1 <- t2 - tmp$lengths[ind_spike[1]] + 1

			if(t2==length(pans)) break
			
			if(t1>1)v1 <- pans[t1-1] else v1 <- pans[t2+1]
			if(t2 < length(pans))v2 <- pans[t2+1]
			
			pans[t1:t2] <- mean(v1,v2)
			pan[indc[t1:t2]] <- mean(v1,v2)
			tmp <- rle(pans>0)
			ind_spike <- which(tmp$lengths<4)
		}

		
		# Decide if region growing is required
		Grow <- FALSE
		if(all(pans<0))Grow <- TRUE
		
		while(Grow){
			Smax <- Smax + mean(ps)
			indc <- pan_profile(phic,rc,Smax,r,ps,psi)
			rs <- rc[indc]
			pans <- pan[indc] - shad_thr
			phis <- phic[indc]
			l <- length(pans)

			if((all(pans[(l-2):l]>0))|(Smax>20*hmax)) Grow <- FALSE
		}
		
		if(Smax>=20*hmax) next
		
		smoothingSpline <- smooth.spline(rs, pans, spar=0.2)

		pans <- smoothingSpline$y
		
		p1 <- min(pans)
		if(p1>0){
		# there is no shadow associated with this object
			next
		}
		
		pt0 <- which(pans==p1)
		r0 <- rs[pt0]
		
		pt_in <- pt0 - 1 + rle(pans[pt0:length(pans)]<0)$lengths[1]
		pt_out <- pt_in + 1
		
		r1 <- rs[pt_in]
		r2 <- rs[pt_out]
		
		p1 <- pans[pt_in]
		p2 <- pans[pt_out]
		
		if(p1<0){
			rtip <- r1 + (r2-r1)*(0-p1)/(p2-p1)
		}else{
			rtip=NA
		}

		xy_tip <- xy0 + rtip*c(cos(psi),-sin(psi))

		Sprime <- sqrt(sum((xy_tip-xy0)^2))
		h1 <- Sprime*cos(psi) / (0.5*tan(theta_sat)*cos(alpha_sat)-tan(theta_sun)*cos(alpha_sun))
		#h2 <- Sprime*sin(psi-pi/2) / (abs(tan(theta_sun)*sin(-sun_az))+abs(0.5*tan(theta_sat)*sin(-sat_az)))

		#val <- mean(h1,h2)
		tree_height[id] <- h1
	}
	
	return(tree_height)
}


shadow_tip<-function(xy0,xy,r,psi){

		xyc <- xy[,]
		xyc[,1] <- xyc[,1] - xy0[1]
		xyc[,2] <- xyc[,2] - xy0[2]
		rc <- sqrt(rowSums(xyc^2))
		phic <- atan2(xyc[,2],xyc[,1])

		hmax <- r*4
		Smax <- hmax*tan(theta_sun)

		indc <- pan_profile(phic,rc,Smax,r,ps,psi)

		if(length(indc)<4){
			return(xy0)
		}
		
		rs <- rc[indc]
		pans <- pan[indc] - shad_thr
		phis <- phic[indc]

		if(all(pans>0)) return(xy0)
		
		# Delete positive points at the small r
		if(pans[1]>0){
			tmp <- rle(pans>0)
			pt <- tmp$lengths[1]
			pans[1:pt] <- pans[pt+1]
			pan[indc[1:pt]] <- pans[pt+1]
		}
		
		#Delete short spikes in data
		tmp <- rle(pans>0)

		ind_spike <- which(tmp$lengths<3)
		while(length(ind_spike)>0){

			t2 <- sum(tmp$lengths[1:ind_spike[1]])
			t1 <- t2 - tmp$lengths[ind_spike[1]] + 1

			if(t2==length(pans)) break
			
			if(t1>1)v1 <- pans[t1-1] else v1 <- pans[t2+1]
			if(t2 < length(pans))v2 <- pans[t2+1]
			
			pans[t1:t2] <- mean(v1,v2)
			pan[indc[t1:t2]] <- mean(v1,v2)
			tmp <- rle(pans>0)
			ind_spike <- which(tmp$lengths<4)
		}

		
		# Decide if region growing is required
		Grow <- FALSE
		if(all(pans<0))Grow <- TRUE
		
		# Correction
		if(pans[length(pans)]<0) Grow <- TRUE
		
		while(Grow){
			Smax <- Smax + mean(ps)
			indc <- pan_profile(phic,rc,Smax,r,ps,psi)
			rs <- rc[indc]
			pans <- pan[indc] - shad_thr
			phis <- phic[indc]
			l <- length(pans)

			if((all(pans[(l-2):l]>0))|(Smax>20*hmax)) Grow <- FALSE
		}
		
		if(Smax>=20*hmax) return(xy0)
		
		smoothingSpline <- smooth.spline(rs, pans, spar=0.2)

		pans <- smoothingSpline$y
		
		p1 <- min(pans)
		if(p1>0){
		# there is no shadow associated with this object
			return(xy0)
		}
		
		pt0 <- which(pans==p1)
		r0 <- rs[pt0]
		
		pt_in <- pt0 - 1 + rle(pans[pt0:length(pans)]<0)$lengths[1]
		pt_out <- pt_in + 1
		
		r1 <- rs[pt_in]
		r2 <- rs[pt_out]
		
		p1 <- pans[pt_in]
		p2 <- pans[pt_out]
		
		if(p1<0){
			rtip <- r1 + (r2-r1)*(0-p1)/(p2-p1)
		}else{
			rtip=NA
		}

		val <- xy0 + rtip*c(cos(psi),sin(psi))
		
	return(val)
}

show_shadow <- function(Trees,color){

	if(nrow(Trees)<1) return(-1)

	for(id in 1:nrow(Trees)){		
		xy0 <- c(Trees$x[id],Trees$y[id]) 
		r <- 0.5*Trees$d[id]
		h <- Trees$h[id]

		dr <- 0.5*h*tan(theta_sat)
		dx <- dr*sin(sat_az) 
		dy <- dr*cos(sat_az)

		xy1 <- xy0 + c(dx,dy)

		lines(c(xy0[1],xy1[1]),c(xy0[2],xy1[2]),col="blue")

		xshad <- c(xy1[1],xy1[1] - h*tan(theta_sun)*sin(sun_az))
		yshad <- c(xy1[2],xy1[2] - h*tan(theta_sun)*cos(sun_az))
		
		lines(xshad,yshad,col="red")

		rd <- rep(r,Nd)
		
		xd <- xy1[1] + rd*cos(phid)
		yd <- xy1[2] + rd*sin(phid)
		
		lines(xd,yd,col=color,lwd=1)
				
	}
	points(Trees$xtrue,Trees$ytrue,col="white",pch=16)

}

zero_cross<- function(y){
	tmp <- rle(y>0)
	n_cross <- length(tmp$lengths)-1
	
	x <- array(0,c(n_cross,2))
	#x[,1] <- tmp$lengths[1:n_cross]
	x[,1] <- cumsum(tmp$lengths[1:n_cross])
	x[,2] <- x[,1]+1
	
	val <- list(ncross=n_cross,cross=x)
	
	return(val)
}


# Compares two angles; first argument is a vector 
angle_diff<- function(phi,psi){

	val <- phi - psi
	val <- (val + pi) %% (2*pi) - pi

	return(val)
}

	

shadow_profile<-function(xy_center,s,xyt,alpha,ps){

	if(FALSE){
		xy_center<-xy1
		alpha <- alpha_shad
		xyt <- xy[shad,]
		xy <- xy.pan
		s <- S
		ps <- mean(ps.pan)
	}

	#dpsi <- 45*pi/180

	xy_rel <- xyt
	if(length(xyt)>2){
		xy_rel[,1] <- xy_rel[,1] - xy_center[1]
		xy_rel[,2] <- xy_rel[,2] - xy_center[2]

		r_rel <- sqrt(rowSums(xy_rel^2))
		phi_rel <- atan2(xy_rel[,2],xy_rel[,1])
	} else{
		xy_rel[1] <- xy_rel[1] - xy_center[1]
		xy_rel[2] <- xy_rel[2] - xy_center[2]

		r_rel <- sqrt(sum(xy_rel^2))
		phi_rel <- atan2(xy_rel[2],xy_rel[1])
	
	}

	s_max <- s

	dpsi_corr <- array(dpsi,length(r_rel))
	dpsi_corr <- dpsi_corr * ps/r_rel

	dphi <- angle_diff(phi_rel,alpha)
	
	ind_rel <- which((abs(dphi)<=dpsi_corr))
	if(length(ind_rel)==0) return(0)

	#points(xyt[ind_rel,],pch=16,col="yellow",cex=0.5)

	shad_length <- max(r_rel[ind_rel])

	# Refine by zero crossing estimate
	xy_rel <- xy
	xy_rel[,1] <- xy_rel[,1] - xy_center[1]
	xy_rel[,2] <- xy_rel[,2] - xy_center[2]

	r_rel <- sqrt(rowSums(xy_rel^2))
	phi_rel <- atan2(xy_rel[,2],xy_rel[,1])

	dpsi_corr <- array(dpsi,length(r_rel))
	dpsi_corr <- dpsi_corr * ps/r_rel

	dphi <- angle_diff(phi_rel,alpha)
	ind_rel <- which((abs(dphi)<=dpsi_corr)&(r_rel)<=shad_length+4*ps)

	#points(xy[ind_rel,],pch=16,col="white",cex=0.5)

	rs   <- r_rel[ind_rel]
	indr <- order(rs)
	ind_rel <- ind_rel[indr]
	rs <- rs[indr]
	pans <- pan[ind_rel] - shad_thr

	if(FALSE){
		windows()
		plot(rs,pans)
		lines(rs,pans)
		abline(v=s)
		abline(v=shad_length,lty=2)
		abline(h=0)
	}

	ind <- which(abs(rs-shad_length)==min(abs(rs-shad_length)))
	r0 <- mean(rs[ind])
	p0 <- pans[ind]

	if(p0>0){
		pans <- pans[1:ind]
		rs <- rs[1:ind]
		
	} else{
		pans <- pans[ind:length(pans)]
		rs <- rs[ind:length(rs)]
	}

	tmp <- zero_cross(pans)

	if(tmp$ncross==0){
		# Error
		# Do something
		return(0)
	} else{
		
		ind <- which.max(tmp$cross[,1])
		#ind <- which.min(tmp$cross[,1])
		tmp2 <- tmp$cross[ind,]

		r1<- rs[tmp2[1]]
		r2<- rs[tmp2[2]]
		p1<- pans[tmp2[1]]
		p2<- pans[tmp2[2]]
		
		rtip <- r1 + (r2-r1)*(0-p1)/(p2-p1)

		#abline(v=rtip,col="blue")
	}

	return(rtip)
}


f_pollock<-function(x,n){
	return((1-x^n)^(1/n))
}

Ey <- function(x,y,n){
	y_est <- f_pollock(x,n)
	return(mean((y_est-y)^2))
}

Ex <- function(x,y,n){
	x_est <- f_pollock(y,n)
	return(mean((x_est-x)^2))
}

E <- function(x,y,n){
	x_est <- f_pollock(y,n)
	y_est <- f_pollock(x,n)
	return(mean((x_est-x)^2+(y_est-y)^2))
}

dEdn<- function(x,y,n){
	ind <- order(x,decreasing=FALSE)
	x<-x[ind]
	y<-y[ind]
	
#	ind <- which((x>0)&(y>0))
	
	x_est <- f_pollock(y,n)
	y_est <- f_pollock(x,n)

	ind <- which((x>0)&(y>0)&(x_est>0)&(y_est>0))
	x<-x[ind]
	y<-y[ind]
	x_est<-x_est[ind]
	y_est<-y_est[ind]
	
	var <- (y_est-y)*(y_est*log(y_est)+((x/y_est)^n)*log(x))
	var <- var + (x_est-x)*(x_est*log(x_est)+((y/x_est)^n)*log(y))

	return(sum(var)*(-2)/n)
}

pollock_n <- function(x,y){
	n_arr <- seq(1,3.0,0.05)
	E_arr <- array(0,length(n_arr))

	for(i in 1:length(n_arr)){
		n <- n_arr[i]
		E_arr[i] <- E(x,y,n)
	}

	k <- which(E_arr==min(E_arr))
	n_est <- mean(n_arr[k])

	return(n_est)
}

pollock_n_v2 <- function(x,y){
	n_arr <- seq(1,3.0,0.05)
	E_arr <- array(0,length(n_arr))

	for(i in 1:length(n_arr)){
		n <- n_arr[i]
		E_arr[i] <- E(x,y,n)
	}

	k <- which(E_arr==min(E_arr))
	n_est <- mean(n_arr[k])

	x_est <- f_pollock(y,n)
	y_est <- f_pollock(x,n)
	ind <- which((x>0)&(y>0)&(x_est>0)&(y_est>0))
	return(c(n_est,min(E_arr),length(ind)))
}


pollock_gd_n <- function(x,y){
	# Initial value
	n<-2
	
	E_arr <- array(0,0)
	n_arr <- array(0,0)
	
	#windows()
	#plot(x,y,xlim=c(0,1),ylim=c(0,1))
	#xdisp <- seq(0,1,0.01)
	
	converged <-0
	
	while(TRUE){

		E1<-E(x,y,n)
		eps<-0.01
		dn <- -eps*dEdn(x,y,n)
		
		n<-n+dn
		E2<-E(x,y,n)
		E_arr <- c(E_arr,E2)
		n_arr <- c(n_arr,n)
		if((E1-E2<=0.0001)) converged <- converged + 1 else converged <- 0
		if(converged>=3) break
		
#		lines(xdisp,f_pollock(xdisp,n_est))
	}

#	windows()
#	par(mfrow=c(1,3))
#	plot(E_arr)
#	plot(n_arr)	
#	plot(x,y,xlim=c(0,1),ylim=c(0,1))
#	xdisp <- seq(0,1,0.01)
#	lines(xdisp,f_pollock(xdisp,n))
	
	return(n)
}

pollock_gd_n_v2 <- function(x,y){
	# Initial value
	n<-2
	
	E_arr <- array(0,0)
	n_arr <- array(0,0)
	
	#windows()
	#plot(x,y,xlim=c(0,1),ylim=c(0,1))
	#xdisp <- seq(0,1,0.01)
	
	converged <-0
	
	while(TRUE){

		E1<-E(x,y,n)
		eps<-0.01
		dn <- -eps*dEdn(x,y,n)
		
		n<-n+dn
		E2<-E(x,y,n)
		E_arr <- c(E_arr,E2)
		n_arr <- c(n_arr,n)
		if((E1-E2<=0.0001)) converged <- converged + 1 else converged <- 0
		if(converged>=3) break
		
#		lines(xdisp,f_pollock(xdisp,n_est))
	}

#	windows()
#	par(mfrow=c(1,3))
#	plot(E_arr)
#	plot(n_arr)	
#	plot(x,y,xlim=c(0,1),ylim=c(0,1))
#	xdisp <- seq(0,1,0.01)
#	lines(xdisp,f_pollock(xdisp,n))

	x_est <- f_pollock(y,n)
	y_est <- f_pollock(x,n)
	ind <- which((x>0)&(y>0)&(x_est>0)&(y_est>0))
	
	rmse <- E2
	
	return(c(n,rmse,length(ind)))
}


uv_tangent <- function(n,cx_arr,theta,af){
	val <- data.frame(array(0,c(length(cx_arr),2)))
	names(val) <- c("ut","vt")

	if(af==0){
		val$vt <- 1
		val$ut <- 0
		return(val)
	}
	
	der0 <- -1/(af*tan(theta))

	#cx_arr <- seq(0,1-0.01,0.01)
	vt_arr <- array(0,length(cx_arr))

	for(i in 1:length(cx_arr)){
		cx   <- cx_arr[i]
		vmax <- sqrt(1-cx^2)
		v    <- seq(0,vmax,0.001)
		u    <- (1-(cx^2+v^2)^(n/2))^(1/n)

		dudv <- -v*((1-u^n)^(1-2/n))*(u^(1-n))
		
		vt <- v[which.min(abs(dudv-der0))]
		if(length(vt)>0)vt_arr[i] <- vt
	}

	ut_arr <- (1-(cx_arr^2+vt_arr^2)^(n/2))^(1/n)
	val$vt <- vt_arr
	val$ut <- ut_arr
	
	return(val)
}

uv_tangent_polar <- function(n,phi_arr,theta,af){
	val <- data.frame(array(0,c(length(phi_arr),2)))

	der0 <- -1/(af*tan(theta))

	names(val) <- c("ut","vt")

	cx_arr <- cos(phi_arr)
	vt_arr <- array(0,length(phi_arr))

	for(i in 1:length(phi_arr)){
		cx   <- cos(phi_arr[i])
		vmax <- sqrt(1-cx^2)
		v    <- seq(0,vmax,0.001)
		u    <- (1-(cx^2+v^2)^(n/2))^(1/n)

		dudv <- -v*((1-u^n)^(1-2/n))*(u^(1-n))
		
		vt <- v[which.min(abs(dudv-der0))]
		if(length(vt)>0)vt_arr[i] <- vt
	}

	ut_arr <- (1-(cx_arr^2+vt_arr^2)^(n/2))^(1/n)
	val$vt <- vt_arr
	val$ut <- ut_arr
	
	return(val)
}


ndvi_shape<- function(Trees,A){
	var <- array(0,nrow(Trees))
	
	rmse <- array(0,nrow(Trees))
	
	xy <- coordinates(A)
	ps <- A@grid@cellsize
	
	for(id in 1:nrow(Trees)){
		val <- 0
		xy0 <- as.numeric(Trees[id,1:2])
		r <- as.numeric(0.5*Trees$d[id])

		xyc <- xy[,]
		xyc[,1] <- xyc[,1] - xy0[1]
		xyc[,2] <- xyc[,2] - xy0[2]
		rc <- sqrt(rowSums(xyc^2))
		phic <- atan2(xyc[,2],xyc[,1])

		ind <- which(rc<=r)

		x<-rc[ind]
		y<-A@data$ndvi[ind]

		x<-x/r
		y<-y-min(y)
		y<-y/max(y)

		#windows()
		#plot(x,y,xlim=c(0,1),ylim=c(0,1))

		#n_est <- pollock_n(x,y)
		n_est <- pollock_gd_n(x,y)
		var[id] <- n_est
		rmse[id] <- E(x,y,var[id])

		#xdisp <- seq(0,1,0.01)
		#lines(xdisp,f_pollock(xdisp,n_est))
		#title(round(n_est,digits=3))
	}

	return(var)
}

# Matching of two sets of blobs
# based on RANSAC algorithm
# assumes affine shift geometric transformation model
# r_search is the search radius for matching candidates
# req_precision is the distance threshold for matching acceptance
# n_samp is the size of random sample in RANSAC

blob_match <- function(Blobs_A,Blobs_B,r_search,req_precision,n_samp){
	
	if(any(is.null(Blobs_A),is.null(Blobs_B))){
		cat(as.character(Sys.time()),"pol",i.pol,"Error: no blobs to match","\n")
		return(-1)
	}
	
	if(nrow(Blobs_A)*nrow(Blobs_B)==0){
		cat(as.character(Sys.time()),"pol",i.pol,"Error: no blobs to match","\n")
		return(-1)
	}

	row.names(Blobs_A) <- 1:nrow(Blobs_A)
	row.names(Blobs_B) <- 1:nrow(Blobs_B)

	# Initial transformation parameters
	dx <- 0
	dy <- 0

	# Search radius, in metres
	r_search <- 40

	matched_pairs <- array(0,c(0,2))

	for(i in 1:nrow(Blobs_A)){

		blob <- Blobs_A[i,]

		tmp <- geo_transform(Blobs_B,dx,dy)

		tmp$r <- sqrt((tmp$x - blob$x)^2 + (tmp$y - blob$y)^2)

		ind <- order(tmp$r)
		tmp <- tmp[ind,]

		ind <- which(tmp$r<=r_search)
		if(length(ind)>0){ 
			tmp <- tmp[ind,]

			#tmp$d <- abs(tmp$d - blob$d)
			tmp$dratio <- pmin(tmp$d,blob$d)/pmax(tmp$d,blob$d)

			#tmp$ndvi <- abs(tmp$ndvi - blob$ndvi)
			#tmp$magn <- abs(tmp$magn-blob$magn)

			#order(tmp$d,decreasing=TRUE)
			
			ind <- which(tmp$dratio>=(1-0.1)*max(tmp$dratio))

			if(length(ind)>0)
			{
				ind_match <- as.numeric(row.names(tmp[ind,]))
				matched_pairs <- rbind(matched_pairs,cbind(i,ind_match))
			}
		}
	}

	# Optimizer of geometric model

	N_p <- nrow(matched_pairs)
	if(N_p<=3){
		cat(as.character(Sys.time()),"pol",i.pol,"Error: too few pair candidates found","\n")
		return(-1)
	}
		
	samples <- array(0,c(0,5))

	# random sample
	Nsamp <- min(N_p,n_samp)

	if(N_p<n_samp) sarr <- 1:N_p else sarr <- sample(1:N_p,Nsamp,replace=TRUE)

	for(iter in 1:Nsamp){
	#for(iter in 1:N_p){

		s <- sarr[iter]
		#s <- iter

		i <- matched_pairs[s,1]
		j <- matched_pairs[s,2]

		# Estimate model parameters
		xya <- Blobs_A[i,1:2]
		xyb <- Blobs_B[j,1:2]

		dxy <- as.numeric(xya - xyb)
		dx <- dxy[1]
		dy <- dxy[2]

		rest_pairs <- matched_pairs[-s,]
		rest_pairs <- cbind(rest_pairs,0)
		tmp <- geo_transform(Blobs_B,dx,dy)
			
		for(l in 1:nrow(rest_pairs)){

			i <- matched_pairs[l,1]
			j <- matched_pairs[l,2]

			xya <- Blobs_A[i,1:2]
			xyb <- tmp[j,1:2]

			res <- sqrt(sum((xya-xyb)^2))
			
			rest_pairs[l,3] <- res

		}

		ind <- which(rest_pairs[,3]<req_precision)

		n_cons <- length(ind)
		if(n_cons>0){
			res <- mean(rest_pairs[ind,3])

			samples <- rbind(samples, c(s, n_cons, dx,dy,res)) 
		}

	}

	# Check error handling
	if(max(samples[,2])<2){
		cat(as.character(Sys.time()),"pol",i.pol,"Error: matching failed","\n")
		next
	}

	ind0 <- which(samples[,2]==max(samples[,2]))
	if(length(ind0)>1) ind <- ind0[which.min(samples[ind0,5])] else ind <- ind0[1]

	s <- samples[ind,1]

	i <- matched_pairs[s,1]
	j <- matched_pairs[s,2]

	dx <- samples[ind,3]
	dy <- samples[ind,4]

	tmp <- geo_transform(Blobs_B,dx,dy)

	res_arr <- array(0,nrow(matched_pairs))	

	for(l in 1:nrow(matched_pairs)){

		i <- matched_pairs[l,1]
		j <- matched_pairs[l,2]

		xya <- Blobs_A[i,1:2]
		xyb <- tmp[j,1:2]

		res <- sqrt(sum((xya-xyb)^2))
		
		res_arr[l] <- res

	}

	# req_precision is the distance threshold
	ind <- which(res_arr<req_precision)

	n_cons <- length(ind)

	val <- matched_pairs[ind,]
	return(val)
}



blob_match_v2 <- function(Blobs_A,Blobs_B,r_search,req_precision,n_samp){
	
	if(any(is.null(Blobs_A),is.null(Blobs_B))){
		#cat(as.character(Sys.time()),"pol",i.pol,"Error: no blobs to match","\n")
		return(-1)
	}
	
	if((nrow(Blobs_A)==0) | (nrow(Blobs_B)==0)){
		#cat(as.character(Sys.time()),"pol",i.pol,"Error: no blobs to match","\n")
		return(-1)
	}

	row.names(Blobs_A) <- 1:nrow(Blobs_A)
	row.names(Blobs_B) <- 1:nrow(Blobs_B)

	# Initial transformation parameters
	dx <- 0
	dy <- 0


	A <- array(0,c(nrow(Blobs_A),3))
	for(k in 1:3) A[,k] <- Blobs_A[,k]

	B <- array(0,c(nrow(Blobs_B),3))
	for(k in 1:3) B[,k] <- Blobs_B[,k]


	# Find possible candidates based on the distance and diameter
	#z1 <- match_cand(Blobs_A,Blobs_B, 0,0,40)
	#matched_pairs <- match_cand(Blobs_A,Blobs_B, dx,dy,r_search)

	#Using C++ version
	matched_pairs <- match_candC(A, B, dx, dy, r_search)


	# Optimizer of geometric model

	N_p <- nrow(matched_pairs)
	if(N_p<=3){
		#cat(as.character(Sys.time()),"pol",i.pol,"Error: too few pair candidates found","\n")
		return(-1)
	}
		
	samples <- array(0,c(0,5))

	# random sample
	Nsamp <- min(N_p,n_samp)

	if(N_p<n_samp) sarr <- 1:N_p else sarr <- sample(1:N_p,Nsamp,replace=TRUE)

	#samples <- my_fun(A,B,sarr,matched_pairs,req_precision,Nsamp)
	samples <- optimizer_geomC(A,B,sarr,matched_pairs,req_precision,Nsamp)
	
	# Check error handling
	if(max(samples[,2])<2){
		#cat(as.character(Sys.time()),"pol",i.pol,"Error: matching failed","\n")
		return(-1)
	}

	ind0 <- which(samples[,2]==max(samples[,2]))
	if(length(ind0)>1) ind <- ind0[which.min(samples[ind0,5])] else ind <- ind0[1]

	samples[ind,]
	
	s <- samples[ind,1]

	i <- matched_pairs[s,1]
	j <- matched_pairs[s,2]

	dx <- samples[ind,3]
	dy <- samples[ind,4]

	val <- select_matchC(A, B, matched_pairs, dx, dy, req_precision)

	return(val)
}


transform_xy <- function(xy,xy0,alpha){
	
	if(length(xy)>2){
		x <- xy[,1]
		y <- xy[,2]
	} else{
		x <- xy[1]
		y <- xy[2]
	}
	
	x0<-xy0[1]
	y0 <- xy0[2]
	phi_arr <- atan2(y,x)
	r_arr   <- sqrt(x^2+y^2)
	
	phi_arr <- phi_arr + alpha
	
	xp <- x0 + r_arr*cos(phi_arr)
	yp <- y0 + r_arr*sin(phi_arr)
	
	if(length(xy)>2){
		return(cbind(xp,yp))
	} else return(c(xp,yp))
}

project_pollock_draw <-function(h,hf,R,x0,y0,theta,n,alpha){

	h1 <- h*hf*h0
	h2 <- h*h0 - h1

	af <- h2/R

	a  <- h2
	b  <- R

	xA <- 0
	yA <- h1*tan(theta) - R

	tmp <- uv_tangent(n,0,theta,af)
	z_tan <- as.numeric(tmp$ut)*h2
	r_tan <- as.numeric(tmp$vt)*R

	v    <- as.numeric(tmp$vt)
	u    <- as.numeric(tmp$ut)

	dudv <- -v*((1-u^n)^(1-2/n))*(u^(1-n))

	der <- dudv * af
	der0 <- -1/(tan(theta))
	if(abs(der-der0)>0.1){
		rp<-R
		zp<-0
	}

	xB <- 0
	yB <- r_tan
	zB <- h1+z_tan

	xD <- 0
	yD <- yB+zB*tan(theta)

	xC <- (xA+xD)/2
	yC <- (yA+yD)/2

	top_model <- 0.5*(r_tan-R+(z_tan+h1+h1)*tan(theta))
	Dt_model <- r_tan+R+z_tan*tan(theta)

	Dt_arr <- yD-yA
	top_arr <- yC-0

	# plot xy plane projection
	# plot the near side

	plot(0,0,xlim=x0+c(-1.5*R,1.5*R),ylim=y0+c(-2*R,2*R+(h1+h2)*tan(theta)),xlab="x",ylab="y",col=NA,asp=1)

	xarr <- seq(-R,R,0.25)
	yarr <- sqrt(R^2-xarr^2)

	xyarr<-transform_xy(cbind(xarr,yarr),c(x0,y0),alpha-pi/2-pi)
	lines(xyarr,lty=2)

	yarr <- -sqrt(R^2-xarr^2)

	xyarr<-transform_xy(cbind(xarr,yarr),c(x0,y0),alpha-pi/2-pi)

	lines(xyarr,lty=2)

	abline(h=0)
	abline(v=0)

	yarr <- h1*tan(theta) - sqrt(-xarr^2+R^2)

	xyarr<-transform_xy(cbind(xarr,yarr),c(x0,y0),alpha-pi/2-pi)
	lines(xyarr)

	xarr <- seq(0,R,0.25)
	yarr <- h1*tan(theta) - sqrt(-xarr^2+R^2)

	xarr <- seq(-R,R,0.25)

	tmp <- uv_tangent(n,xarr/R,theta,af)
	zp <- as.numeric(tmp$ut)*h2
	rp <- as.numeric(tmp$vt)*R

	yarr <- rp + (h1+zp)*tan(theta)

	xyarr<-transform_xy(cbind(xarr,yarr),c(x0,y0),alpha-pi/2-pi)
	lines(xyarr)

	xy <- cbind(rbind(0,xA,xB,xD),rbind(0,yA,yB,yD))

	xyp <- transform_xy(xy,c(x0,y0),alpha-pi/2-pi)

	points(xyp,pch=16)
	text(xyp,labels=c("O","A","B","D"),pos=c(4,3,3,3))

}

project_pollock <- function(h1,h2,R,x0,y0,theta,n,alpha){

	xy_shad <- array(0,c(0,2))

	af <- h2/R

	a  <- h2
	b  <- R

	xA <- 0
	yA <- h1*tan(theta) - R

	tmp <- uv_tangent(n,0,theta,af)
	z_tan <- as.numeric(tmp$ut)*h2
	r_tan <- as.numeric(tmp$vt)*R

	xB <- 0
	yB <- r_tan
	zB <- h1+z_tan

	xD <- 0
	yD <- yB+zB*tan(theta)

	xC <- (xA+xD)/2
	yC <- (yA+yD)/2

	top_model <- 0.5*(r_tan-R+(z_tan+h1+h1)*tan(theta))
	Dt_model <- r_tan+R+z_tan*tan(theta)

	Dt_arr <- yD-yA
	top_arr <- yC-0

	xarr <- seq(-R,R,0.25)
	#yarr <- sqrt(R^2-xarr^2)
	yarr <- -sqrt(R^2-xarr^2)

	xyarr <- transform_xy(cbind(xarr,yarr),c(x0,y0),alpha-pi/2-pi)
	
	yarr <- h1*tan(theta) - sqrt(-xarr^2+R^2)

	xyarr<-transform_xy(cbind(xarr,yarr),c(x0,y0),alpha-pi/2-pi)

	xy_shad <- rbind(xy_shad,xyarr)

	xarr <- seq(0,R,0.25)
	yarr <- h1*tan(theta) - sqrt(-xarr^2+R^2)

	xarr <- seq(-R,R,0.25)

	tmp <- uv_tangent(n,xarr/R,theta,af)
	zp <- as.numeric(tmp$ut)*h2
	rp <- as.numeric(tmp$vt)*R

	yarr <- rp + (h1+zp)*tan(theta)

	xyarr<-transform_xy(cbind(xarr,yarr),c(x0,y0),alpha-pi/2-pi)
	
	xy_shad <- rbind(xy_shad,xyarr[nrow(xyarr):1,])

	xy <- cbind(rbind(0,xA,xB,xD),rbind(0,yA,yB,yD))

	xy_shad <- rbind(xy_shad,xy_shad[1,])
	
	return(xy_shad)
}

project_pollock_quantitative <- function(h1,h2,R,x0,y0,theta,n){

	af <- h2/R

	a  <- h2
	b  <- R

	#xA <- 0
	#yA <- h1*tan(theta) - R

	tmp <- uv_tangent(n,0,theta,af)
	z_tan <- as.numeric(tmp$ut)*h2
	r_tan <- as.numeric(tmp$vt)*R

	#v    <- as.numeric(tmp$vt)
	#u    <- as.numeric(tmp$ut)

	#dudv <- -v*((1-u^n)^(1-2/n))*(u^(1-n))

	#der <- dudv * af
	#der0 <- -1/(tan(theta))
	#if(abs(der-der0)>0.1){
	#	rp<-R
	#	zp<-0
	#}

	#xB <- 0
	#yB <- r_tan
	#zB <- h1+z_tan

	#xD <- 0
	#yD <- yB+zB*tan(theta)

	#xC <- (xA+xD)/2
	#yC <- (yA+yD)/2

	top_model <- 0.5*(r_tan-R+(z_tan+h1+h1)*tan(theta))
	Dt_model <- r_tan+R+z_tan*tan(theta)

	#Dt_arr <- yD-yA
	#top_arr <- yC-0
	
	return(c(top_model,Dt_model))
}

eval_shad_v4 <- function(params,xy0,R,theta_sat,alpha_sat,theta_sun,alpha_sun){

	h  <- params[1]*h0
	hf <- params[2]
	n  <- params[3]

	h1 <- h*hf
	h2 <- h - h1

	tmp <- project_pollock_quantitative(h1,h2,R,xy0[1],xy0[2],theta_sat,n)
	xy1 <- xy0 + tmp[1]*c(cos(alpha_sat),sin(alpha_sat))

	shad_pol <- project_pollock(h1,h2,R,xy1[1],xy1[2],theta_sun,n,alpha_sun)

	# Get pixels inside the projected shadow region
	tmp.in <- point.in.polygon(xy[,1],xy[,2],shad_pol[,1],shad_pol[,2])
	ind_in <- shad[which(tmp.in[shad]==1)]
	
	ind_out <- shad[which(tmp.in[shad]==0)]
	
	ind <- which(tmp.in==1)
	
	val <- 2.0*length(ind_in) - length(ind_out) - 2.0*(length(ind)-length(ind_in))
	#val <- 2*length(ind_in) - length(ind_out) - length(ind)
	
	# normalize by area of the shadow mask
	val <- val/length(shad)
	
	return(val)
}	

interpolate_max <- function(x,y){

	n <- length(x)

	if(n<7){
		if(n==1) return(x)
		if(n==2) return(max(x))
		if(n==3){
			if(Debug) plot(x,y)
			return(parabol_interp_v2(x,y))
		}
		return(parabol_interp_v2(x,y))
	}

	# Fit polynomial 4th order
	pol_fit <- lm(y ~ x+I(x^2)+I(x^3)+I(x^4))
	coeff <- pol_fit$coefficients

	if(Debug){
		plot(x,y)
		xx <- seq(min(x),max(x),length.out=250)
		lines(xx,predict(pol_fit,data.frame(x=xx)),col="blue")
	}

	# get max value
	tmp <- polyroot(coeff[2:5]*(1:4))
	ind <- which((Re(tmp)>=min(x))&(Re(tmp)<=max(x))&(abs(Im(tmp))<1e-08))
	tmp <- Re(tmp[ind])
	
	if(length(tmp)==0){
		ind <- which(y==max(y))
		return(mean(x[ind]))
	}
	
	if(length(tmp)==1) val <- tmp else val <- tmp[which.max(predict(pol_fit,data.frame(x=tmp)))]

	return(val)
}	

grid_optim_v3 <- function(params,obj_fun,dpar=dpar,lower=lower,upper=upper){
	h  <- params[1]
	hf <- params[2]
	n  <- params[3]

	dh  <- dpar[1]
	dhf <- dpar[2]
	dn  <- dpar[3]
	
	h_min  <- lower[1]
	hf_min <- lower[2]
	n_min  <- lower[3]
	
	h_max  <- upper[1]
	hf_max <- upper[2]
	n_max  <- upper[3]
	
	h_arr <- seq(h_min,h_max,dh)
	obj_arr <- array(0,length(h_arr))

	# One at a time: h
	for(i1 in 1:length(h_arr)){
		h <- h_arr[i1]
		obj_arr[i1] <- obj_fun(c(h,hf,n),xy0,R,theta_sat,alpha_sat,theta_sun,alpha_sun)
	}

	if(Debug) windows()
	
	h <- determine_max_v2(h_arr,obj_arr)
	
	if(Debug){
		#windows()
		#plot(h_arr,obj_arr)
		title(main="h")
		abline(v=h,lty=2)
	}

	if(h==0){
		hf <-0
		n <- 2
		return(c(h,hf,n,obj_fun(c(h,hf,n),xy0,R,theta_sat,alpha_sat,theta_sun,alpha_sun)))
	}

	#Now hf
	hf_arr <- seq(hf_min,hf_max,dhf)
	obj_arr <- array(0,length(hf_arr))

	# One at a time: hf
	for(i1 in 1:length(hf_arr)){
		hf <- hf_arr[i1]
		obj_arr[i1] <- obj_fun(c(h,hf,n),xy0,R,theta_sat,alpha_sat,theta_sun,alpha_sun)
	}

	if(Debug) windows()
	
	hf <- determine_max_v2(hf_arr,obj_arr)

	if(Debug){
		#windows()
		#plot(hf_arr,obj_arr)
		title(main="hf")
		abline(v=hf,lty=2)
	}

	#Now n
	n_arr <- seq(n_min,n_max,dn)
	obj_arr <- array(0,length(n_arr))

	for(i1 in 1:length(n_arr)){
		n <- n_arr[i1]
		obj_arr[i1] <- obj_fun(c(h,hf,n),xy0,R,theta_sat,alpha_sat,theta_sun,alpha_sun)
	}
	
	if(Debug) windows()

	n <- determine_max_v2(n_arr,obj_arr)

	if(Debug){
		#windows()
		#plot(n_arr,obj_arr)
		title(main="n")
		abline(v=n,lty=2)	
	}

	return(c(h,hf,n,obj_fun(c(h,hf,n),xy0,R,theta_sat,alpha_sat,theta_sun,alpha_sun)))
}

grid_optim_v4 <- function(params,obj_fun,dpar=dpar,lower=lower,upper=upper){
	h  <- params[1]
	hf <- params[2]
	n  <- params[3]

	dh  <- dpar[1]
	dhf <- dpar[2]
	dn  <- dpar[3]
	
	h_min  <- lower[1]
	hf_min <- lower[2]
	n_min  <- lower[3]
	
	h_max  <- upper[1]
	hf_max <- upper[2]
	n_max  <- upper[3]
	
	h_arr <- seq(h_min,h_max,dh)
	obj_arr <- array(0,length(h_arr))

	# One at a time: h
	for(i1 in 1:length(h_arr)){
		h <- h_arr[i1]
		obj_arr[i1] <- obj_fun(c(h,hf,n),xy0,R,theta_sat,alpha_sat,theta_sun,alpha_sun)
	}

	if(Debug) windows()

	if(min(obj_arr)<max(obj_arr)){
		obj_arr <- obj_arr-min(obj_arr,na.rm=TRUE)
		obj_arr <- obj_arr/max(obj_arr,na.rm=TRUE)
		
		h <- determine_max_v2(h_arr,obj_arr)

		ind <- which(obj_arr >= 0.9*max(obj_arr))
		h_min <- min(h_arr[ind])-dh/2
		if(h_min<0) h_min <- 0
		h_max <- max(h_arr[ind])+dh/2
		if(h_max>h0) h_max <- h0
	}
	
	if(Debug){
		#windows()
		plot(h_arr,obj_arr)
		title(main="h")
		abline(v=h,lty=2)
		abline(v=h_min,lty=2,col="blue")
		abline(v=h_max,lty=2,col="red")
	}

	if(h==0){
		hf <-0
		n <- 2
		return(c(h,hf,n,obj_fun(c(h,hf,n),xy0,R,theta_sat,alpha_sat,theta_sun,alpha_sun),h_min,h_max,hf_min,hf_max,n_min,n_max))
	}

	#Now hf
	hf_arr <- seq(hf_min,hf_max,dhf)
	obj_arr <- array(0,length(hf_arr))

	# One at a time: hf
	for(i1 in 1:length(hf_arr)){
		hf <- hf_arr[i1]
		obj_arr[i1] <- obj_fun(c(h,hf,n),xy0,R,theta_sat,alpha_sat,theta_sun,alpha_sun)
	}

	if(Debug) windows()

	if(min(obj_arr)<max(obj_arr)){
		obj_arr <- obj_arr-min(obj_arr,na.rm=TRUE)
		obj_arr <- obj_arr/max(obj_arr,na.rm=TRUE)
		
		hf <- determine_max_v2(hf_arr,obj_arr)

		ind <- which(obj_arr >= 0.9*max(obj_arr))
		hf_min <- min(hf_arr[ind])-dhf/2
		if(hf_min<0) hf_min <- 0
		hf_max <- max(hf_arr[ind])+dhf/2
		if(hf_max>0.9) hf_max <- 0.9
	}

	if(Debug){
		#windows()
		plot(hf_arr,obj_arr)
		title(main="hf")
		abline(v=hf,lty=2)
		abline(v=hf_min,lty=2,col="blue")
		abline(v=hf_max,lty=2,col="red")
	}

	#Now n
	n_arr <- seq(n_min,n_max,dn)
	obj_arr <- array(0,length(n_arr))

	for(i1 in 1:length(n_arr)){
		n <- n_arr[i1]
		obj_arr[i1] <- obj_fun(c(h,hf,n),xy0,R,theta_sat,alpha_sat,theta_sun,alpha_sun)
	}

	if(Debug) windows()

	if(min(obj_arr)<max(obj_arr)){
		obj_arr <- obj_arr-min(obj_arr,na.rm=TRUE)
		obj_arr <- obj_arr/max(obj_arr,na.rm=TRUE)

		n <- determine_max_v2(n_arr,obj_arr)

		ind <- which(obj_arr >= 0.9*max(obj_arr))
		n_min <- min(n_arr[ind])-dn/2
		if(n_min<1.1) n_min <- 1.1
		n_max <- max(n_arr[ind])+dn/2
		if(n_max>3.0) n_min <- 3.0
	}

	if(Debug){
		#windows()
		plot(n_arr,obj_arr)
		title(main="n")
		abline(v=n,lty=2)
		abline(v=n_min,lty=2,col="blue")
		abline(v=n_max,lty=2,col="red")
	}

	return(c(h,hf,n,obj_fun(c(h,hf,n),xy0,R,theta_sat,alpha_sat,theta_sun,alpha_sun),h_min,h_max,hf_min,hf_max,n_min,n_max))
}


draw_shadow <- function(xy0,R,h,hf,n,theta_sat,alpha_sat,theta_sun,alpha_sun){

	h1 <- h*hf
	h2 <- h-h1

	tmp <- project_pollock_quantitative(h1,h2,R,xy0[1],xy0[2],theta_sat,n)
	xy1 <- xy0 + tmp[1]*c(cos(alpha_sat),sin(alpha_sat))

	rd <- rep(R,Nd)

	xd <- xy1[1] + rd*cos(phid)
	yd <- xy1[2] + rd*sin(phid)

	points(xy1[1],xy1[2],col="blue",pch=16)
	lines(xd,yd,col="blue",lwd=1)

	shad_pol <- project_pollock(h1,h2,R,xy1[1],xy1[2],theta_sun,n,alpha_sun)
	lines(shad_pol,col="white")
}

draw_shad_fit <- function(title_str){
	#windows()
	image(Pan_sub,col=gray((0:255)/255),axes=TRUE)
	display_all_blobs(Blobs[id,],"green")
	title(main=title_str)

	# Draw the shadow outline for h,hf,n
	draw_shadow(xy0,R,h*h0,hf,n,theta_sat,alpha_sat,theta_sun,alpha_sun)
}

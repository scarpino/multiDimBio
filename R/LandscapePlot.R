LandscapePlot <- function(Data, Groups = NULL, PDF = FALSE, LocPlot = FALSE, 
    control = c(75, 1, 30)) {
    
    if (is.matrix(Data) == FALSE) {
        cat("Data must be a matrix", "\n", "\n")
    }
    
    if (ncol(Data) > 9) {
        cat("Does not plot more than 9 traits", "\n", "\n")
    }
    
    
    mu1 <- 0  # setting the expected value of x1
    mu2 <- 0  # setting the expected value of x2
    s11 <- 10  # setting the variance of x1
    s12 <- 15  # setting the covariance between x1 and x2
    s22 <- 10  # setting the variance of x2
    rho <- 0.5  # setting the correlation coefficient between x1 and x2
    x1 <- seq(-10, 10, length = 61)  # generating the vector series x1
    x2 <- x1  # copying x1 to x2
    
    f <- function(x1, x2) {
        term1 <- 1/(2 * pi * sqrt(s11 * s22 * (1 - rho^2)))
        term2 <- -1/(2 * (1 - rho^2))
        term3 <- (x1 - mu1)^2/s11
        term4 <- (x2 - mu2)^2/s22
        term5 <- -2 * rho * ((x1 - mu1) * (x2 - mu2))/(sqrt(s11) * 
            sqrt(s22))
        term1 * exp(term2 * (term3 + term4 - term5))
    }  # setting up the function of the multivariate normal density
    M <- outer(x1, x2, f)  #
    
    GRS <- unique(Groups)
    if (length(GRS) == 0) {
        GRS <- "Group"
        Groups <- rep(GRS, length = nrow(Data))
    }
    
    for (j in 1:length(GRS)) {
        use.j <- which(Groups == GRS[j])
        DAT.j <- Data[use.j, ]
        DAT.j <- as.matrix(DAT.j)
        if (ncol(DAT.j) > 1) {
            mu.j <- colMeans(DAT.j, na.rm = TRUE)
        } else {
            mu.j <- mean(DAT.j, na.rm = TRUE)
        }
        
        if (length(mu.j) == 1) {
            M.plot <- M * mu.j
            posx <- c(0)
            posy <- c(0.5)
            labels <- c("1")
        }
        if (length(mu.j) == 2) {
            MA <- cbind(M * 0, M * mu.j[2])
            MB <- cbind(M * mu.j[2], M * 0)
            M.plot <- rbind(MA, MB)
            posx <- c(0.75, -0.75)
            posy <- c(1, 0.5)
            labels <- c("1", "2")
        }
        if (length(mu.j) == 3) {
            M.plot <- matrix(0, nrow = 61 * 2, ncol = 61 * 3)
            M.plot[1:61, 62:122] <- M * mu.j[1]
            M.plot[62:122, 1:61] <- M * mu.j[2]
            M.plot[62:122, 123:183] <- M * mu.j[3]
            M.fill <- matrix(, ncol = 183, nrow = 61)
            M.plot <- rbind(M.fill, M.plot)
            posx <- c(0.5, -0.75, 0.5)
            posy <- c(0.5, 0, -0.15)
            labels <- c("1", "2", "3")
        }
        if (length(mu.j) == 4) {
            MA <- cbind(M * mu.j[1], M * mu.j[2])
            MB <- cbind(M * mu.j[3], M * mu.j[4])
            M.plot <- rbind(MA, MB)
            posx <- c(-0.5, 0.75, -1, 0.2)
            posy <- c(1, 0.75, 0, -0.5)
            labels <- c("1", "2", "3", "4")
        }
        if (length(mu.j) == 5) {
            M.plot <- matrix(0, nrow = 61 * 3, ncol = 61 * 3)
            M.plot[1:61, 62:122] <- M * mu.j[1]
            M.plot[62:122, 1:61] <- M * mu.j[2]
            M.plot[62:122, 123:183] <- M * mu.j[3]
            M.plot[123:183, 1:61] <- M * mu.j[4]
            M.plot[123:183, 123:183] <- M * mu.j[5]
            posx <- c(0.5, -0.3, 0.3, -0.75, 0)
            posy <- c(1, 0.4, 0.4, 0, 0)
            labels <- c("1", "2", "3", "4", "5")
        }
        if (length(mu.j) == 6) {
            M.plot <- matrix(0, nrow = 61 * 4, ncol = 61 * 3)
            M.plot[1:61, 62:122] <- M * mu.j[1]
            M.plot[62:122, 1:61] <- M * mu.j[2]
            M.plot[62:122, 123:183] <- M * mu.j[3]
            M.plot[123:183, 1:61] <- M * mu.j[4]
            M.plot[123:183, 123:183] <- M * mu.j[5]
            M.plot[184:244, 62:122] <- M * mu.j[6]
            M.fill <- matrix(, ncol = 61, nrow = 244)
            M.plot <- cbind(M.fill, M.plot)
            posx <- c(0.5, -0.3, 0.3, -0.75, 0, -0.75)
            posy <- c(1, 0.6, 0.4, 0.25, 0.25, 0)
            labels <- c("1", "2", "3", "4", "5", "6")
        }
        if (length(mu.j) == 7) {
            M.plot <- matrix(0, nrow = 61 * 4, ncol = 61 * 5)
            M.plot[1:61, 123:183] <- M * mu.j[1]
            M.plot[62:122, 62:122] <- M * mu.j[2]
            M.plot[62:122, 184:244] <- M * mu.j[3]
            M.plot[123:183, 1:61] <- M * mu.j[4]
            M.plot[123:183, 245:305] <- M * mu.j[5]
            M.plot[184:244, 62:122] <- M * mu.j[6]
            M.plot[184:244, 184:244] <- M * mu.j[7]
            M.fill <- matrix(, ncol = 305, nrow = 61)
            M.plot <- rbind(M.fill, M.plot)
            posx <- c(0.5, -0.3, 0.5, -0.75, 0.75, -0.75, 0.2)
            posy <- c(0.6, 0.4, 0.25, 0.25, -0.2, -0.1, -0.65)
            labels <- c("1", "2", "3", "4", "5", "6", "7")
        }
        if (length(mu.j) == 8) {
            M.plot <- matrix(0, nrow = 61 * 4, ncol = 61 * 4)
            M.plot[1:61, 62:122] <- M * mu.j[1]
            M.plot[1:61, 123:183] <- M * mu.j[2]
            M.plot[62:122, 1:61] <- M * mu.j[3]
            M.plot[62:122, 184:244] <- M * mu.j[4]
            M.plot[123:183, 1:61] <- M * mu.j[5]
            M.plot[123:183, 184:244] <- M * mu.j[6]
            M.plot[184:244, 62:122] <- M * mu.j[7]
            M.plot[184:244, 123:183] <- M * mu.j[8]
            posx <- c(0.2, 0.6, -0.3, 0.7, -0.75, 0.75, -0.75, 
                -0.1)
            posy <- c(0.5, 0.5, 0.4, 0, 0.25, -0.35, -0.1, -0.2)
            labels <- c("1", "2", "3", "4", "5", "6", "7", "8")
        }
        if (length(mu.j) == 9) {
            M.plot <- matrix(0, nrow = 61 * 5, ncol = 61 * 6)
            M.plot[1:61, 184:244] <- M * mu.j[1]
            M.plot[62:122, 62:122] <- M * mu.j[2]
            M.plot[62:122, 245:305] <- M * mu.j[3]
            M.plot[123:183, 1:61] <- M * mu.j[4]
            M.plot[123:183, 306:366] <- M * mu.j[5]
            M.plot[184:244, 62:122] <- M * mu.j[6]
            M.plot[184:244, 245:305] <- M * mu.j[7]
            M.plot[245:305, 123:183] <- M * mu.j[8]
            M.plot[245:305, 184:244] <- M * mu.j[9]
            M.fill <- matrix(, ncol = 366, nrow = 61)
            M.plot <- rbind(M.fill, M.plot)
            posx <- c(0.6, -0.3, 0.6, -0.75, 0.75, -0.75, 0.4, 
                -0.75, -0.1)
            posy <- c(0.5, 0.4, 0.1, 0.25, -0.35, -0.1, -0.45, 
                -0.2, -0.3)
            labels <- c("1", "2", "3", "4", "5", "6", "7", "8", 
                "9")
        }
        if (length(mu.j) > 9) {
            cat("Warning! Only the first 9 axes are shown", "\n", 
                "\n")
            M.plot <- matrix(0, nrow = 61 * 5, ncol = 61 * 6)
            M.plot[1:61, 184:244] <- M * mu.j[1]
            M.plot[62:122, 62:122] <- M * mu.j[2]
            M.plot[62:122, 245:305] <- M * mu.j[3]
            M.plot[123:183, 1:61] <- M * mu.j[4]
            M.plot[123:183, 306:366] <- M * mu.j[5]
            M.plot[184:244, 62:122] <- M * mu.j[6]
            M.plot[184:244, 245:305] <- M * mu.j[7]
            M.plot[245:305, 123:183] <- M * mu.j[8]
            M.plot[245:305, 184:244] <- M * mu.j[9]
            M.fill <- matrix(, ncol = 366, nrow = 61)
            M.plot <- rbind(M.fill, M.plot)
            posx <- c(0.6, -0.3, 0.6, -0.75, 0.75, -0.75, 0.4, 
                -0.75, -0.1)
            posy <- c(0.5, 0.4, 0.1, 0.25, -0.35, -0.1, -0.45, 
                -0.2, -0.3)
            labels <- c("1", "2", "3", "4", "5", "6", "7", "8", 
                "9")
        }
        
        timestamp <- as.character(as.integer(Sys.time()))
        print(timestamp)
        if (PDF == TRUE) {
            pdf(paste(Groups[j], timestamp, "Functional Landscape.pdf"))
            persp(seq(0, 1, length.out = nrow(M.plot)), seq(0, 
                1, length.out = nrow(M.plot)), M.plot, box = FALSE, 
                col = "#F0F0F0", border = "#969696", theta = control[1], 
                r = control[2], phi = control[3])
            dev.off()
        } else {
            jpeg(paste(Groups[j], timestamp, "Functional Landscape.jpeg"), 
                width = 3840, height = 3840, quality = 100)
            persp(seq(0, 1, length.out = nrow(M.plot)), seq(0, 
                1, length.out = nrow(M.plot)), M.plot, box = FALSE, 
                col = "#F0F0F0", border = "#969696", theta = control[1], 
                r = control[2], phi = control[3])
            dev.off()
        }  #end if/else if PDF==TRUE
        if (LocPlot == TRUE) {
            pdf(paste(Groups[j], timestamp, "Functional Landscape-NAMES.pdf"))
            drawScene(surfaceTriangles(1:nrow(M.plot), 1:ncol(M.plot), 
                M.plot, color = "gray"), screen = list(z = 40, 
                x = -60))
            text(posx, posy, labels)
            dev.off()
        }
    }  #end for j
    
}  #end FUNCTION

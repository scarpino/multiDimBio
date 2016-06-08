#' A Function to perform step-wise discriminant analysis using F statistics
#' 
#' Select data using a F tests
#' 
#' 
#' @param Data A (non-empty), numeric matrix of data values
#' @param Group A (non-empty), vector indicating group membership.
#' Length(unique(Group))==2
#' @param target The number of desired traits. Target cannot be greater than
#' the number of columns in Data
#' @param p.adj.method The method used to control for false discovery.  The
#' default setting is 'holm'
#' @param Missing.Data The method used to handle missing data.  The default,
#' 'Complete' will use completeData to impute missing data, setting
#' Missing.Data='Remove' will remove all individuals with missing data.
#' FSelect cannot handle missing data.
#' @return FSelect returns list containing at least the following components:
#' 
#' \item{Selected}{ An ordered list indicating which columns were selected.  }
#' \item{F.Selected}{ An ordered list containing the F statistics for each
#' column indicated in Selected.  } \item{PrF}{ An ordered list containing the
#' p values for each column indicated in Selected. } \item{PrNotes}{ A string
#' indicating which method was used to control for multiple comparisons }
#' \item{model}{ An lm object with the final model results.  }
#' @seealso \code{\link{completeData}}
#' @references Costanza M, Afifi A (1979). Comparison of Stopping Rules in
#' Forward Stepwise Discriminant Analysis. Journal of the American Statistical
#' Association, pp. 777 - 78
#' 
#' Habbema J, Hermans J (1977). Selection of Variables in Discriminant Analysis
#' by F - Statistics and Error Rate. Technometrics, 19(4), 487 - 493.
#' 
#' Jennrich R (1977). Stepwise discriminant analysis, volume 3. New York Wiley
#' Sons.
#' @examples
#' 
#' data(Nuclei)
#' data(Groups)
#' npcs<-floor(ncol(Nuclei)/5)
#' 
#' dat.comp <- completeData(data = Nuclei, n_pcs = npcs) 
#' groups.use <- c(1,2) 
#' use.dat <- which(Groups==groups.use[1]|Groups==groups.use[2])
#' 
#' dat.use <- Nuclei[use.dat,]
#' GR.use <- Groups[use.dat]
#' 
#' #not run
#' #FSelect(DAT.use,GR.use,3)
#' 
#' @export FSelect
FSelect <- function(Data, Group, target, p.adj.method = "holm", 
    Missing.Data = "Complete") {
    
    Group <- as.factor(Group)
    Group <- as.numeric(Group)
    
    # Missing data
    if (Missing.Data == "Complete" & length(which(is.na(Data) == 
        TRUE)) > 0) {
        cat("Missing data imputed using completeData to exlude missing data set Missing.Data=Remove", 
            "\n", "\n")
        Data <- completeData(Data, n_pcs = floor(ncol(Data)/5), 
            cut.trait = 1, cut.ind = 1)$complete_dat
    }
    
    if (Missing.Data == "Remove" & length(which(is.na(Data) == 
        TRUE)) > 0) {
        cat("Individuals with missing data removed, to impute missing data set Missing.Data=Complete", 
            "\n", "\n")
        rm <- na.omit(Data)
        Data <- Data[-attr(rm, "na.action"), ]
        Group <- Group[-attr(rm, "na.action")]
    }
    
    F.selected <- c()
    selected <- c()
    T.selected <- c()
    PrFs <- c()
    GRs <- unique(Group)
    if (length(GRs) > 2) {
        print("This isnt going to work, n.Groups needs to equal 2")
    }
    n1 <- length(which(Group == GRs[1]))
    n2 <- length(which(Group == GRs[2]))
    df1 <- n1 - 1
    df2 <- n2 - 1
    counter = 0
    if (min(Data) < 1e-04) {
        TOL <- min(Data)
    } else {
        TOL <- 1e-04
    }
    
    while (counter < target) {
        
        cols <- 1:ncol(Data)
        
        if (counter == 0) {
            Fs <- c()
            Ts <- c()
            for (c in cols) {
                m.i <- lda(Group ~ Data[, c], tol = TOL)
                scores <- predict(m.i)
                F.i <- partialF(m.i, Group, 0)
                # F.i<-V.w(scores$x[,1],Group)/V.b(scores$x[,1],Group)
                T.i <- sum(m.i$svd)
                
                Fs <- c(Fs, F.i)
                Ts <- c(Ts, T.i)
            }  #end for c in 1:ncol(Data)
            select <- which(Fs == max(Fs))
            PrF <- df(max(Fs), df1, df2)
            PrF <- p.adjust(PrF, method = p.adj.method, c)
            selected <- c(selected, cols[select])
            F.selected <- c(F.selected, max(Fs))
            T.selected <- c(T.selected, Ts[select])
            PrFs <- c(PrFs, PrF)
        } else {
            c.in <- cols[-selected]
            
            Fs <- c()
            Ts <- c()
            for (c in c.in) {
                m.i <- lda(Group ~ Data[, c(selected, c)], tol = TOL)
                scores <- predict(m.i)
                F.i <- partialF(m.i, Group, T.selected[counter])
                # F.i<-V.w(scores$x[,1],Group)/V.b(scores$x[,1],Group)
                T.i <- sum(m.i$svd)
                
                Fs <- c(Fs, F.i)
                Ts <- c(Ts, T.i)
            }  #end for c in c.in
            select <- which(Fs == max(Fs))
            if (max(Fs) < 1) {
                PrF <- 1
            } else {
                PrF <- df(max(Fs), df1, df2)
                PrF <- p.adjust(PrF, method = p.adj.method, c)
            }  #end if/else max(Fs)
            selected <- c(selected, c.in[select])
            F.selected <- c(F.selected, max(Fs))
            T.selected <- c(T.selected, Ts[select])
            PrFs <- c(PrFs, PrF)
            
        }  #end if/else counter
        
        counter <- counter + 1
        
    }  #end while loop
    # running the final model
    m.final <- lda(Group ~ Data[, selected], tol = TOL)
    
    notes <- paste("PrF has been", p.adj.method, "adjusted for", 
        counter, "comparisons", sep = " ")
    
    results <- list("Selected" = selected, "F.Selected" = F.selected, 
        "PrF" = PrFs, "PrNotes" = notes, "model" = m.final)
    
    return(results)
}  #end F.select


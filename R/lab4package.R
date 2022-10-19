library(ggplot2)
library(gridExtra)
library(methods)


#' linreg
#'
#' linreg() is a RC class used for create a linear regression
#'
#' @param formula a formula
#' @param data a data.frame
#' @return attributes and methods of a new RC class
#' @import methods
#' @import ggplot2
#' @import gridExtra
#' @export linreg
#' @export




linreg <- setRefClass("linreg",
                      fields = list(
                        formula = "formula",
                        data = "data.frame",
                        dataname = "character",
                        x = "matrix",
                        y = "matrix",
                        beta_hat ="matrix",
                        y_hat = "matrix",
                        e_hat = "matrix",
                        df_1 = "numeric",
                        MSE = "matrix",
                        The_variance_of_beta = "numeric",
                        t_value_for_beta = "matrix",
                        p_value = "matrix",
                        sum = "data.frame"),
                      methods = list(
                        initialize =  function(formula,data){
                          formula <<- formula
                          data <<- data
                          dataname <<- deparse(substitute(data))
                          x <<- model.matrix(formula,data = data)
                          y <<- as.matrix(data[,all.vars(formula)][1])
                          beta_hat <<- solve(t(x) %*% x) %*% t(x) %*% y
                          y_hat <<- x %*% beta_hat
                          e_hat <<- y - (x %*% beta_hat)
                          df_1 <<- nrow(x) - nrow(beta_hat)
                          MSE <<- (t(e_hat) %*% e_hat) / df_1
                          The_variance_of_beta <<- MSE[1] * diag(solve(t(x) %*% x))
                          t_value_for_beta <<- beta_hat / sqrt(The_variance_of_beta)
                          p_value <<- pt(t_value_for_beta, df_1)
                          sum <<- data.frame(beta_hat,sqrt(The_variance_of_beta),t_value_for_beta,p_value)
                          colnames(sum) <<- c("coefficients","standard error","t_value","P_value")
                        },

                        print =  function(){
                          cat("linreg(formula = ", deparse(formula), ", data = ",dataname,")", sep = "")
                          cat("\n\nCoefficients:\n",sep = "")
                          print.default(format(coef()))

                        },

                        resid = function(){
                          return(as.vector(e_hat))
                        },

                        pred = function(){
                          return(y_hat)
                        },

                        coef = function(){
                          coef_vector <- c()
                          names_coef <- c()
                          len <- length(beta_hat)
                          for(i in 1:len){
                            coef_vector[i] <- beta_hat[i,1]
                          }
                          names_coef <- rownames(beta_hat)
                          names(coef_vector) <- names_coef
                          return(coef_vector)
                        },

                        summary = function(){
                          sum[,4] <<- "***"
                          print.data.frame(sum,digits = 3)
                          cat("\nResidual standard error:",sqrt(MSE),"on",df_1,"degrees of freedom",sep = " ")
                        },

                        plot = function(){
                          data_1 <- data.frame(y_hat,e_hat)
                          data_2 <- data.frame(y_hat,sqrt(abs(e_hat / sqrt(MSE[1,1]))))

                          p1 <- ggplot(data = data_1,aes(x = data_1[,1],y = data_1[,2])) + geom_point(size=3,shape=1) +
                            geom_smooth(formula = y ~ x,method="lm",color="red",se=FALSE) +
                            xlab(paste("Fitted values\n","lm(",deparse(formula),")")) + ylab("Residuals") +
                            ggtitle("Residuals VS Fitted")

                          p2 <- ggplot(data = data_2,aes(x = data_2[,1],y = data_2[,2])) + geom_point(size=3,shape=1) +
                            geom_smooth(formula = y ~ x,method="lm",color="red",se=FALSE) +
                            xlab(paste("Fitted valuse\n","lm(",deparse(formula),")")) + ylab(expression(sqrt("|Standardized residules|"))) +
                            ggtitle("Scale-Location")

                          grid.arrange(p1,p2)



                        }


                      )
)

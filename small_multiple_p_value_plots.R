## Remove all previous variables
rm(list=ls(all=TRUE))
library(ggplot2)


draw.grid <- function(input.params)
{
  lower.sweep.value <- input.params[1]
  upper.sweep.value <- input.params[2]
  step.size <- input.params[3]
  
  ################################################################################
  ## Define the intervals in the parameter space that we are going to sweep
  ################################################################################
  pbar.opt <- pbar.ctrl <- seq(lower.sweep.value, upper.sweep.value, step.size) ## Proportions of optimized and control groups
  
  powers <- 0:7 ## We will raise 10 to this vector to get a log-sweep list
  n.opt <- n.ctrl <- 10^powers ## Sample sizes of optimized and control groups

  ## Use these to test against the book's solved example:
  # pbar.opt <- .69
  # pbar.ctrl <- .76
  # n.opt <- 300
  # n.ctrl <- 350
  
  ## Use these to test against the stattrek example, which has a right answer of p=.034
  # pbar.opt <- .38
  # pbar.ctrl <- .51
  # n.opt <- 100
  # n.ctrl <- 200
  
  ## Build a df with a single row for each combination of our params, covering every permutation
  parameter.space <- expand.grid(n.opt = n.opt, n.ctrl = n.ctrl, pbar.opt = pbar.opt, pbar.ctrl = pbar.ctrl)
  
  ################################################################################
  ################################################################################
  
  
  
  ################################################################################
  ## Check to see if our proportion can be approximated by ordinary Z stats
  ## n*p must be >= 5 and also nq >= 5 where q = 1 - p
  ################################################################################
  
  # # First, check for the optimized proportion
  # parameter.space <- transform(parameter.space, np.opt = n.opt * pbar.opt)
  # parameter.space<- transform(parameter.space, nq.opt = n.opt * (1 - pbar.opt))
  # 
  # ## Second, check for the control proportion
  # parameter.space <- transform(parameter.space, np.ctrl = n.ctrl * pbar.ctrl)
  # parameter.space<- transform(parameter.space, nq.ctrl = n.ctrl * (1 - pbar.ctrl))
  # 
  # # Now remove all rows where np < 5 or nq < 5
  # enough.events <- 5
  # parameter.space <- subset(parameter.space, np.opt >= enough.events & np.ctrl >= enough.events & nq.opt >= enough.events & nq.ctrl >= enough.events)
  # parameter.space <- parameter.space[!(parameter.space$np.opt <= 5),]
  ################################################################################
  ################################################################################
  
  
  
  ################################################################################
  ## Build a function that will calculate the p values (and several intermediate
  ## results, also) and then return a df containing them
  ################################################################################
  
  calculate.p.values <- function(parameter.space)
    {
    n1 = parameter.space['n.opt']
    n2 = parameter.space['n.ctrl']
    pbar1 = parameter.space['pbar.opt']
    pbar2 = parameter.space['pbar.ctrl']
    
    ## Calculate pooled proportion
    p.hat <- (pbar1*n1 + pbar2*n2) / (n1 + n2)
    colnames(p.hat)[1] <- "p.hat" ## Fix the column name
    
    ## Get estimated standard error of difference of two proportions
    sigma.hat <- sqrt( p.hat * (1 - p.hat) * (1 / n1 + 1 / n2) ) ## Approximate, use for sample p's aka pbar
    colnames(sigma.hat)[1] <- "sigma.hat" ## Fix the column name
    
    ## Calculate the test statistic, which in this case is a Z-stat
    calculated.z.score <- ( (pbar1 - pbar2) - 0)/( sigma.hat ) ## Set second term in numerator = 0 because null hypothesis is that there is no difference between the proportions
    colnames(calculated.z.score)[1] <- "calculated.z.score" ## Fix the column name
    
    ## Now find the p-values by sampling R's built-in normal distribution. All previous results were just stepping stones on the way to this one.
    p.value <- 2 * pnorm(as.matrix(calculated.z.score), lower.tail = T)
    colnames(p.value)[1] <- "p.value" ## Fix the column name
    
    ## Bind all those results together with the input parameters and return them
    results <- cbind(parameter.space, sigma.hat, calculated.z.score, p.value)
    return(results)
    }
  
  ################################################################################
  ################################################################################
  
  
  ################################################################################
  ## Do cosmetic improvements
  ################################################################################
  results <- calculate.p.values(parameter.space) ## Call the function defined above
  ## If you make the next line <= 1, it will show the horizontal line plots along the diagonal of the faceted grid
  results <- subset(results, p.value < 1) ## Wash out the lower left half of the plots, which make the faceted graphic busy
  
  ## Fix the names so the appear more attractively in the plot
  colnames(results)[which(names(results) == "pbar.opt")] <- "Optimized.Proportion"
  colnames(results)[which(names(results) == "pbar.ctrl")] <- "Control.Proportion"
  ################################################################################
  ################################################################################
  
  
  ################################################################################
  ## Draw the plot
  ################################################################################
  
  ## Now start plotting
  plot <- ggplot(results, aes(x = n.opt, y = p.value) ) +
                  geom_point(shape = 5, aes( colour = factor(n.ctrl) ) ) +
                  geom_line(aes(colour = factor(n.ctrl))) +
                  scale_x_log10(name = "Optimized Group Size") 
  #                 scale_y_continuous(limits=c(0, .25))
  #                 geom_smooth(aes(group = factor(n.ctrl), colour = factor(n.ctrl)))
                
  # Add faceting
  plot <- plot + facet_grid(Optimized.Proportion ~ Control.Proportion 
                            , scales = 'fixed'
                            , space = 'free'
                            , labeller = 'label_both'
  )
  
  ## Fix the legend name
  plot <- plot + scale_colour_discrete(name  = "Size of Control Group")
  ################################################################################
  ################################################################################
  
  ################################################################################
  ## Make filenames and directories
  ################################################################################
  ## Make a descriptive string to stick in the filename
  range.identifier <- paste(": prob range", lower.sweep.value, "to", upper.sweep.value, "by", step.size, sep =" ")
  
  ## Figure out filename to use
  filename.to.save = paste("Faceted P Values", range.identifier, '.pdf', sep = "")
  
  ## Build directories to save in
  directory.to.save.in = file.path(getwd(), "significance_charts", as.character(Sys.Date()))
  dir.create(directory.to.save.in, recursive = TRUE, showWarnings = FALSE)
  
  ################################################################################
  ################################################################################
  
  ## Display the plot
  pdf(file.path(directory.to.save.in, filename.to.save), width = 20, height = 10)
  print(plot)
  dev.off()
  return(paste("Succesfully plotted", filename.to.save))
}

## Choose some intervals to sweep and bind them into a df
a <- c(.5,.6,.02)
b <- c(0,1,.2)
c <- c(0,.5,.1)
d <- c(.75, .76, .002)
e <- c(.6,.65, .01)
df = data.frame(a,b,c,d,e)

## Apply the draw.grid function to each row of the df and print a success message for each one
print(apply(df, 2, draw.grid))
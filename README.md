# dash - Dirichlet adaptive shrinkage for compositional data with applications in smoothing


- Authors:    [Kushal K Dey](https://github.com/kkdey),   [Dongyue Xie](https://github.com/DongyueXie), [Zhengrong Xing](https://github.com/zrxing), [Matthew Stephens](http://stephenslab.uchicago.edu/)


## Description

Compositional data is observed in various settings - composition of chemicals in experiment 
under lab settings, composition of bases at a particular position of the DNA seqeunce. Often, 
instead of relative composition or compositional percentages, we observe the counts of the
compositional categories in the data. Under such a setting, we present a new strategy called 
Dirichlet Adaptive Shrinkage (*dash*) that adaptively estimates the true underlying composition
effectively. A special case of this model is the Beta adaptive shrinkage for two compositional 
categories and such an adaptive approach can also be used in adaptive multiscale modeling of 
timecourse data. We present a R package **dashr** that comprises of all these methods.

For model formulation and the model fitting algorithm along with a discussion on
choices of parameters, please check out our [vignette](vignettes/dash.Rmd)


##  Example Applications of dash 

We first provide an example application of **dash** in shrinking the base compositional probabilities in a Trasncription Factor Binding Site (TFBS) for a particular transcription factor. 

```{r,warning=FALSE,message=FALSE,fig.width=7,fig.height=7}
library(Logolas)
library(grid)
library(gridBase)
library(ecostructure)
library(ggplot2)
library(Biobase)
library(dash)
```

We provide a simulation example of a position frequency matrix (PFM). 

```{r,warning=FALSE,message=FALSE,fig.width=7,fig.height=7}

xmat <- cbind(c(5, 0, 2, 0),
              c(1, 1, 0, 1),
              c(100, 100, 50, 100),
              c(20, 50, 100, 10),
              c(10, 10, 200, 20),
              c(50, 54, 58, 53),
              c(1,1,1,3),
              c(2, 4, 1, 1))
rownames(xmat) <- c("A", "C", "G", "T")
colnames(xmat) <- paste0("pos-", 1:dim(xmat)[2])
xmat_norm <- apply(xmat, 2, function(x) return(x/sum(x)))

xmat
```

We fit the Dirichlet adaptive shrinkage (dash) model to the position frequency matrix generated above.


```{r,warning=FALSE,message=FALSE,fig.width=7,fig.height=7}
out <- dash(xmat, optmethod = "mixEM", verbose=FALSE, bf=TRUE)
```


We present the logo plot representations of the PWM matrix obtained by normalizing the sample PFM matrix and the one after applying dash. We use the R package [Logolas](kkdey.github.io/Logolas-pages) to visually represent the logos. 

```{r}
grid.newpage()
layout.rows <- 1
layout.cols <- 2
top.vp <- viewport(layout=grid.layout(layout.rows, layout.cols,
                                      widths=unit(rep(6,layout.cols), rep("null", 2)),
                                      heights=unit(c(20,50), rep("lines", 2))))

plot_reg <- vpList()
l <- 1
for(i in 1:layout.rows){
  for(j in 1:layout.cols){
    plot_reg[[l]] <- viewport(layout.pos.col = j, layout.pos.row = i, name = paste0("plotlogo", l))
    l <- l+1
  }
}


plot_tree <- vpTree(top.vp, plot_reg)

color_profile = list("type" = "per_row", 
                     "col" = RColorBrewer::brewer.pal(4,name ="Spectral"))

pushViewport(plot_tree)
seekViewport(paste0("plotlogo", 1))
logomaker(xmat_norm,color_profile = color_profile,
          frame_width = 1,
          pop_name = "pre dash PWM",
          newpage = F)

seekViewport(paste0('plotlogo',2))
logomaker(out$posmean,color_profile = color_profile,
          frame_width = 1,
          pop_name = "post dash PWM",
          newpage = F)
```

<img src="vignettes/test/dash_app.png" alt="Structure Plot" height="700" width="1000">

## Contact

For any queries or questions or bugs, please open an issue in this Github page or write to Kushal K Dey at [kshldey@gmail.com](kshldey@gmail.com). 

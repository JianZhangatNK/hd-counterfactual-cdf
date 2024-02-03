## Example
rm(list = ls())
source('counterCDF.R')
source('counterCDF_stats.R')
source('counterCDF_graph.R')
n=500
B=500
p=50
f=10
ydx = simulate_data(n, p)
y = ydx$y
X = ydx$X
D = ydx$D
y.grids = seq(-2, 2.5, by=0.5)
object0 = counterCDF(y, basisF(X), D, y.grids, B, f)
## Counterfactual CDF estimates at y.grids:
object1 = object0$Estimate
object1

## OB decomposition of distribution functions
object3 = OB_decomp(object1)
object3

## ATT
OUTPUT_STATS = counterCDF_stats(object0)
OUTPUT_STATS$ATT

## Graphs
counterCDF_graph(object0)




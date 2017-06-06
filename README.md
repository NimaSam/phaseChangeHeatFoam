# phaseChangeHeatFoam
This is a solver for boiling and condensation which is written based on OpenFOAM-220 solver (interFoam).
## key features
* track interface using VOF method which is improved with Lafaurie smoothing filter
* consider two mass transfer models (Lee and Tanasawa)
* four validation test cases are provided:  
  (i) Stefan problem, 
  (ii) two-dimensional film boiling,
  (iii) the film condensation on a horizontal plate, 
  (iv) the laminar film condensation over a vertical plate

## Installation
```bash
cd Application
./Allwmake
```

## Documentation
It includes published papers and Ph.D Thesis [in Persian] 


## References
* [N. Samkhaniani, M.R. Ansari, "Numerical simulation of bubble condensation using CF-VOF", Progress in Nuclear Energy 89, 120-131](http://www.sciencedirect.com/science/article/pii/S0149197016300269)
* [N. Samkhaniani, M.R. Ansari, "Numerical simulation of superheated vapor bubble rising in stagnant liquid ", Heat Mass Transfer (2017). doi:10.1007/s00231-017-2031-6](https://link.springer.com/article/10.1007/s00231-017-2031-6)
* [N. Samkhaniani, M.R. Ansari, "The evaluation of the diffuse interface method for phase change simulations using OpenFOAM", Heat Transfer–Asian Research. 2017;00:1–31. doi:10.1002/htj.21268](http://onlinelibrary.wiley.com/doi/10.1002/htj.21268/full)

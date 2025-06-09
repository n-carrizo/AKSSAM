
# AKSSAM
> Automatic Knot Selection in Smooth Additive Models
> https://github.com/n-carrizo/AKSSAM
> Author: Nicolás Carrizosa

This repository contains the implementation in R of AKSSAM, which is an algorithm designed to perform knot selection in generalized additive models (GAMs).
   
The files included in this repository are listed below: 

- `Functions.r`: Contains the main functions required for **AKKSAM**, implemented under the name `GAM.asplines3.wood2`.

- `Simulations/` folder:
  - `Simulations.qmd`: Contains the code to reproduce the simulation studies presented.
  - `simulations.rdata`: Pre-computed `.RData` file to avoid the computational cost of re-running `Simulations.qmd`.

- `Real Data/` folder:
  - `Australia.qmd`: Contains the code used to reproduce the analysis of the real-world scenario.

## Contact

This project was developed by *Nicolás Carrizosa* (https://github.com/n-carrizo) as part of a research project within the Universidad Carlos III de Madrid and which constitutes as well his Master's Thesis (TFM).

It benefited from the support of the grant PRTR-CNS2023, funded by MICIU/AEI /10.13039/501100011033 and by European Union NextGenerationEU/PRTR and is part of the project *Modelos Aditivos con Restricciones para la Optimización Global*.

## Disclaimer

This is a preliminary version of the algorithm. It is unstable and may exhibit issues. It has not been optimized correctly yet.

## License

This repository is licensed under the [MIT License](LICENSE).
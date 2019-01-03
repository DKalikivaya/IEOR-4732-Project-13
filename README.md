# IEOR-4732-Project-13
Model Calibration for Equities: Test few models to see which one works the best for single name stocks (APPL, MSFT etc). We are finding a model or models to fit the surface and compare different fits against each other. 

## Stocks
Apple, Microsoft, Facebook

## Option type
Call and Put

## Models
Heston, VGSSD, VGSA

## Optimisation 
Grid Search, Nelder Mead and Gradient based

# Modules
1. DataProcessing: This module reads the stoock option data and returns the market prices
2. Generic Model: This module has functions to choose generic algorithm(Eg:'Grid Search') and generic model(Eg:'Heston')
3. Calibration: This module has loss function to evaluate the model prices vs. market prices
4. Project 13: This is the main generic function to calibrate the models

## References
IEOR 4732 Model Calibration


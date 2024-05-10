# AFSC MSE Operating Model
**Joshua A. Zahner (_UAF_)**, Ben Williams (_NOAA_)

The AFSC Operating Model is a generalized, age structured, multi-sex, multi-fleet, spatially explicit fisheries simulation model. It is intended to be used by fisheries scientists working in the North Pacific to quickly build fisheries simulation models to test various hypotheses about stock dynamics. While designed with the structure of many data-rich stock assessments used for manangement in Alaska (under the North Pacfic Fisheries Management Council) in mind, it could be adapted for use in other regions. 

### General design
`afscOM` is an age structured, multi-sex, multi-fleet, spatially explicit simulation model, that can closely replicate the dynamics of a wide vareity of stocks. 

Users can define the specific demographic parameters (natural mortality, maturity, weight-at-age, selectivity, etc.) that determine stock dynamics, as well as how those parameters vary across time, age, sex, region, and fleet. The package then provides a simple `project()` function that accepts an initial population state, the demographic parameters, and additional information about future recruitment and catch levels, and simulates the population forward in time. 

Users can, optionally, also request the model to generate observations (with provided levels of observation error) from the population.

_The model has not been thoroughly tested with multiple spatial regions as of May 2023._

---
### NOAA License

Software code created by U.S. Government employees is not subject to
copyright in the United States (17 U.S.C. §105). The United
States/Department of Commerce reserve all rights to seek and obtain
copyright protection in countries other than the United States for
Software authored in its entirety by the Department of Commerce. To this
end, the Department of Commerce hereby grants to Recipient a
royalty-free, nonexclusive license to use, copy, and create derivative
works of the Software outside of the United States.

<img src="https://raw.githubusercontent.com/nmfs-general-modeling-tools/nmfspalette/main/man/figures/noaa-fisheries-rgb-2line-horizontal-small.png" height="75" alt="NOAA Fisheries">

[U.S. Department of Commerce](https://www.commerce.gov/) | [National
Oceanographic and Atmospheric Administration](https://www.noaa.gov) |
[NOAA Fisheries](https://www.fisheries.noaa.gov/)


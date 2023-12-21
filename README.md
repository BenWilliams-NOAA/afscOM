# AFSC MSE Operating Model
# General design
2 - fleet
2 - survey
2 - sex

naming conventions:  
generally use snake_case

t,...,T = time (typically annual step)  
a,...,A = age, starts at >=1, A is a plus group  
l,...,L = length_bin, L is a plus group  
s = sex, 1 or 2
 if two then 1 = female, 2 = male - consider making a 3 for unsexed (e.g., unsexed size comps)


_obs = observations  
_sd = standard deviation (or cv)  
_ind = index of 1 (use) or 0 (do not use), length T

catch_obs  
catch_sd

srv_obs  
srv_sd  
srv_ind

patc_obs  
patc_ind  
patc_iss

pltc_obs  
pltc_ind  
pltc_iss

pats_obs  
pats_ind  
pats_iss

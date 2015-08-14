#First, generate a list of all variables in the ntuple that have a nonzero RMS and their ranges using stackergen. Variables with 0 RMS across all distributions are warned about and removed by adding a '#' to the front
#We do this twice: Once to write out the variables on a linear set of axes and once more taking the base-10 log of the x axis (This means crop can find maxima more efficiently for some variables)
stackergen inputfiles.txt allvarslin.txt 20 lin lin
stackergen inputfiles.txt allvarslog.txt 20 log lin

#Next merge the two output files and remove cuts for which a log scale makes no sense (PID can be -ve)
cat allvarslog.txt | grep -v "PID" > allvars.txt
cat allvarslin.txt | grep -v "PT" >> allvars.txt

#now strip out variables that we don't want to optimise on or plot. In this case I only want LOKI variables as I can use them in DaVinci:
cat allvars.txt | grep LOKI > goodvars.txt

#There is no point cutting differently on identical final states: mergevars turns kaonplus_X and kaonminus_X variables into min(kaonplus_X,kaonminus_X) and max(kaonplus_X,kaonminus_X). The more powerful of the two will then be chosen by sepper: 
mergevars goodvars.txt goodvars_merged.txt kaonplus kaonminus

#here I'm removing some more variables that I know won't improve the optimisation, either because they're 100% correlated with another variable or because they're mass windows that I don't want to optimise on:
cat goodvars_merged.txt \
	| grep -v "D_s_LOKI_M" \
	| grep -v "D_s_LOKI_MM" \
	| grep -v "phi_LOKI_M" \
	| grep -v "phi_LOKI_MM" \
	| grep -v "_LOKI_M)" \
	| grep -v "_LOKI_MM)" \
	> goodvars_merged_final.txt
#note: In your analysis you will probably want to do something like this when you understand your variables better. For the moment you can forget about the cuts I remove above and not worry about it.

#Sepper orders cuts based on how likely they are to have good discrimination. If a variable is very similar in both signal and background the separation power is small. If the two distributions differ by a large amount the separation power will be large. Variables with little or no separation power are warned about and removed by adding a '#' to the front
sepper inputfiles.txt goodvars_merged_final.txt goodvars_merged_sorted.txt

#Next we skim off the 10 most powerful variables from the top of the output file: 
head -21 goodvars_merged_sorted.txt > bestvars.txt

#Care should be taken with correlated varibles: Crop handles these sensibly, but it's worth taking a look at the correlation matrix and the 2D scatters to see if you're trying to optimise redundant vars. Look at the output in  correlations/Correlations.pdf and correlations/correlations_pageXXXX.pdf 
# The number at the end (0.9) means "if there's a linear correlation between two variables of greater than 0.9, kill one of them"
corr inputfiles.txt bestvars.txt bestvars.txt 4 4 2048 1536 correlations colz 0.9

#Stacker plots your variables stacked or layered in signal and background, multiple entries per page and saves them to the folder histograms. I've normalised the entries to unit area in inputfilesnorm.txt and we want layered plots to look for overlap in the signal/background distributions:
stacker inputfilesnorm.txt correlations/uniquevars.txt 4 4 2048 1536 histograms hist,nostack

#eff plots the efficiency in each variable of two datasets: One with a cut applied, one without. In the file inputfilescut.txt I've aplpied a cut to all the input datasets. This is a good way of looking at the acceptance effects of a cut on other distributions. The plots are output to the folder efficiencies
eff inputfilescut.txt inputfiles.txt correlations/uniquevars.txt 4 4 2048 1536 efficiencies e1

#varstocuts turns your variable list into a cut list that crop understands
varstocuts inputfiles.txt correlations/uniquevars.txt bestcuts.txt

#finally, crop is run to optimise your cuts
crop inputfiles.txt bestcuts.txt 10 SB crop_output.root | tee crop_output.log

# simpletools
Simpletools: Handy command line tools for ntuple manipulation and analysis
(c) Conor Fitzpatrick 2007-2015

If you find these tools useful either in whole or in part, please cite the documentation:

simpletools: Handy command line tools for ntuple manipulation and analysis
LHCb Internal Note 2009-029
Conor Fitzpatrick, University of Edinburgh

The most recent version of simpletools, including the above documentation will always be available from
/afs/cern.ch/user/c/cofitzpa/public/simpletools/

Bug reports, feature requests, comments and patches are always welcome.
conor.fitzpatrick@cern.ch

changelog
v2.0n - Made cmake compatible with C++11, removed TCint dependency to work with ROOT 6. 

v2.0m - Fix thanks to Konstantin Schubert: "Error in <TString::Replace>: first argument out of bounds: pos = -1, Length = 12" no longer occurs when there is no slash in the path. 

v2.0l - Made stacker use legends instead of stats boxes
      - Plots ins stacker now produced in .C, .png formats too

v2.0k - Moved to use of cmake thanks to B. Couturier. 
	-To compile now: 
	 cd build 
	 cmake ../
	 make && make install

v2.0j - Added support for array variables in cropcutensemble
v2.0i - Swapped correlated variable that gets removed- first should go, not second as the first is actually lower down the list if you've run sepper on it. Also added improved makefle ROOTSYS path suggested by B. Couturier

v2.0h - Punzi FoM added at request of Stephen Ogilvy and John Beddow. Significance is hard-coded at 3.0, make sure you change it to what you want in include/cropFoMs.h before compiling. To use, comment out cropFoM_SoRSB and comment in cropFoM_Punzi in src/crop.cc

v2.0g - rangefinder.cc added. Feed it a textfile containing some non-standard cutstrings and it will determine the operating range for them, writing a new varfile. 

v2.0f - subtle bug in crop found (thanks to D. Hill): cuts weren't correctly being determined as "greater than" (>) or "less than" (<) after ordering. Fixed.

v2.0e - mergevars can now merge 3 vars and is sped up 

v2.0d - Added the ability to corr to ignore 100% correlated variables. It now outputs an additional file that contains only variables that are less than 100% correlated. 

v2.0c - Added a new tool, tuplescrambler. This takes as input an ntuple and randomly shuffles it row-wise based on a random seed. This is handy for toy studies where sequential sampling is needed but the input ntuple may not be fully random. 
      - Added a check to tools that alter ntuples so that the paths in the final ntuple are the same as those in the original (no more doubling up of the ntuple directory structure)

v2.0a - MAJOR REWRITE
	should be backwards compatible, if not let me know
	-all code is now OO, file lists for crop/stacker are unified. To update them run updatedatstore on the old files. 
	-crop now can use a beta_s FoM. Just alter the cropoptimisationengine line in crop.cc to use cropFoM_betas instead prior to compiling. 
	-new tools! 
		-varstocuts turns a stacker varfile into a crop cutfile automatically 
		-sepper computes the separation power of the contents of a varfile, 
		-corr produces 2D scatters from 2 varfiles 
		-mergevars turns variables X and Y into min(X,Y) and max(X,Y) in order to handle multiple instances of the same final state sensibly
	-makefile modified to handle shared libs. It's a bit hackish though.

v1.0j - crop can now use per event weights or sWeights. Instead of specifying a weight for each file, you can specify a string composed of ntuple variables, eg: (signal_sWeight*0.05) use this by giving the command line argument -s (-sb for batch mode)
      - stacker, stackergen and eff can also handle sWeights

v1.0h - stacker, stackergen now require an additional argument to be present in the filelist: preprocessing cuts in standard cut format. For ntuples you don't want to apply a preproc. cut to, just put "1" in there instead of something more complex like "pt>100". 
	stackergen: Variables read from the ttrees are now checked to make sure they have a non-zero RMS. If all trees have a zero rms, they aren't really worth plotting, so the var. line is output with a preceeding '#' meaning stacker ignores it. 
	crop: Added a patch submitted by Gareth Rogers that fixes creation of the temporary preprocessing files: Instead of creating files in /tmp as was default, it looks to the $TMPDIR variable to define its temporary directory. 

v1.0g - added a tool "kstest.cc" to compute the Kolmogorov p-value matrices of datasets passed to stacker. Uses the same filelist and varlist as stacker

v1.0f - added a check to ensure the specified ntuple exists in cutapplier. Thanks to Ruggero Turra for spotting this.

v1.0e - added ability to disable gui for batchmode/headless operation. Add '-b' as the last argument to crop

v1.0d - crop: added a plot of the evolution of S/sqrt(S+B) as a function of the number of cut reoptimisations, 
	as requested by Ulrich Kerzel. 

v1.0c - crop: added updating gui to show the S/sqrt(S+B), RoC, signal eff and background rejection of the last optimised cut.
	these plots are now written as a single canvas to the output file. 

v1.0b - added verbose error messages to parsing code of crop and stacker: 
	Misconfigured config files now tell you where they are misconfigured instead of just segfaulting.  
	Thanks to Gareth Rogers for this suggestion. 

v1.0a - initial public release

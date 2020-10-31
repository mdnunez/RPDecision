# s2afc_RTjoyplots.R - Plots overlapping reaction time distributions
#
# Copyright (C) 2018 Michael D. Nunez, <mdnunez1@uci.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

## Record of Revisions
#   Date           Programmers               Description of change
#   ====        =================            =====================
#  11/20/17        Michael Nunez                 Original code
#  11/29/17        Michael Nunez                 Color edits
#  11/30/17        Michael Nunez          Load nondecision time posteriors
#  09/17/18        Michael Nunez         Regenerate figure based on new model


## Necessary packages
library(ggplot2)
library(ggjoy)
library(viridis)
library(R.matlab)

loadloc = '/home/michael/data10/Stroke2AFC/young'
posloc = '/home/michael/data10/Stroke2AFC/jagsout'
saveloc = '/home/michael/Dropbox/Research/s2afc/Paper'

## Code

# Read in the reaction times
reactiontimes = readMat(paste(loadloc,
  '/reactiontime.mat',sep=""))
action_selection = reactiontimes[1]$rt[1,1,1]$as
execution_only = reactiontimes[1]$rt[2,1,1]$eo

permissing_as = sum(!is.finite(action_selection))/length(action_selection)
permissing_eo = sum(!is.finite(execution_only))/length(execution_only)

# Read the non-decision time posterior distributions
jagsout = readMat(paste(posloc,'/jagsmodelbasic_lapse_8_17_16_51.mat',sep=""))


nobs = dim(action_selection)[2]
nsubs = dim(action_selection)[1]

singletrial <- data.frame(matrix(ncol = 3, nrow = 2*nsubs*nobs))
colnames(singletrial) <- c("RT", "DT", "Condition")

sublabs <- vector()
letters = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N")

for (sub in seq(1,nsubs)) {
	singletrial$RT[((sub-1)*nobs+1):(sub*nobs)] = action_selection[sub,]
	singletrial$Condition[((sub-1)*nobs+1):(sub*nobs)] = sprintf("%dSub AS",sub)
}
for (sub in seq(1,nsubs)) {
	singletrial$RT[((sub-1 + nsubs)*nobs+1):((sub+nsubs)*nobs)] = execution_only[sub,]
	singletrial$Condition[((sub-1 + nsubs)*nobs+1):((sub+nsubs)*nobs)] = sprintf("%dSub EO",sub)
  sublabs[(sub*2 -1):(sub*2)] = c(letters[sub],letters[sub])
}

# Create Reaction Time plot

gg <- ggplot(singletrial,aes(x=RT,y=Condition,fill=Condition)) +
  scale_fill_cyclical(values = c("blue","red"), guide="legend",
    labels = c("Action Selection (AS)", "Execution Only (EO)")) +
  geom_joy(scale=5, alpha=.67) + xlim(0,2000) + 
  xlab('Reaction Time (ms)') + ylab('Participant') +
  scale_y_discrete(labels=sublabs) +
  theme(axis.text=element_text(size=14),axis.title=element_text(size=14,face='bold'), legend.text = element_text(size=14),
    legend.title=element_text(size=14,face='bold'))

nondecision_pos <- matrix(ncol = nsubs*2, nrow = 15000)

# Plot non-decision time posteriors on top of reaction time distributions

track = 1
for (sub in seq(1,nsubs)) {
  for (cond in seq(1,2)) {
    temp = as.vector(jagsout[4]$chains[56+track,1,1]$tersub)
    nondecision_pos[,track] = temp
    tempmin = quantile(temp,.025)*1000
    tempmax = quantile(temp,.975)*1000
    tempmedian = quantile(temp,.5)*1000
    ypos = (sub>=10)*(-18) + (sub<10)*(10) + track
    # ypos = track
    if (cond == 1) {
      color = 'black'
      singletrial$DT[((sub-1)*nobs+1):(sub*nobs)] = action_selection[sub,] - tempmedian
    } else {
      color = 'black'
      singletrial$DT[((sub-1 + nsubs)*nobs+1):((sub+nsubs)*nobs)] = execution_only[sub,] - tempmedian
    }
    gg = gg + annotate("segment", x=tempmin, y=ypos, xend= tempmax, yend= ypos, colour = color,size=3,alpha=.5)
    # gg = gg + annotate("text", x=tempmedian-10, y=ypos-.5, label="*", size=10, alpha=.5, colour= color)
    track = track + 1
  }
}

# Create Decision Time plot

gg2 <- ggplot(singletrial,aes(x=DT,y=Condition,fill=Condition)) +
  scale_fill_cyclical(values = c("blue","red"), guide="legend",
    labels = c("Action Selection (AS)", "Execution Only (EO)")) +
  geom_joy(scale=5, alpha=.67) + xlim(0,2000) + 
  xlab('Decision Time (ms)') + ylab('Participant') +
  scale_y_discrete(labels=sublabs) +
  theme(axis.text=element_text(size=14),axis.title=element_text(size=14,face='bold'), legend.text = element_text(size=14),
    legend.title=element_text(size=14,face='bold'))

## Save plots

png(paste(saveloc,'/reaction_time_distributions.png',sep=""),units="in",width=10,height=10,res=300)
plot(gg)
dev.off()

png(paste(saveloc,'/decision_time_distributions.png',sep=""),units="in",width=10,height=10,res=300)
plot(gg2)
dev.off()


par(mfrow = c(2,5))

ost1 = read.table("OST1.single_mut.txt", header = F)
ost1_bar = barplot(ost1$V2/sum(ost1$V2), ylim = c(0, 0.3), ylab = "Proportion of polymorphic sites", main = "OST1")
text(x = ost1_bar, y= ost1$V2/sum(ost1$V2), label = paste(ost1$V1), pos = 3, cex = 0.8, col = "red")

GDP1 = read.table("GDP1.single_mut.txt", header = F)
GDP1_bar = barplot(GDP1$V2/sum(GDP1$V2), ylim = c(0, 0.3), ylab = "Proportion of polymorphic sites", main = "GDP1")
text(x = GDP1_bar, y= GDP1$V2/sum(GDP1$V2), label = paste(GDP1$V1), pos = 3, cex = 0.8, col = "red")

PFY1 = read.table("PFY1.single_mut.txt", header = F)
PFY1_bar = barplot(PFY1$V2/sum(PFY1$V2), ylim = c(0, 0.3), ylab = "Proportion of polymorphic sites", main = "PFY1")
text(x = PFY1_bar, y= PFY1$V2/sum(PFY1$V2), label = paste(PFY1$V1), pos = 3, cex = 0.8, col = "red")

RNR1 = read.table("RNR1.single_mut.txt", header = F)
RNR1_bar = barplot(RNR1$V2/sum(RNR1$V2), ylim = c(0, 0.3), ylab = "Proportion of polymorphic sites", main = "RNR1")
text(x = RNR1_bar, y= RNR1$V2/sum(RNR1$V2), label = paste(RNR1$V1), pos = 3, cex = 0.8, col = "red")

RNR2 = read.table("RNR2.single_mut.txt", header = F)
RNR2_bar = barplot(RNR2$V2/sum(RNR2$V2), ylim = c(0, 0.3), ylab = "Proportion of polymorphic sites", main = "RNR2")
text(x = RNR2_bar, y= RNR2$V2/sum(RNR2$V2), label = paste(RNR2$V1), pos = 3, cex = 0.8, col = "red")

STM1 = read.table("STM1.single_mut.txt", header = F)
STM1_bar = barplot(STM1$V2/sum(STM1$V2), ylim = c(0, 0.3), ylab = "Proportion of polymorphic sites", main = "STM1")
text(x = STM1_bar, y= STM1$V2/sum(STM1$V2), label = paste(STM1$V1), pos = 3, cex = 0.8, col = "red")

TDH1 = read.table("TDH1.single_mut.txt", header = F)
TDH1_bar = barplot(TDH1$V2/sum(TDH1$V2), ylim = c(0, 0.3), ylab = "Proportion of polymorphic sites", main = "TDH1")
text(x = TDH1_bar, y= TDH1$V2/sum(TDH1$V2), label = paste(TDH1$V1), pos = 3, cex = 0.8, col = "red")

TDH2 = read.table("TDH2.single_mut.txt", header = F)
TDH2_bar = barplot(TDH2$V2/sum(TDH2$V2), ylim = c(0, 0.3), ylab = "Proportion of polymorphic sites", main = "TDH2")
text(x = TDH2_bar, y= TDH2$V2/sum(TDH2$V2), label = paste(TDH2$V1), pos = 3, cex = 0.8, col = "red")

TDH3 = read.table("TDH3.single_mut.txt", header = F)
TDH3_bar = barplot(TDH3$V2/sum(TDH3$V2), ylim = c(0, 0.3), ylab = "Proportion of polymorphic sites", main = "TDH3")
text(x = TDH3_bar, y= TDH3$V2/sum(TDH3$V2), label = paste(TDH3$V1), pos = 3, cex = 0.8, col = "red")

VMA7 = read.table("VMA7.single_mut.txt", header = F)
VMA7_bar = barplot(VMA7$V2/sum(VMA7$V2), ylim = c(0, 0.3), ylab = "Proportion of polymorphic sites", main = "VMA7")
text(x = VMA7_bar, y= VMA7$V2/sum(VMA7$V2), label = paste(VMA7$V1), pos = 3, cex = 0.8, col = "red")

STM1_freq = read.table("STM1.freq", header = F)
TDH1_freq = read.table("TDH1.freq", header = F)
TDH2_freq = read.table("TDH2.freq", header = F)
TDH3_freq = read.table("TDH3.freq", header = F)
GDP1_freq = read.table("GDP1.freq", header = F)
VMA7_freq = read.table("VMA7.freq", header = F)
OST1_freq = read.table("OST1.freq", header = F)
RNR1_freq = read.table("RNR1.freq", header = F)
RNR2_freq = read.table("RNR2.freq", header = F)
PFY1_freq = read.table("PFY1.freq", header = F)

par(mfrow=c(1,1))
boxplot(OST1_freq$V2, GDP1_freq$V2, PFY1_freq$V2, RNR1_freq$V2, RNR2_freq$V2, STM1_freq$V2, TDH1_freq$V2, TDH2_freq$V2, TDH3_freq$V2, VMA7_freq$V2, names = c("OST1", "GDP1", "PFY1", "RNR1", "RNR2", "STM1", "TDH1", "TDH2", "TDH3", "VMA7"), outline=FALSE, ylab = "Alternative allele frequency")

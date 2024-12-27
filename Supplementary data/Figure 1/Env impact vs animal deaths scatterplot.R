# Scatterplot for the CO2 emissions of animal products (Poore & Nemecek, 2018) and the animal welfare impacts (Faunalytics, 2022)
# (c) Soemano Zeijlmans, 2024

#Create PNG output
png("EnvImpactvsAnimalDeaths.png", width = 1600, height = 1800, res=300)

#Specify font
windowsFonts(A = windowsFont("Garamond"))

#Create dataframe with impact for both externalities
gfg_data <- data.frame(x =c(84.4554505,9.86582363,3.151719611,23.87758122,12.30568178,4.669470953,13.63244514,26.86586234), #environmental impact
                       y = c(0.00575,1.07232,0.00011,0.00111,0.83369,0.29087,5.95405,249.47945), #total lives/kg
                       lab=c('Beef','Chicken','Milk','Cheese','Pork','Eggs','Farmed fish','Farmed shrimp')) 
gfg_data       

#Plot the data
plot(gfg_data$x,                                 
     gfg_data$y,
     log = "y",
     xlab = expression("CO"[2]*"-equivalent emissions (kg/kg)"), 
     ylab = "Total animal deaths (animals/kg, log scale)",
     main = "Greenhouse gas emissions and animal\ndeaths for different animal proteins",
     family = "A",
     font=1,
     ylim = c(1e-04, 200),
     xlim = c(0, 90)) 

#Add text
text(gfg_data$x,                               
     gfg_data$y, 
     labels = gfg_data$lab, 
     pos = 4,
     family = "A",
     font=1,
     cex = 0.7) 

#Close PNG connection
dev.off()

#End of script

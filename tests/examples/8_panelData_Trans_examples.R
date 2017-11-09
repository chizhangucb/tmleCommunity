#***************************************************************************************
# Examples of Panel Data Transformation
data(Grunfeld, package = "plm")
data(comSample.wmT.bA.bY_list)
pData <- comSample.wmT.bA.bY_list$comSample.wmT.bA.bY
#***************************************************************************************

# 1. Fixed effect transformation ("within") where individual effect is introduced 
pData.FE <- panelData_Trans(data = pData, xvar = c("E1", "E2", "W1", "W2", "W3", "A"),
                            yvar = "Y", index = "id", effect = "individual", 
                            model = "within", transY = TRUE)
# "E1", E2" and "A" are removed since they are constant in community level
names(pData.FE)  # "Y" "W1" "W2" "W3" 
head(pData.FE)

# 2. Same as example 1 but not transforming the outcome variable 
pData.FE.fixY <- panelData_Trans(data = pData, xvar = c("E1", "E2", "W1", "W2", "W3", "A"),
                                 yvar = "Y", index = "id", effect = "individual", 
                                 model = "within", transY = FALSE)
all.equal(pData.FE.fixY$Y, pData$Y)  # TRUE

# 3. Same as example 2 but different yvar and xvar 
pData.FE.fixY.2 <- panelData_Trans(data = pData, xvar = c("E1", "E2", "W1", "W2", "W3"),
                                   yvar = "A", index = "id", transY = FALSE)

# 4. Pooled OLS transformation ("pooling") where individual effect is introduced 
pData.pool <- panelData_Trans(data = pData, xvar = c("E1", "E2", "W1", "W2", "W3", "A"),
                              yvar = "Y", index = "id", effect = "individual", 
                              model = "pooling", transY = TRUE)
names(pData.pool)  # Y" "(Intercept)" "E1" "E2" "W1" "W2" "W3" "A"      

# 5. Random effect transformation ("random") where time effect is introduced 
Grunfeld.RE <- panelData_Trans(yvar = "inv", xvar = c("value", "capital"), data = Grunfeld, 
                               effect = "time", model = "random", index = "year")

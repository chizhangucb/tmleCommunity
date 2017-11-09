#***************************************************************************************
# Examples of Panel Data Transformation
data(comSample.wmT.bA.bY_list)
pData <- comSample.wmT.bA.bY_list$comSample.wmT.bA.bY
#***************************************************************************************

#***************************************************************************************
# 1.1 Using "within" model (fixed effect transformation) 
#***************************************************************************************
pData.FE <- panelData_Trans(data = pData, xvar = c("E1", "E2", "W1", "W2", "W3", "A"),
                            yvar = "Y", index = "id", effect = "individual", 
                            model = "random", transY = TRUE)
# "E1", E2" and "A" are removed since they are constant in community level
names(comSample.wmT.bA.bY.FE)  # "Y" "W1" "W2" "W3" 
head(comSample.wmT.bA.bY.FE)

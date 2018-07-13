\dontrun{
  
# load sample stData object
data("coprecip")
attach(coprecip)

# compute correlations directly from stData object
cors = teleCor(coprecip)
str(cors)

# manually specify data matrices (result is identical to above)
cors = teleCor(Y = Y, Z = Z, coords.s = coords.s, coords.r = coords.r)

# view pointwise correlations between Denver and ocean locations
plot(cors, coords.s = c(-105, 39.73))

}

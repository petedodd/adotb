## making a country key
library(here)
library(data.table)

UN <- c("Angola",
        "Bangladesh",
        "Brazil",
        "Central African Republic",
        "China",
        "Congo",
        "Democratic People's Republic of Korea",
        "Democratic Republic of the Congo",
        "Ethiopia",
        "Gabon",
        "India",
        "Indonesia",
        "Kenya",
        "Lesotho",
        "Liberia",
        "Mongolia",
        "Mozambique",
        "Myanmar",
        "Namibia",
        "Nigeria",
        "Pakistan",
        "Papua New Guinea",
        "Philippines",
        "Sierra Leone",
        "South Africa",
        "Thailand",
        "Uganda",
        "United Republic of Tanzania",
        "Viet Nam",
        "Zambia")

ihme <- c(
  "Angola",
  "Bangladesh",
  "Brazil",
  "Central African Republic",
  "China",
  "Congo",
  "Democratic People's Republic of Korea",
  "Democratic Republic of the Congo",
  "Ethiopia",
  "Gabon",
  "India",
  "Indonesia",
  "Kenya",
  "Lesotho",
  "Liberia",
  "Mongolia",
  "Mozambique",
  "Myanmar",
  "Namibia",
  "Nigeria",
  "Pakistan",
  "Papua New Guinea",
  "Philippines",
  "Sierra Leone",
  "South Africa",
  "Thailand",
  "Uganda",
  "United Republic of Tanzania",
  "Viet Nam",
  "Zambia"
)


isoz <- c("AGO","BGD","BRA",
          "CAF","CHN","COG",
          "PRK","COD","ETH",
          "GAB","IND","IDN",
          "KEN","LSO","LBR",
          "MNG","MOZ","MMR",
          "NAM","NGA","PAK",
          "PNG","PHL","SLE",
          "ZAF","THA","UGA",
          "TZA","VNM","ZMB")



ckey <- data.table(iso3=isoz,ihme=ihme,UN=UN)

## ckey[,.(iso3,ihme)]
## ckey[,.(iso3,UN)]

ckey[,newcountry:=UN]
ckey[UN=="Democratic People's Republic of Korea",newcountry:="DPR Korea"]
ckey[UN=="Democratic Republic of the Congo",newcountry:="DR Congo"]
ckey[UN=="United Republic of Tanzania",newcountry:="UR Tanzania"]
ckey[UN=="Papua New Guinea",newcountry:="PNG"]


save(ckey,file=here('rawdata/ckey.Rdata'))

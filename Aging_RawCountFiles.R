#This is the collection of Muscle Aging Raw Count Files to load in
#From https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE226117
#The files are saved on my computer at /Users/brussm/Documents/RStudioProjects/Rapa_PwR/GSE226117_RAW
#Load files into RStudio environment and proceed with analysis.
#The files are from several experimental Groups.
#1. Female c57bl/6 mice ages (6mo, 12mo, 18mo, 21mo, 24mo, 27mo)
    #Gastrocnemius (gas), Soleus (sol), Tibialis Anterior (ta) and Diaphragm (dia)
#2. Male c57bl/6 mice...
#3. Male Wistar Rats
#---------------------------#
getwd()

#////////////////////////////////#
#-----Female Gastrocnemius--------#
#////////////////////////////////#
file_list_gas <- c(
  # 6mo_gas
  "GSM7064465_X1385216_raw_counts.txt", "GSM7064466_X1385217_raw_counts.txt",
  "GSM7064467_X1385218_raw_counts.txt", "GSM7064468_X1385219_raw_counts.txt",
  "GSM7064469_X1385220_raw_counts.txt", "GSM7064470_X1385221_raw_counts.txt",
  "GSM7064471_X1385303_raw_counts.txt", "GSM7064472_X1385304_raw_counts.txt",
  
  # 12mo_gas
  "GSM7064421_X1385228_raw_counts.txt", "GSM7064422_X1385229_raw_counts.txt",
  "GSM7064423_X1385230_raw_counts.txt", "GSM7064424_X1385231_raw_counts.txt",
  "GSM7064425_X1385232_raw_counts.txt", "GSM7064426_X1385233_raw_counts.txt",
  "GSM7064427_X1385305_raw_counts.txt", "GSM7064428_X1385306_raw_counts.txt",
  
  # 18mo_gas
  "GSM7064429_X1385240_raw_counts.txt", "GSM7064430_X1385241_raw_counts.txt",
  "GSM7064431_X1385242_raw_counts.txt", "GSM7064432_X1385243_raw_counts.txt",
  "GSM7064433_X1385244_raw_counts.txt", "GSM7064434_X1385245_raw_counts.txt",
  "GSM7064435_X1385309_raw_counts.txt", "GSM7064436_X1385310_raw_counts.txt",
  
  # 21mo_gas
  "GSM7064437_X1385260_raw_counts.txt", "GSM7064438_X1385261_raw_counts.txt",
  "GSM7064439_X1385262_raw_counts.txt", "GSM7064440_X1385263_raw_counts.txt",
  "GSM7064441_X1385264_raw_counts.txt", "GSM7064442_X1385265_raw_counts.txt",
  "GSM7064443_X1385266_raw_counts.txt", "GSM7064444_X1385267_raw_counts.txt",
  
  # 24mo_gas
  "GSM7064445_X1385272_raw_counts.txt", "GSM7064446_X1385273_raw_counts.txt",
  "GSM7064447_X1385274_raw_counts.txt", "GSM7064448_X1385275_raw_counts.txt",
  "GSM7064449_X1385276_raw_counts.txt", "GSM7064450_X1385277_raw_counts.txt",
  "GSM7064451_X1385278_raw_counts.txt", "GSM7064452_X1385335_raw_counts.txt",
  "GSM7064453_X1385336_raw_counts.txt", "GSM7064454_X1385337_raw_counts.txt",
  
  # 27mo_gas
  "GSM7064455_X1385279_raw_counts.txt", "GSM7064456_X1385280_raw_counts.txt",
  "GSM7064457_X1385281_raw_counts.txt", "GSM7064458_X1385282_raw_counts.txt",
  "GSM7064459_X1385283_raw_counts.txt", "GSM7064460_X1385292_raw_counts.txt",
  "GSM7064461_X1385293_raw_counts.txt", "GSM7064462_X1385294_raw_counts.txt",
  "GSM7064463_X1385360_raw_counts.txt", "GSM7064464_X1385361_raw_counts.txt"
)

sample_groups_gas <- list(
  "6" = c("X1385216", "X1385217", "X1385218", "X1385219",
          "X1385220", "X1385221", "X1385303", "X1385304"),
  "12" = c("X1385228", "X1385229", "X1385230", "X1385231",
           "X1385232", "X1385233", "X1385305", "X1385306"),
  "18" = c("X1385240", "X1385241", "X1385242", "X1385243",
           "X1385244", "X1385245", "X1385309", "X1385310"),
  "21" = c("X1385260", "X1385261", "X1385262", "X1385263",
           "X1385264", "X1385265", "X1385266", "X1385267"),
  "24" = c("X1385272", "X1385273", "X1385274", "X1385275",
           "X1385276", "X1385277", "X1385278", "X1385335",
           "X1385336", "X1385337"),
  "27" = c("X1385279", "X1385280", "X1385281", "X1385282",
           "X1385283", "X1385292", "X1385293", "X1385294",
           "X1385360", "X1385361")
)


#////////////////////////////////#
#-----Female Soleus------------#
#////////////////////////////////#

file_list_sol <- c(
  # 6mo_sol
  "GSM7064516_X1384990_raw_counts.txt", "GSM7064517_X1384991_raw_counts.txt",
  "GSM7064518_X1384992_raw_counts.txt", "GSM7064519_X1384993_raw_counts.txt",
  "GSM7064520_X1384994_raw_counts.txt", "GSM7064521_X1384995_raw_counts.txt",
  "GSM7064522_X1385286_raw_counts.txt", "GSM7064523_X1385287_raw_counts.txt",
  
  # 12mo_sol
  "GSM7064473_X1385002_raw_counts.txt", "GSM7064474_X1385003_raw_counts.txt",
  "GSM7064475_X1385004_raw_counts.txt", "GSM7064476_X1385005_raw_counts.txt",
  "GSM7064477_X1385006_raw_counts.txt", "GSM7064478_X1385007_raw_counts.txt",
  "GSM7064479_X1385288_raw_counts.txt", "GSM7064480_X1385289_raw_counts.txt",
  
  # 18mo_sol
  "GSM7064481_X1385014_raw_counts.txt", "GSM7064482_X1385015_raw_counts.txt",
  "GSM7064483_X1385016_raw_counts.txt", "GSM7064484_X1385017_raw_counts.txt",
  "GSM7064485_X1385018_raw_counts.txt", "GSM7064486_X1385019_raw_counts.txt",
  "GSM7064487_X1385313_raw_counts.txt", "GSM7064488_X1385314_raw_counts.txt",
  
  # 21mo_sol
  "GSM7064489_X1385034_raw_counts.txt", "GSM7064490_X1385035_raw_counts.txt",
  "GSM7064491_X1385036_raw_counts.txt", "GSM7064492_X1385037_raw_counts.txt",
  "GSM7064493_X1385038_raw_counts.txt", "GSM7064494_X1385039_raw_counts.txt",
  "GSM7064495_X1385040_raw_counts.txt", "GSM7064496_X1385041_raw_counts.txt",
  
  # 24mo_sol
  "GSM7064497_X1385046_raw_counts.txt", "GSM7064498_X1385048_raw_counts.txt",
  "GSM7064499_X1385049_raw_counts.txt", "GSM7064500_X1385050_raw_counts.txt",
  "GSM7064501_X1385051_raw_counts.txt", "GSM7064502_X1385052_raw_counts.txt",
  "GSM7064503_X1385344_raw_counts.txt", "GSM7064504_X1385345_raw_counts.txt",
  "GSM7064505_X1385346_raw_counts.txt",
  
  # 27mo_sol
  "GSM7064506_X1385053_raw_counts.txt", "GSM7064507_X1385054_raw_counts.txt",
  "GSM7064508_X1385055_raw_counts.txt", "GSM7064509_X1385056_raw_counts.txt",
  "GSM7064510_X1385057_raw_counts.txt", "GSM7064511_X1385058_raw_counts.txt",
  "GSM7064512_X1385059_raw_counts.txt", "GSM7064513_X1385060_raw_counts.txt",
  "GSM7064514_X1385366_raw_counts.txt", "GSM7064515_X1385367_raw_counts.txt"
)

sample_groups_sol <- list(
  "6" = c("X1384990", "X1384991", "X1384992", "X1384993",
          "X1384994", "X1384995", "X1385286", "X1385287"),
  "12" = c("X1385002", "X1385003", "X1385004", "X1385005",
           "X1385006", "X1385007", "X1385288", "X1385289"),
  "18" = c("X1385014", "X1385015", "X1385016", "X1385017",
           "X1385018", "X1385019", "X1385313", "X1385314"),
  "21" = c("X1385034", "X1385035", "X1385036", "X1385037",
           "X1385038", "X1385039", "X1385040", "X1385041"),
  "24" = c("X1385046", "X1385048", "X1385049", "X1385050",
           "X1385051", "X1385052", "X1385344", "X1385345",
           "X1385346"),
  "27" = c("X1385053", "X1385054", "X1385055", "X1385056",
           "X1385057", "X1385058", "X1385059", "X1385060",
           "X1385366", "X1385367")
)


#////////////////////////////////#
#-----Female Tibialis Anterior---#
#////////////////////////////////#

file_list_ta <- c(
  # 6mo_ta
  "GSM7064567_X1385067_raw_counts.txt", "GSM7064568_X1385068_raw_counts.txt",
  "GSM7064569_X1385069_raw_counts.txt", "GSM7064570_X1385070_raw_counts.txt",
  "GSM7064571_X1385071_raw_counts.txt", "GSM7064572_X1385072_raw_counts.txt",
  "GSM7064573_X1385319_raw_counts.txt", "GSM7064574_X1385320_raw_counts.txt",
  
  # 12mo_ta
  "GSM7064524_X1385079_raw_counts.txt", "GSM7064525_X1385080_raw_counts.txt",
  "GSM7064526_X1385081_raw_counts.txt", "GSM7064527_X1385082_raw_counts.txt",
  "GSM7064528_X1385083_raw_counts.txt", "GSM7064529_X1385084_raw_counts.txt",
  "GSM7064530_X1385321_raw_counts.txt", "GSM7064531_X1385322_raw_counts.txt",
  
  # 18mo_ta
  "GSM7064532_X1385091_raw_counts.txt", "GSM7064533_X1385092_raw_counts.txt",
  "GSM7064534_X1385093_raw_counts.txt", "GSM7064535_X1385094_raw_counts.txt",
  "GSM7064536_X1385095_raw_counts.txt", "GSM7064537_X1385096_raw_counts.txt",
  "GSM7064538_X1385326_raw_counts.txt",
  
  # 21mo_ta
  "GSM7064539_X1385111_raw_counts.txt", "GSM7064540_X1385112_raw_counts.txt",
  "GSM7064541_X1385113_raw_counts.txt", "GSM7064542_X1385114_raw_counts.txt",
  "GSM7064543_X1385115_raw_counts.txt", "GSM7064544_X1385116_raw_counts.txt",
  "GSM7064545_X1385117_raw_counts.txt", "GSM7064546_X1385118_raw_counts.txt",
  
  # 24mo_ta
  "GSM7064547_X1385123_raw_counts.txt", "GSM7064548_X1385124_raw_counts.txt",
  "GSM7064549_X1385125_raw_counts.txt", "GSM7064550_X1385126_raw_counts.txt",
  "GSM7064551_X1385127_raw_counts.txt", "GSM7064552_X1385128_raw_counts.txt",
  "GSM7064553_X1385129_raw_counts.txt", "GSM7064554_X1385353_raw_counts.txt",
  "GSM7064555_X1385354_raw_counts.txt", "GSM7064556_X1385355_raw_counts.txt",
  
  # 27mo_ta
  "GSM7064557_X1385130_raw_counts.txt", "GSM7064558_X1385131_raw_counts.txt",
  "GSM7064559_X1385132_raw_counts.txt", "GSM7064560_X1385133_raw_counts.txt",
  "GSM7064561_X1385134_raw_counts.txt", "GSM7064562_X1385135_raw_counts.txt",
  "GSM7064563_X1385136_raw_counts.txt", "GSM7064564_X1385137_raw_counts.txt",
  "GSM7064565_X1385372_raw_counts.txt", "GSM7064566_X1385373_raw_counts.txt"
)

sample_groups_ta <- list(
  "6" = c("X1385067", "X1385068", "X1385069", "X1385070",
          "X1385071", "X1385072", "X1385319", "X1385320"),
  "12" = c("X1385079", "X1385080", "X1385081", "X1385082",
           "X1385083", "X1385084", "X1385321", "X1385322"),
  "18" = c("X1385091", "X1385092", "X1385093", "X1385094",
           "X1385095", "X1385096", "X1385326"),
  "21" = c("X1385111", "X1385112", "X1385113", "X1385114",
           "X1385115", "X1385116", "X1385117", "X1385118"),
  "24" = c("X1385123", "X1385124", "X1385125", "X1385126",
           "X1385127", "X1385128", "X1385129", "X1385353",
           "X1385354", "X1385355"),
  "27" = c("X1385130", "X1385131", "X1385132", "X1385133",
           "X1385134", "X1385135", "X1385136", "X1385137",
           "X1385372", "X1385373")
)
#------------------------------------------------------------#

#////////////////////////////////#
#-----Male Gastrocnemius--------#
#////////////////////////////////#

file_list_m_gas <- c(
  # 6mo_m_gas
  "GSM7064667_X1385222_raw_counts.txt", "GSM7064668_X1385223_raw_counts.txt",
  "GSM7064669_X1385224_raw_counts.txt", "GSM7064670_X1385225_raw_counts.txt",
  "GSM7064671_X1385226_raw_counts.txt", "GSM7064672_X1385227_raw_counts.txt",
  "GSM7064673_X1385301_raw_counts.txt", "GSM7064674_X1385302_raw_counts.txt",
  
  # 12mo_m_gas
  "GSM7064623_X1385234_raw_counts.txt", "GSM7064624_X1385235_raw_counts.txt",
  "GSM7064625_X1385236_raw_counts.txt", "GSM7064626_X1385237_raw_counts.txt",
  "GSM7064627_X1385238_raw_counts.txt", "GSM7064628_X1385239_raw_counts.txt",
  "GSM7064629_X1385307_raw_counts.txt", "GSM7064630_X1385308_raw_counts.txt",
  
  # 18mo_m_gas
  "GSM7064631_X1385246_raw_counts.txt", "GSM7064632_X1385247_raw_counts.txt",
  "GSM7064633_X1385248_raw_counts.txt", "GSM7064634_X1385249_raw_counts.txt",
  "GSM7064635_X1385250_raw_counts.txt", "GSM7064636_X1385251_raw_counts.txt",
  "GSM7064637_X1385311_raw_counts.txt", "GSM7064638_X1385312_raw_counts.txt",
  
  # 21mo_m_gas
  "GSM7064639_X1385252_raw_counts.txt", "GSM7064640_X1385253_raw_counts.txt",
  "GSM7064641_X1385254_raw_counts.txt", "GSM7064642_X1385255_raw_counts.txt",
  "GSM7064643_X1385256_raw_counts.txt", "GSM7064644_X1385257_raw_counts.txt",
  "GSM7064645_X1385258_raw_counts.txt", "GSM7064646_X1385259_raw_counts.txt",
  
  # 24mo_m_gas
  "GSM7064647_X1385268_raw_counts.txt", "GSM7064648_X1385269_raw_counts.txt",
  "GSM7064649_X1385270_raw_counts.txt", "GSM7064650_X1385271_raw_counts.txt",
  "GSM7064651_X1385329_raw_counts.txt", "GSM7064652_X1385330_raw_counts.txt",
  "GSM7064653_X1385331_raw_counts.txt", "GSM7064654_X1385332_raw_counts.txt",
  "GSM7064655_X1385333_raw_counts.txt", "GSM7064656_X1385334_raw_counts.txt",
  
  # 27mo_m_gas
  "GSM7064657_X1385295_raw_counts.txt", "GSM7064658_X1385296_raw_counts.txt",
  "GSM7064659_X1385297_raw_counts.txt", "GSM7064660_X1385298_raw_counts.txt",
  "GSM7064661_X1385299_raw_counts.txt", "GSM7064662_X1385300_raw_counts.txt",
  "GSM7064663_X1385356_raw_counts.txt", "GSM7064664_X1385357_raw_counts.txt",
  "GSM7064665_X1385358_raw_counts.txt", "GSM7064666_X1385359_raw_counts.txt"
)

sample_groups_m_gas <- list(
  "6" = c("X1385222", "X1385223", "X1385224", "X1385225",
          "X1385226", "X1385227", "X1385301", "X1385302"),
  "12" = c("X1385234", "X1385235", "X1385236", "X1385237",
           "X1385238", "X1385239", "X1385307", "X1385308"),
  "18" = c("X1385246", "X1385247", "X1385248", "X1385249",
           "X1385250", "X1385251", "X1385311", "X1385312"),
  "21" = c("X1385252", "X1385253", "X1385254", "X1385255",
           "X1385256", "X1385257", "X1385258", "X1385259"),
  "24" = c("X1385268", "X1385269", "X1385270", "X1385271",
           "X1385329", "X1385330", "X1385331", "X1385332",
           "X1385333", "X1385334"),
  "27" = c("X1385295", "X1385296", "X1385297", "X1385298",
           "X1385299", "X1385300", "X1385356", "X1385357",
           "X1385358", "X1385359")
)
#------------------------------------------------------#

#////////////////////////////////#
#-----Male Soleus------------#
#////////////////////////////////#

file_list_m_sol <- c(
  # 6mo_m_sol
  "GSM7064719_X1384996_raw_counts.txt", "GSM7064720_X1384997_raw_counts.txt",
  "GSM7064721_X1384998_raw_counts.txt", "GSM7064722_X1384999_raw_counts.txt",
  "GSM7064723_X1385000_raw_counts.txt", "GSM7064724_X1385001_raw_counts.txt",
  "GSM7064725_X1385284_raw_counts.txt", "GSM7064726_X1385285_raw_counts.txt",
  
  # 12mo_m_sol
  "GSM7064675_X1385008_raw_counts.txt", "GSM7064676_X1385009_raw_counts.txt",
  "GSM7064677_X1385010_raw_counts.txt", "GSM7064678_X1385011_raw_counts.txt",
  "GSM7064679_X1385012_raw_counts.txt", "GSM7064680_X1385013_raw_counts.txt",
  "GSM7064681_X1385290_raw_counts.txt", "GSM7064682_X1385291_raw_counts.txt",
  
  # 18mo_m_sol
  "GSM7064683_X1385020_raw_counts.txt", "GSM7064684_X1385021_raw_counts.txt",
  "GSM7064685_X1385022_raw_counts.txt", "GSM7064686_X1385023_raw_counts.txt",
  "GSM7064687_X1385024_raw_counts.txt", "GSM7064688_X1385025_raw_counts.txt",
  "GSM7064689_X1385315_raw_counts.txt", "GSM7064690_X1385316_raw_counts.txt",
  
  # 21mo_m_sol
  "GSM7064691_X1385026_raw_counts.txt", "GSM7064692_X1385027_raw_counts.txt",
  "GSM7064693_X1385028_raw_counts.txt", "GSM7064694_X1385029_raw_counts.txt",
  "GSM7064695_X1385030_raw_counts.txt", "GSM7064696_X1385031_raw_counts.txt",
  "GSM7064697_X1385032_raw_counts.txt", "GSM7064698_X1385033_raw_counts.txt",
  
  # 24mo_m_sol
  "GSM7064699_X1385042_raw_counts.txt", "GSM7064700_X1385043_raw_counts.txt",
  "GSM7064701_X1385044_raw_counts.txt", "GSM7064702_X1385045_raw_counts.txt",
  "GSM7064703_X1385338_raw_counts.txt", "GSM7064704_X1385339_raw_counts.txt",
  "GSM7064705_X1385340_raw_counts.txt", "GSM7064706_X1385341_raw_counts.txt",
  "GSM7064707_X1385342_raw_counts.txt", "GSM7064708_X1385343_raw_counts.txt",
  
  # 27mo_m_sol
  "GSM7064709_X1385061_raw_counts.txt", "GSM7064710_X1385062_raw_counts.txt",
  "GSM7064711_X1385063_raw_counts.txt", "GSM7064712_X1385064_raw_counts.txt",
  "GSM7064713_X1385065_raw_counts.txt", "GSM7064714_X1385066_raw_counts.txt",
  "GSM7064715_X1385362_raw_counts.txt", "GSM7064716_X1385363_raw_counts.txt",
  "GSM7064717_X1385364_raw_counts.txt", "GSM7064718_X1385365_raw_counts.txt"
)

sample_groups_m_sol <- list(
  "6" = c("X1384996", "X1384997", "X1384998", "X1384999",
          "X1385000", "X1385001", "X1385284", "X1385285"),
  "12" = c("X1385008", "X1385009", "X1385010", "X1385011",
           "X1385012", "X1385013", "X1385290", "X1385291"),
  "18" = c("X1385020", "X1385021", "X1385022", "X1385023",
           "X1385024", "X1385025", "X1385315", "X1385316"),
  "21" = c("X1385026", "X1385027", "X1385028", "X1385029",
           "X1385030", "X1385031", "X1385032", "X1385033"),
  "24" = c("X1385042", "X1385043", "X1385044", "X1385045",
           "X1385338", "X1385339", "X1385340", "X1385341",
           "X1385342", "X1385343"),
  "27" = c("X1385061", "X1385062", "X1385063", "X1385064",
           "X1385065", "X1385066", "X1385362", "X1385363",
           "X1385364", "X1385365")
)
#---------------------------------------------------#

#////////////////////////////////#
#-----Female Tibialis Anterior---#
#////////////////////////////////#

file_list_m_ta <- c(
  # 6mo_m_ta
  "GSM7064770_X1385073_raw_counts.txt", "GSM7064771_X1385074_raw_counts.txt",
  "GSM7064772_X1385075_raw_counts.txt", "GSM7064773_X1385076_raw_counts.txt",
  "GSM7064774_X1385077_raw_counts.txt", "GSM7064775_X1385078_raw_counts.txt",
  "GSM7064776_X1385317_raw_counts.txt", "GSM7064777_X1385318_raw_counts.txt",
  
  # 12mo_m_ta
  "GSM7064727_X1385085_raw_counts.txt", "GSM7064728_X1385086_raw_counts.txt",
  "GSM7064729_X1385087_raw_counts.txt", "GSM7064730_X1385088_raw_counts.txt",
  "GSM7064731_X1385089_raw_counts.txt", "GSM7064732_X1385090_raw_counts.txt",
  "GSM7064733_X1385323_raw_counts.txt", "GSM7064734_X1385324_raw_counts.txt",
  
  # 18mo_m_ta
  "GSM7064735_X1385097_raw_counts.txt", "GSM7064736_X1385098_raw_counts.txt",
  "GSM7064737_X1385099_raw_counts.txt", "GSM7064738_X1385100_raw_counts.txt",
  "GSM7064739_X1385102_raw_counts.txt", "GSM7064740_X1385327_raw_counts.txt",
  "GSM7064741_X1385328_raw_counts.txt",
  
  # 21mo_m_ta
  "GSM7064742_X1385103_raw_counts.txt", "GSM7064743_X1385104_raw_counts.txt",
  "GSM7064744_X1385105_raw_counts.txt", "GSM7064745_X1385106_raw_counts.txt",
  "GSM7064746_X1385107_raw_counts.txt", "GSM7064747_X1385108_raw_counts.txt",
  "GSM7064748_X1385109_raw_counts.txt", "GSM7064749_X1385110_raw_counts.txt",
  
  # 24mo_m_ta
  "GSM7064750_X1385119_raw_counts.txt", "GSM7064751_X1385120_raw_counts.txt",
  "GSM7064752_X1385121_raw_counts.txt", "GSM7064753_X1385122_raw_counts.txt",
  "GSM7064754_X1385347_raw_counts.txt", "GSM7064755_X1385348_raw_counts.txt",
  "GSM7064756_X1385349_raw_counts.txt", "GSM7064757_X1385350_raw_counts.txt",
  "GSM7064758_X1385351_raw_counts.txt", "GSM7064759_X1385352_raw_counts.txt",
  
  # 27mo_m_ta
  "GSM7064760_X1385138_raw_counts.txt", "GSM7064761_X1385139_raw_counts.txt",
  "GSM7064762_X1385140_raw_counts.txt", "GSM7064763_X1385141_raw_counts.txt",
  "GSM7064764_X1385142_raw_counts.txt", "GSM7064765_X1385143_raw_counts.txt",
  "GSM7064766_X1385368_raw_counts.txt", "GSM7064767_X1385369_raw_counts.txt",
  "GSM7064768_X1385370_raw_counts.txt", "GSM7064769_X1385371_raw_counts.txt"
)

sample_groups_m_ta <- list(
  "6" = c("X1385073", "X1385074", "X1385075", "X1385076",
          "X1385077", "X1385078", "X1385317", "X1385318"),
  "12" = c("X1385085", "X1385086", "X1385087", "X1385088",
           "X1385089", "X1385090", "X1385323", "X1385324"),
  "18" = c("X1385097", "X1385098", "X1385099", "X1385100",
           "X1385102", "X1385327", "X1385328"),
  "21" = c("X1385103", "X1385104", "X1385105", "X1385106",
           "X1385107", "X1385108", "X1385109", "X1385110"),
  "24" = c("X1385119", "X1385120", "X1385121", "X1385122",
           "X1385347", "X1385348", "X1385349", "X1385350",
           "X1385351", "X1385352"),
  "27" = c("X1385138", "X1385139", "X1385140", "X1385141",
           "X1385142", "X1385143", "X1385368", "X1385369",
           "X1385370", "X1385371")
)


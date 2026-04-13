library(tidyverse)
library(wqr)
library(kableExtra)
library(rmarkdown)
library(knitr)
library(lubridate)
library(stringr)
library(blastula)


sclass <- c("Routine Daily", "THMs/HAAs Monthly", "Violation Check Samples")
param <- c("Chlorine Residual Total", "Field-Chlorine Residual Total")
end_date <- floor_date(Sys.Date(), "month") #floor date rounds the current date down to the first of the month
start_date <- end_date - months(1)
dbp_start_date <- end_date - years(1) - months(1)
pm_start <- end_date - months(2)
year <- isoyear(start_date)
percentgrabsamples_start_date <- floor_date(Sys.Date(), "year") - years(7)
currentmonth <- (end_date - months(1))
pmonth <- currentmonth - months(1)
currentmonthlabel <- format(as.Date(currentmonth), '%B')
previousmonthlabel <- format(as.Date(pmonth), '%B')

#install.packages("wqr)

# wqr::get_site_info()
# wqr::get_site_info() -> site_info
# site_info$loc_id
# wqr::get_site_info(in_service == 1)

#partnership sites(from kyles word doc)
# fcrt_sites <- c("1300","1201", "1301", "2400", "2600", "3800", "1101", "1102", "1103", 
#                 "1104", "1106", "1107", "1301", "1302", "1303", "1304", "1305", "1401",
#                 "1601", "1602", "1603", "1604", "1608", "1701", "1702", "1704", "1705",
#                 "1707", "1708", "1710", "1711", "1712", "1713", "1715", "1716", "1718", 
#                 "1719", "1720", "1721", "2401", "2402", "2403", "2501", "2502", "2601", 
#                 "2602", "2603", "2701", "2703", "2704", "2706", "2707", "2708", 
#                 "2709", "2712", "2713", "2716", "2718", "3801", "3803", "3804", "3901",
#                 "3904", "3905", "3906", "3907", "3909", "3911", "3912", "3913",
#                 "3914", "7101", "7207", "7301", "7302", "7401", "7502", "7601")
# fcrt_sites <- c(fcrt_sites, paste0(fcrt_sites, c("U")), paste0(fcrt_sites, c("D")))
# 
# crt_sites <- c("4001", "4005", "5004", "6001", "6002", "7102", "7103", "7200", "7201",
#                "7300", "7303")
# crt_sites <- c(crt_sites, paste0(crt_sites, c("U")), paste0(crt_sites, c("D")))
# 
# rmvd_sites <- c("2702", "3910")
# rmvd_sites <- c(rmvd_sites, paste0(rmvd_sites, c("U")), paste0(rmvd_sites, c("D")))
# 
# allsites <-  c(fcrt_sites, crt_sites, rmvd_sites)
# 

#from SRA_db
sradb <- read_sradb_table("si_site_info")

allsites <- sradb$loc_id
in_service <- sradb$loc_id[which(sradb$in_service==TRUE)]
rtcr <- c(sradb$loc_id[which(sradb$rtcr_site==TRUE)], paste0(sradb$loc_id[which(sradb$rtcr_site==TRUE)], c("U")), paste0(sradb$loc_id[which(sradb$rtcr_site==TRUE)], c("D")))

high_service_sites <- c(4001, 4005, 5004, 6001, 6002)
high_service_sites <- c(high_service_sites, paste0(high_service_sites, c("U")), paste0(high_service_sites, c("D")))
dbp_sites <- sradb$loc_id[which(sradb$dbp_site==TRUE)]

#all sites prior to covid--------------------------------------------------------
precovidsites <- c(1101, 1102, 1103, 1104, 1106, 1107, 1202, 1301, 1302, 1303, 1304,
                   1305, 1401, 1601, 1602, 1603, 1604, 1608, 1701, 1702, 1703, 1704, 
                   1705, 1707, 1708, 1709, 1710, 1711, 1712, 1713, 1715, 1716, 1717,
                   1718, 1719, 1720, 1721, 2401, 2402, 2403, 2501, 2502, 2601, 2602,
                   2603, 2701, 2702, 2703, 2704, 2706, 2707, 2708, 2709, 2712, 2713, 
                   2715, 2716, 2718, 3801, 3803, 3804, 3901, 3902, 3904, 3905, 3906, 
                   3907, 3908, 3909, 3910, 3911, 3912, 3913, 3914, 4001, 4004, 4005, 
                   4006, 4103, 4201, 4202, 4501, 4502, 4503, 4903, 5004, 5103, 5201, 
                   5202, 5501, 5502, 5903, 6001, 6002, 6103, 6201, 6202, 6501, 6502,
                   6903, 7101, 7102, 7103, 7104, 7200, 7201, 7204, 7205, 7207, 7300, 
                   7301, 7302, 7303, 7401, 7502, 7601, 8602)

precovid_non_comp <- c(4001, 4004, 4005, 4006, 4103, 4201, 4202, 4501, 4502, 4503,
                       4903, 5004, 5103, 5201, 5202, 5501, 5502, 5903, 6001, 6002, 
                       6103, 6201, 6202, 6501, 6502, 6903, 7102, 7103, 7104, 7200, 
                       7201, 7204, 7205, 7300, 7303)

precovid_oos <- c(1202, 1703, 1709, 1717, 2702, 2715, 3902, 3908, 3910, 7104, 7204, 7205,
                  8602)

precovid_dep_comp_sites <- setdiff(precovidsites, precovid_non_comp)
precovid_dep_comp_sites <- setdiff(precovid_dep_comp_sites, precovid_oos) #sites used for compliance prior to covid
precovid_dep_comp_sites <- c(precovid_dep_comp_sites, paste0(precovid_dep_comp_sites, c("U")), paste0(precovid_dep_comp_sites, c("D")))
# including the sites we no longer sample from will not change counts as long as we do not sample from them



#current month df---------------------------------------------------------------

currentmonth <- read_LIMS(
  start_date = start_date,
  end_date = end_date,
  site = allsites,
  sample_class = sclass,
  parameter = c(param, "Coliforms Total (Colilert)")
) %>% filter(is.na(project_no) | project_no == "THM_HAA Monthly")



countcurrent <- currentmonth %>% 
  pivot_wider(
    id_cols = c(lims_number, site, date_time),
    names_from = parameter,
    values_from = result
  ) %>% 
  rename(
    'field_chlorine'='Field-Chlorine Residual Total', 'lab_chlorine' = 'Chlorine Residual Total', 'total_coliforms' = 'Coliforms Total (Colilert)'
  ) %>% 
  filter(
    site %in% c(rtcr) & !is.na(field_chlorine)|
      site %in% c(setdiff(allsites, rtcr)) & !is.na(lab_chlorine)|
      site %in% c(rtcr) & is.na(field_chlorine) & !is.na(lab_chlorine)
  )%>% 
  mutate(
    type = ifelse(site %in% c(rtcr) & !is.na(field_chlorine) & !is.na(total_coliforms), "Field-Chlorine Residual Total/RTCR", 
                  ifelse(site %in% c(rtcr) & !is.na(field_chlorine) & is.na(total_coliforms), "Field-Chlorine Residual Total Only",
                  ifelse(site %in% c(setdiff(allsites, rtcr)) & !is.na(lab_chlorine), "Chlorine Residual Total",
                         ifelse(site %in% rtcr & is.na(field_chlorine) & !is.na(lab_chlorine), "z", NA)))),
    
    result_to_use = ifelse(type == "z", lab_chlorine,
                           ifelse(type == "Field-Chlorine Residual Total/RTCR", field_chlorine,
                                  ifelse( type == "Chlorine Residual Total", lab_chlorine,
                                          ifelse( type == "Field-Chlorine Residual Total Only", field_chlorine, NA)))) 
  ) %>% 
  group_by(type
  ) %>% 
  summarise(count = n(),
            count_lessthan0.5 = sum(result_to_use<0.5),
            count_lessthan1 = sum(result_to_use<1)
  )

countcurrent <- countcurrent %>% 
  mutate(
    count_lessthan0.5 = ifelse(is.na(count_lessthan0.5), 0, count_lessthan0.5),
    count_lessthan1 = ifelse(is.na(count_lessthan1), 0, count_lessthan1),
    percent = (count_lessthan0.5/count)*100
  )

countcurrent<-countcurrent %>% 
  add_row(
    type = "Total", count = sum(countcurrent$count), count_lessthan0.5 = sum(countcurrent$count_lessthan0.5), count_lessthan1 = sum(countcurrent$count_lessthan1), percent = (sum(countcurrent$count_lessthan0.5)/sum(countcurrent$count))*100
  ) %>% 
  arrange(type) 
#PROBLEM WITH CURRENT SCRIPT IS THAT IT INCLUDES OCCT IN THE Z SAMPLES. FIGURE OUT HOW TO EXCLUDE THEM. Jan(1) and sep (4) - enddate
#have to arrange table this way to pick results out of it. Just in case there is a month where there is no instance of lab chlorine ONLY at a field site
#z represents field sites(distribution_sites, fieldstorsite, covid_sites) where no field data was taken but lab data was taken.
#typically these results are caused from DBP sites (2701 specifically) that are excluded when you query for field-chlorine individually (because there was only lab taken) but are included when you query for both parameters. 


coliformcount1 <- ifelse("Field-Chlorine Residual Total Only" %in% countcurrent$type, countcurrent$count[2], 0)
coliformcount2 <- ifelse("Field-Chlorine Residual Total Only" %in% countcurrent$type, countcurrent$count[3], 0)
drr_count <- ifelse("Field-Chlorine Residual Total Only" %in% countcurrent$type, countcurrent$count[3], countcurrent$count[2]) + coliformcount1
total_count <- ifelse("Field-Chlorine Residual Total Only" %in% countcurrent$type, countcurrent$count[4], countcurrent$count[3])

numless0.5 <- ifelse("Field-Chlorine Residual Total Only" %in% countcurrent$type, countcurrent$count_lessthan0.5[4], countcurrent$count_lessthan0.5[3])
numless1 <- ifelse("Field-Chlorine Residual Total Only" %in% countcurrent$type, countcurrent$count_lessthan1[4], countcurrent$count_lessthan1[3])
percentless0.5 <- ifelse("Field-Chlorine Residual Total Only" %in% countcurrent$type, countcurrent$percent[4], countcurrent$percent[3])

phrase <- ifelse("Field-Chlorine Residual Total Only" %in% countcurrent$type,paste(paste(drr_count, "field chlorine residual samples were collected in accordance with DRR requirements."), paste(coliformcount2, "field chlorine residual samples were collected in accordance with RTCR requirements.")), paste(drr_count, "field chlorine residual samples were collected in accordance with DRR/RTCR requirements."))




#site count --------------------------------------------------------------------
num_site<-currentmonth %>% 
  group_by(site) %>% 
  summarise(
    count = n()
  )

site_count <- length(unique(num_site$site))

#previous month df--------------------------------------------------------------
previousmonth <- read_LIMS(
  start_date = pm_start,
  end_date = start_date,
  site = allsites,
  sample_class = sclass,
  parameter = param
) %>% 
  filter(
    !is.na(result),
    is.na(project_no) | project_no == "THM_HAA Monthly"
  )

#previous month counts----------------------------------------------------------
# countprevious <- previousmonth %>% 
#   filter(
#     !is.na(result),
#     site %in% c(distribution_sites, fieldstorsite, covid_sites) & parameter == "Field-Chlorine Residual Total"| site %in% c(labstorsite, high_service_sites) & parameter == "Chlorine Residual Total"
#   ) %>% 
#   group_by(parameter
#   ) %>% 
#   summarise(
#     count = n(),
#     count_lessthan0.5 = sum(result<0.5),
#     count_lessthan1 = sum(result<1)
#     )
# previoustotal <- sum(countprevious$count)


countprevious <- previousmonth %>% 
  pivot_wider(
    id_cols = c(lims_number, site, date_time),
    names_from = parameter,
    values_from = result
  ) %>% 
  rename(
    'field_chlorine'='Field-Chlorine Residual Total', 'lab_chlorine' = 'Chlorine Residual Total'
    ) %>% 
  filter(
    site %in% rtcr & !is.na(field_chlorine)|
      site %in% c(setdiff(allsites, rtcr)) & !is.na(lab_chlorine)|
      site %in% c(setdiff(allsites, rtcr)) & is.na(field_chlorine) & !is.na(lab_chlorine)
  ) %>% 
  mutate(
    type = ifelse(site %in% rtcr & !is.na(field_chlorine), "Field-Chlorine Residual Total", 
                  ifelse(site %in% c(setdiff(allsites, rtcr)) & !is.na(lab_chlorine), "Chlorine Residual Total",
                  ifelse(site %in% rtcr & is.na(field_chlorine) & !is.na(lab_chlorine), "z", NA))),
    result_to_use = ifelse(type == "zt", lab_chlorine,
                           ifelse(type == "Field-Chlorine Residual Total", field_chlorine,
                                  ifelse( type == "Chlorine Residual Total", lab_chlorine, NA))) 
     ) %>% 
  group_by(type
  ) %>% 
  summarise(count = n(),
            count_lessthan0.5 = sum(result_to_use<0.5),
            count_lessthan1 = sum(result_to_use<1)
            )
countprevious <- countprevious %>% 
  mutate(
    count_lessthan0.5 = ifelse(is.na(count_lessthan0.5), 0, count_lessthan0.5),
    count_lessthan1 = ifelse(is.na(count_lessthan1), 0, count_lessthan1),
    percent = (count_lessthan0.5/count)*100
  )

countprevious <- countprevious %>% 
  add_row(
    type = "Total", count = sum(countprevious$count), count_lessthan0.5 = sum(countprevious$count_lessthan0.5), count_lessthan1 = sum(countprevious$count_lessthan1), percent = (sum(countprevious$count_lessthan0.5)/sum(countprevious$count))*100
  ) %>% 
  arrange(type)
#have to arrange table this way to pick results out of it. Just in case there is a month where there is no instance of lab chlorine ONLY at a field site
#z stands for lab taken, not field at a field site

previous_rtcrdrr_count <- countprevious$count[2]
previous_total_count <- countprevious$count[3]
previous_numless0.5 <-  countprevious$count_lessthan0.5[3]
previous_numless1 <-  countprevious$count_lessthan1[3]
previous_percentless0.5 <-  countprevious$percent[3]

updown<- ifelse(total_count>previous_total_count, "up from", 
                ifelse(total_count<previous_total_count, "down from", 
                       "equal to"))



# #all time df--------------------------------------------------------------------
# alltime <- read_LIMS(
#   start_date = "2015-01-01",
#   end_date = end_date,
#   site = allsites,
#   sample_class = sclass,
#   parameter = param
# )
# 
# #all time table-----------------------------------------------------------------
# check <- alltime %>% 
#   filter(
#     !is.na(result),
#     site %in% c(distribution_sites, fieldstorsite, covid_sites) & parameter == "Field-Chlorine Residual Total" | site %in% c(labstorsite, high_service_sites) & parameter == "Chlorine Residual Total"
#   ) %>% 
#   group_by(year(date_time),month(date_time)) %>% 
#   summarise(count = n(),
#             count_lessthan0.5 = sum(result<0.5),
#             count_lessthan1 = sum(result<1)) %>% 
#   mutate(percent = (count_lessthan0.5/count)*100) %>% 
#   rename('year' = 'year(date_time)', 'month' = 'month(date_time)') %>% 
#   pivot_wider(
#     id_cols = year,
#     names_from = month,
#     values_from = percent
#   )
# 
# feb <- read_LIMS(
#   start_date= "2022-02-01",
#   end_date = "2022-03-01",
#   sample_class = sclass,
#   site = c(distribution_sites, fieldstorsite, covid_sites),
#   parameter = "Chlorine Residual Total",
#   project_no = "THM_HAA Monthly"
# )
#DBPs---------------------------------------------------------------------------

dbp_df<- read_LIMS(
  site = dbp_sites,
  parameter = c("5 Haloacetic acids", "Total THMs"),
  sample_class = sclass,
  start_date = dbp_start_date,
  end_date = end_date
) %>% 
  filter(!is.na(result))

#creating dbp table-------------------------------------------------------------
tthm<- dbp_df %>%
  filter(
    parameter == "Total THMs"
  ) %>% 
  mutate(
    date_time = format(as.Date(date_time, format = "%Y-%m-%d"), "%Y-%m")
  ) %>% 
  group_by(date_time) %>% 
  summarise(tthm_min = min(result, na.rm = T),
            tthm_max = max(result, na.rm = T)) %>% 
  arrange(desc(date_time))

haa5<- dbp_df %>% 
  filter(
    parameter == "5 Haloacetic acids"
  ) %>% 
  mutate(
    date_time = format(as.Date(date_time, format = "%Y-%m-%d"), "%Y-%m")
  ) %>% 
  group_by(date_time) %>% 
  summarise(haa5_min = min(result, na.rm = T),
            haa5_max = max(result, na.rm = T)) %>% 
  arrange(desc(date_time))

dbp_table<-inner_join(
  tthm, 
  haa5,
  by = "date_time"
)

dbp_kable<- dbp_table %>% 
  mutate(across(matches("tthm_min|tthm_max|haa5_min|haa5_max"), ~round(.x, 1))) %>% 
  mutate(
    tthm_min = cell_spec(tthm_min, background = ifelse(tthm_min>80, "red", ifelse(tthm_min>60, "yellow", "white"))),
    tthm_max = cell_spec(tthm_max, background = ifelse(tthm_max>80, "red", ifelse(tthm_max>60, "yellow", "white"))),
    haa5_min = cell_spec(haa5_min, background = ifelse(haa5_min>60, "red", ifelse(haa5_min>45, "yellow", "white"))),
    haa5_max = cell_spec(haa5_max, background = ifelse(haa5_max>60, "red", ifelse(haa5_max>45, "yellow", "white")))
  ) %>% 
  kbl(caption = "Monthly Minimum and Maximum Total Trihalomethanes and 5 Haloacetic Acids", col.names = c("Date", "Total THM Minimum", "Total THM Maximum", "5 HAA Minimum (ug/L)", "5 HAA Maximum (ug/L"), escape = FALSE, align = rep('c', 5))%>%
  kable_classic(full_width = T) #kable table summarizing a years worth of DBP data
#dbp blurbs---------------------------------------------------------------------

tthm_blurb <- ifelse( 
  dbp_table$tthm_min[1]>dbp_table$tthm_min[2] && dbp_table$tthm_max[1]>dbp_table$tthm_max[2], 
  "Total THM values increased compared to last month.", 
  ifelse(
    dbp_table$tthm_min[1]>dbp_table$tthm_min[2] && dbp_table$tthm_max[1]<dbp_table$tthm_max[2],
    "Total THM Minimum increased and Total THM Maximum decreased compared to last month.",
    ifelse(
      dbp_table$tthm_min[1]<dbp_table$tthm_min[2] && dbp_table$tthm_max[1]>dbp_table$tthm_max[2],
      "Total THM Minimum decreased and Total THM Maximum increased compared to last month.",
      ifelse(
        dbp_table$tthm_min[1]<dbp_table$tthm_min[2] && dbp_table$tthm_max[1]<dbp_table$tthm_max[2],
        "Total THM values decreased compared to last month."
      )
    )
  ))

haa5_blurb <- ifelse( 
  dbp_table$haa5_min[1]>dbp_table$haa5_min[2] && dbp_table$haa5_max[1]>dbp_table$haa5_max[2], 
  "Total HAA5 values increased compared to last month.", 
  ifelse(
    dbp_table$haa5_min[1]>dbp_table$haa5_min[2] && dbp_table$haa5_max[1]<dbp_table$haa5_max[2],
    "Total HAA5 Minimum increased and Total HAA5 Maximum decreased compared to last month.",
    ifelse(
      dbp_table$haa5_min[1]<dbp_table$haa5_min[2] && dbp_table$haa5_max[1]>dbp_table$haa5_max[2],
      "Total HAA5 Minimum decreased and Total HAA5 Maximum increased compared to last month.",
      ifelse(
        dbp_table$haa5_min[1]<dbp_table$haa5_min[2] && dbp_table$haa5_max[1]<dbp_table$haa5_max[2],
        "Total HAA5 values decreased compared to last month."
      )
    )
  ))

#counts for number of dbp mcl exceedances in the month------------------------------

mcl_df<-read_LIMS(
  start_date = start_date,
  end_date = end_date,
  parameter = c("Total THMs", "5 Haloacetic acids"),
  sample_class = sclass,
  site = rtcr
) #df used to get counts over the mcls. had to use a different date range then dbp_df

mcl<- mcl_df %>% 
  group_by(parameter) %>% 
  summarise(
    greater80 = sum(result>=80),
    greater60 = sum(result>=60)
  )
  
haa5exceed<-mcl$greater60[1]
tthmexceed<-mcl$greater80[2]

#summary table------------------------------------------------------------------
# 
column_dsc1 <- c(paste("Total DRR Samples for",paste(currentmonthlabel, year)), "# of samples < 1.00 mg/L", "# of samples >0.5 mg/L, <1.00 mg/L", "# of samples >0.15 mg/L, <0.50 mg/L", "# of samples <0.15 mg/L")
column_dsc2 <- c(" ","Below PWD distribution system goal of 1.00 mg/L",
                 "Meets PfSW goal of 0.50 mg/L. Does not meet PWD distribution system goal 1.00 mg/L",
                 "Meets PADEP minimum distribution residual required (0.15 mg/L). Does not meet: PfSW goal (0.50 mg/L), PWD distribution system goal (1.00 mg/L)",
                 "Does not meet: PADEP minimum distribution residual required (0.15 mg/L), PfSW goal (0.50 mg/L), PWD distribution system goal (1.00 mg/L)")

summarytable <- read_LIMS(
  site = rtcr,
  start_date = start_date,
  end_date = end_date,
  sample_class = sclass,
  parameter = "Field-Chlorine Residual Total"
) %>%
  filter(
    is.na(result)==FALSE,
    is.na(project_no)==TRUE|project_no == "THM_HAA Monthly"
  )%>%
  group_by(month(date_time)) %>%
  summarise(count_lessthan1 = sum(result<1),
            countbetween1and0.5 = sum(result<1 & result>0.5),
            countbetween0.15and0.5 = sum(result<=0.5 & result>0.15),
            count_lessthan0.15 = sum(result<=0.15))
counts_vector <- c(drr_count, summarytable$count_lessthan1[1], summarytable$countbetween1and0.5[1], summarytable$countbetween0.15and0.5[1], summarytable$count_lessthan0.15[1])

summarytable <- data.frame(column_dsc1, counts_vector, column_dsc2 ) %>%
  kbl(col.names = c("", "Count", "Description"), escape = FALSE, align = rep('c', 5))%>%
  column_spec(c(3), width = "20em") %>%
  kable_classic(c("striped", "condensed"))
# 
# perc19 <- c(0.00,0.00,0.00,0.00,0.13, 0.29, 0.51, 1.17, 2.71, 0.91, 0.00, 0.00)
# perc20 <- c(0.00,0.00,0.00,0.00, 0.00, 0.00,0.57, 1.86, 2.78, 0.65,0.00,0.00)
# perc21 <- c(0.00,0.00,0.00,0.00, 0.00, 0.00,1.85, 0.00, 0.00, 0.47, 0.00, 0.00)
# perc22 <- c(0.00, 0.00, 0.00, 0.00,0.00,0.00,0.15,1.07, 0.87,0.15,0.00,0.00)
# 
# yearlabel <- as.character(isoyear(end_date))
# otable <- data.frame(
#   month = c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"),
#   X2019 = perc19,
#   X2020 = perc20,
#   X2021 = perc21,
#   X2022 = perc22
# )
#
# # write.csv(otable, "otable.csv")
# 
# otable <- read.csv("otable.csv", colClasses = c("NULL", NA, NA, NA, NA, NA))
# 
# percentless0.5 <- "here"
# 
# dfname <- "otable"
# 
# if (currentmonthlabel == "January") {
#   otable <- within(assign(dfname, get(dfname)), assign(yearlabel, ifelse(month == "January", percentless0.5, NA))) 
# } else{
#   if(currentmonthlabel == "February"){
#     otable[2,6] <- percentless0.5
#   }else{
#     if(currentmonthlabel == "March"){
#       otable[3,6] <- percentless0.5
#     }else{
#       if(currentmonthlabel == "April"){
#         otable[4,6] <- percentless0.5
#       }else{
#         if(currentmonthlabel == "May"){
#           otable[5,6] <- percentless0.5
#         }else{
#           if(currentmonthlabel == "June"){
#             otable[6,6] <- percentless0.5
#           }else{
#             if(currentmonthlabel == "July"){
#               otable[7,6] <- percentless0.5
#             }else{
#               if(currentmonthlabel == "August"){
#                 otable[8,6] <- percentless0.5
#               }else{
#                 if(currentmonthlabel == "September"){
#                   otable[9,6] <- percentless0.5
#                 }else{
#                   if(currentmonthlabel == "October"){
#                     otable[10,6] <- percentless0.5
#                   }else{
#                     if(currentmonthlabel == "November"){
#                       otable[11,6] <- percentless0.5
#                     }else{
#                       if(currentmonthlabel == "December"){
#                         otable[2,6] <- percentless0.5
#                       }
#                     }
#                   }
#                 }
#               }
#             }
#           }
#         }
#       }
#     }
#   }
# }
# 
# write.csv(otable, "otable.csv")
# 
# 
# 
# 
# 
# df <- data.frame("A" = c(1, 2, 3, 4), "B" = c("a", "c", "d", "b") )
# Variable<-"C"
# dfname<-"df"
# df<-within ( assign(dfname  , get(dfname) ),
#              assign(yearlabel, NA          )
#            )
# 
# 
# #overall_table------------------------------------------------------------------


##have to uncomment this to add data. Commented so that I could mess around with the code without changing anything
# add_dsop_monthly_counts(
#   year = year,
#   month = month(start_date),
#   n_total = total_count, #total_count
#   n_drr = drr_count, #drr_count
#   n_less_point5 = numless0.5, #numless0.5
#   n_less_one = numless1, #numless1
#   comments =NA
# )

monthly_table <- get_dsop_monthly_counts(
  start_year = year - 7,
  end_year = year) %>% 
  arrange(desc(year), desc(month)) %>%
  mutate(month = month(month, label=TRUE),
         p_less_point5 = round(p_less_point5, digits = 2)) %>% 
  rename(Year = "year") %>% 
  pivot_wider(
    id_cols = Year,
    names_from = month,
    values_from = p_less_point5
  )

monthly_kbl <-  monthly_table %>% mutate(
  Jan = as.numeric(Jan),
  Feb = as.numeric(Feb),
  Mar = as.numeric(Mar),
  Apr = as.numeric(Apr),
  May = as.numeric(May),
  Jun = as.numeric(Jun),
  Jul = as.numeric(Jul),
  Aug = as.numeric(Aug),
  Sep = as.numeric(Sep),
  Oct = as.numeric(Oct),
  Nov = as.numeric(Nov),
  Dec = as.numeric(Dec)
)

monthly_kbl <- monthly_kbl %>% pivot_longer(cols= c('Jul', 'Jun', 'May', 'Apr', 'Mar', 'Feb', 'Jan', 'Dec', 'Nov', 'Oct', 'Sep', 'Aug'),
                                            names_to = "month",
                                            values_to = 'percent')

monthly_kbl <- monthly_kbl %>% 
  mutate(
    percent = cell_spec(percent, background = ifelse(!is.na(percent) & percent >5.0, "red", ifelse(!is.na(percent) & percent > 4.0 & percent < 4.99, "orange", ifelse(!is.na(percent) & percent >3.0 & percent <3.99, "yellow", "white"))))
  )

monthly_kbl <- monthly_kbl %>% 
  pivot_wider(
    id_cols = "Year",
    names_from = "month",
    values_from = "percent"
  ) %>%
  select("Year", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec") %>% 
  kbl(caption = "Percent of Grab Samples < 0.50 mg/L TCL2 Residual - Entire PWD Distribution System - by Month", escape = FALSE, align = rep('c', 5))%>%
  kable_classic(full_width = T) #kable table summarizing a years worth of DBP data


# monthly_kbl <- monthly_table %>% 
#   mutate(across(matches("Year|Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec"), ~round(.x, 1))) %>% 
#   mutate(
#     Jan = cell_spec(Jan, background = ifelse(Jan>5.0, "red", ifelse(Jan>4.0 & Jan<4.99, "orange", ifelse(Jan>3.0 & Jan<3.99, "yellow", "white")))),
#     Feb = cell_spec(Feb, background = ifelse(Feb>5.0, "red", ifelse(Feb>4.0 & Feb<4.99, "orange", ifelse(Feb>3.0 & Feb<3.99, "yellow", "white")))),
#     Mar = cell_spec(Mar, background = ifelse(Mar>5.0, "red", ifelse(Mar>4.0 & Mar<4.99, "orange", ifelse(Mar>3.0 & Mar<3.99, "yellow", "white")))),
#     Apr = cell_spec(Apr, background = ifelse(Apr>5.0, "red", ifelse(Apr>4.0 & Apr<4.99, "orange", ifelse(Apr>3.0 & Apr<3.99, "yellow", "white")))),
#     May = cell_spec(May, background = ifelse(May>5.0, "red", ifelse(May>4.0 & May<4.99, "orange", ifelse(May>3.0 & May<3.99, "yellow", "white")))),
#     Jun = cell_spec(Jun, background = ifelse(Jun>5.0, "red", ifelse(Jun>4.0 & Jun<4.99, "orange", ifelse(Jun>3.0 & Jun<3.99, "yellow", "white")))),
#     Jul = cell_spec(Jul, background = ifelse(Jul>5.0, "red", ifelse(Jul>4.0 & Jul<4.99, "orange", ifelse(Jul>3.0 & Jul<3.99, "yellow", "white")))),
#     Aug = cell_spec(Aug, background = ifelse(Aug>5.0, "red", ifelse(Aug>4.0 & Aug<4.99, "orange", ifelse(Aug>3.0 & Aug<3.99, "yellow", "white")))), 
#     Sep = cell_spec(Sep, background = ifelse(Sep>5.0, "red", ifelse(Sep>4.0 & Sep<4.99, "orange", ifelse(Sep>3.0 & Sep<3.99, "yellow", "white")))),
#     Oct = cell_spec(Oct, background = ifelse(Oct>5.0, "red", ifelse(Oct>4.0 & Oct<4.99, "orange", ifelse(Oct>3.0 & Oct<3.99, "yellow", "white")))),
#     Nov = cell_spec(Nov, background = ifelse(Nov>5.0, "red", ifelse(Nov>4.0 & Nov<4.99, "orange", ifelse(Nov>3.0 & Nov<3.99, "yellow", "white")))),
#     Dec = cell_spec(Dec, background = ifelse(Dec>5.0, "red", ifelse(Dec>4.0 & Dec<4.99, "orange", ifelse(Dec>3.0 & Dec<3.99, "yellow", "white"))))
#   ) %>% 
#   kbl(caption = "Percent of Grab Samples < 0.50 mg/L TCL2 Residual - Entire PWD Distribution System - by Month", escape = FALSE, align = rep('c', 5))%>%
#   kable_classic(full_width = T) #kable table summarizing a years worth of DBP data




# dsop_drr_df_to_add <- currentmonth %>% select(site, lims_number, sample_class, parameter, result, date_time, project_no)
# 
# dsop_drr <- get_dsop_drr_results(
# add_dsop_drr_results(dsop_drr_df_to_add)
#   start_date = "2023-01-01",
#   end_date = "2023-08-01"
# )

#Partnership Check Counts

library(tidyverse)
library(wqr)
library(kableExtra)
library(rmarkdown)
library(knitr)
library(lubridate)
library(stringr)
library(blastula)

#Vectors------------
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

distribution_sites <- c(1101, 1102, 1103, 1104, 1106, 1107, 1201, 1300, 1301, 1302,
                        1303, 1304, 1305, 1401, 1601, 1602, 1603, 1604, 1608, 1701,
                        1702, 1704, 1705, 1707, 1708, 1710, 1711, 1712, 1713, 1715, 
                        1716, 1718, 1719, 1720, 1721, 2400, 2401, 2402, 2403, 2501,
                        2502, 2600, 2601, 2602, 2603, 2701, 2703, 2704, 2706,
                        2707, 2708, 2709, 2712, 2713, 2716, 2718, 3800, 3801, 3803, 
                        3804, 3901, 3904, 3905, 3906, 3907, 3909, 3911, 3912,
                        3913, 3914
) #removed 2702, 3910. may have to include them to get accurate historical data?
distribution_sites <- c(distribution_sites, paste0(distribution_sites, c("U")), paste0(distribution_sites, c("D")))

rmvd_sites <- c(2702, 3910)
rmvd_sites <- c(rmvd_sites, paste0(rmvd_sites, c("U")), paste0(rmvd_sites, c("D")))

fieldstorsite <- c(7101, 7207, 7301, 7302, 7401, 7502, 7601)
fieldstorsite <- c(fieldstorsite, paste0(fieldstorsite, c("U")), paste0(fieldstorsite, c("D")))

labstorsite <- c(7102, 7103, 7200, 7201, 7300, 7303)
labstorsite <- c(labstorsite, paste0(labstorsite, c("U")), paste0(labstorsite, c("D")))

high_service_sites <- c(4001, 4005, 5004, 6001, 6002)
high_service_sites <- c(high_service_sites, paste0(high_service_sites, c("U")), paste0(high_service_sites, c("D")))

covid_sites <- c(1201, 1301, 2400, 2600, 3800)
covid_sites <- c(covid_sites, paste0(covid_sites, c("U")), paste0(covid_sites, c("D")))

allsites <- c(distribution_sites, fieldstorsite, labstorsite, high_service_sites, covid_sites, rmvd_sites)


#Cl2 DF----------------
currentmonth <- read_LIMS(
  start_date = start_date,
  end_date = end_date,
  site = allsites,
  sample_class = sclass,
  parameter = c(param, "Coliforms Total (Colilert)")
) %>% filter(is.na(project_no) | project_no == "THM_HAA Monthly")


countcurrent_df <- currentmonth %>% 
  pivot_wider(
    id_cols = c(lims_number, site, date_time),
    names_from = parameter,
    values_from = result
  ) %>% 
  rename(
    'field_chlorine'='Field-Chlorine Residual Total', 'lab_chlorine' = 'Chlorine Residual Total', 'total_coliforms' = 'Coliforms Total (Colilert)'
  ) %>% 
  filter(
    site %in% c(distribution_sites, fieldstorsite, covid_sites) & !is.na(field_chlorine)|
      site %in% c(labstorsite, high_service_sites) & !is.na(lab_chlorine)|
      site %in% c(distribution_sites, fieldstorsite, covid_sites) & is.na(field_chlorine) & !is.na(lab_chlorine)
  )%>% 
  mutate(
    type = ifelse(site %in% c(distribution_sites, fieldstorsite, covid_sites) & !is.na(field_chlorine) & !is.na(total_coliforms), "Field-Chlorine Residual Total/RTCR", 
                  ifelse(site %in% c(distribution_sites, fieldstorsite, covid_sites) & !is.na(field_chlorine) & is.na(total_coliforms), "Field-Chlorine Residual Total Only",
                         ifelse(site %in% c(labstorsite, high_service_sites) & !is.na(lab_chlorine), "Chlorine Residual Total",
                                ifelse(site %in% c(distribution_sites, fieldstorsite, covid_sites) & is.na(field_chlorine) & !is.na(lab_chlorine), "z", NA)))),
    
    result_to_use = ifelse(type == "z", lab_chlorine,
                           ifelse(type == "Field-Chlorine Residual Total/RTCR", field_chlorine,
                                  ifelse( type == "Chlorine Residual Total", lab_chlorine,
                                          ifelse( type == "Field-Chlorine Residual Total Only", field_chlorine, NA)))) 
  ) 


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
    site %in% c(distribution_sites, fieldstorsite, covid_sites) & !is.na(field_chlorine)|
      site %in% c(labstorsite, high_service_sites) & !is.na(lab_chlorine)|
      site %in% c(distribution_sites, fieldstorsite, covid_sites) & is.na(field_chlorine) & !is.na(lab_chlorine)
  )%>% 
  mutate(
    type = ifelse(site %in% c(distribution_sites, fieldstorsite, covid_sites) & !is.na(field_chlorine) & !is.na(total_coliforms), "Field-Chlorine Residual Total/RTCR", 
                  ifelse(site %in% c(distribution_sites, fieldstorsite, covid_sites) & !is.na(field_chlorine) & is.na(total_coliforms), "Field-Chlorine Residual Total Only",
                         ifelse(site %in% c(labstorsite, high_service_sites) & !is.na(lab_chlorine), "Chlorine Residual Total",
                                ifelse(site %in% c(distribution_sites, fieldstorsite, covid_sites) & is.na(field_chlorine) & !is.na(lab_chlorine), "z", NA)))),
    
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


#DBP DF----------
dbp_df<- read_LIMS(
  site = c(distribution_sites, fieldstorsite, covid_sites),
  parameter = c("5 Haloacetic acids", "Total THMs"),
  sample_class = sclass,
  start_date = dbp_start_date,
  end_date = end_date
) %>% 
  filter(!is.na(result))

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

dbp_current <- dbp_df %>% filter(format(as.Date(date_time), '%B')==currentmonthlabel)




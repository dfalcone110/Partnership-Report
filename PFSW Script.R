#lims_dsn
#lims_uid 
#lims_pwd 
#owqm_dsn 
#owqm_uid 
#owqm_pwd 
#ops_dsn 
#ops_uid 
#ops_pwd
#sradb_dsn
#sradb_uid
#sradb_pwd

# file.edit("~/.Renviron")


library(tidyverse)
library(wqr)
library(kableExtra)
library(rmarkdown)
library(knitr)
library(lubridate)
library(stringr)
library(blastula)
library(tinytex)
#-------------------
#Variables----------
#-------------------

sclass <- c("Routine Daily", "THMs/HAAs Monthly", "Violation Check Samples")
param <- c("Chlorine Residual Total", "Field-Chlorine Residual Total")
end_date <- as.Date("2025-10-01") #floor_date(Sys.Date(), "month") #floor date rounds the current date down to the first of the month
start_date <- end_date - months(1)
dbp_start_date <- end_date - years(1) - months(1)
pm_start <- end_date - months(2)
year <- isoyear(start_date)
percentgrabsamples_start_date <- floor_date(Sys.Date(), "year") - years(7)
currentmonth <- (end_date - months(1))
pmonth <- currentmonth - months(1)
currentmonthlabel <- format(as.Date(currentmonth), '%B')
previousmonthlabel <- format(as.Date(pmonth), '%B')

previousdate <- start_date - months(1)
previousyear <- isoyear(previousdate)

#---------------------
#Variables from sradb
#---------------------

sradb <- read_sradb_table("si_site_info")

allsites <- c(sradb$loc_id, paste0(sradb$loc_id, c("U")), paste0(sradb$loc_id, c("D")))
in_service <- sradb$loc_id[which(sradb$in_service==TRUE)]
rtcr <- c(sradb$loc_id[which(sradb$rtcr_site==TRUE)], paste0(sradb$loc_id[which(sradb$rtcr_site==TRUE)], c("U")), paste0(sradb$loc_id[which(sradb$rtcr_site==TRUE)], c("D")))

high_service_sites <- c(4001, 4005, 5004, 6001, 6002)
high_service_sites <- c(high_service_sites, paste0(high_service_sites, c("U")), paste0(high_service_sites, c("D")))
dbp_sites <- sradb$loc_id[which(sradb$dbp_site==TRUE)]

#---------------------
#Current Month dataframe
#---------------------
currentmonth <- read_LIMS(
  start_date = start_date,
  end_date = end_date,
  site = allsites,
  sample_class = sclass,
  parameter = c(param, "Coliforms Total (Colilert)")
) %>% filter(is.na(project_no) | project_no == "THM_HAA Monthly") 
#Currently excluding OCCT samples which we collect Lab chlorine at but not Field-Chlorine

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
#Chlorine Residual Total are sites not found in rtcr list that had a lab chlorine sample taken
#Field-Chlorine Residual Total/RTCR are RTCR sites where a FCl2 and matching coliform sample were taken
#Field-Chlorine Residual Only are RTCR sites where an FCL2 sample was taken without a matching coliform


coliformcount1 <- ifelse("Field-Chlorine Residual Total Only" %in% countcurrent$type, countcurrent$count[2], 0)
coliformcount2 <- ifelse("Field-Chlorine Residual Total Only" %in% countcurrent$type, countcurrent$count[3], 0)
drr_count <- ifelse("Field-Chlorine Residual Total Only" %in% countcurrent$type, countcurrent$count[3], countcurrent$count[2]) + coliformcount1
total_count <- ifelse("Field-Chlorine Residual Total Only" %in% countcurrent$type, countcurrent$count[4], countcurrent$count[3])

numless0.5 <- ifelse("Field-Chlorine Residual Total Only" %in% countcurrent$type, countcurrent$count_lessthan0.5[4], countcurrent$count_lessthan0.5[3])
numless1 <- ifelse("Field-Chlorine Residual Total Only" %in% countcurrent$type, countcurrent$count_lessthan1[4], countcurrent$count_lessthan1[3])
percentless0.5 <- ifelse("Field-Chlorine Residual Total Only" %in% countcurrent$type, countcurrent$percent[4], countcurrent$percent[3])

phrase <- ifelse("Field-Chlorine Residual Total Only" %in% countcurrent$type,paste0(paste0("In ", currentmonthlabel, ", " ,drr_count, " field chlorine residual samples were collected in accordance with DRR requirements; "), paste0(coliformcount2, " field chlorine residual samples were collected in accordance with RTCR requirements. ")), paste0("In ",currentmonthlabel, ", " ,drr_count, " field chlorine residual samples were collected in accordance with DRR/RTCR requirements."))

#----------------------------
#Getting total number of sites sampled in the month
#----------------------------

num_site<-currentmonth %>% 
  group_by(site) %>% 
  summarise(
    count = n()
  )

site_count <- length(unique(num_site$site))
#count for num sites


#----------------------------
#Previous Month dataframe
#----------------------------

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


#----------------------------
#DBPs
#----------------------------



#DBP DF----------------------

dbp_df<- read_LIMS(
  site = dbp_sites,
  parameter = c("5 Haloacetic acids", "Total THMs"),
  sample_class = sclass,
  start_date = dbp_start_date,
  end_date = end_date
) %>% 
  filter(!is.na(result))

#tthm table-----------------

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

#haa5 table------------------

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

#combined table---------------

dbp_table<-inner_join(
  tthm, 
  haa5,
  by = "date_time"
) %>% separate(date_time,
               into = c("year", "month"),
               sep = "-") %>% 
  mutate(month = month.abb[as.numeric(month)],
         date_time = paste0(year, "-",month)) %>% select(date_time, tthm_min, tthm_max, haa5_min, haa5_max)


dbp_table2 <- dbp_table %>% 
  mutate(
    tthm_min = ifelse(tthm_min>=60 & tthm_min<80, paste0("\\colorbox{yellow}{",tthm_min, "}"), ifelse(
    tthm_min>=80, paste0("\\colorbox{red}{",tthm_min, "}"), tthm_min
  )),
        tthm_max = ifelse(tthm_max>=60 & tthm_max<80, paste0("\\colorbox{yellow}{",tthm_max, "}"), ifelse(
    tthm_max>=80, paste0("\\colorbox{red}{",tthm_max, "}"), tthm_max
  )),
  haa5_min = ifelse(haa5_min>=45 & haa5_min<60, paste0("\\colorbox{yellow}{",haa5_min, "}"), ifelse(
    haa5_min>=60, paste0("\\colorbox{red}{",haa5_min, "}"), haa5_min
  )),
        haa5_max = ifelse(haa5_max>=45 & haa5_max<60, paste0("\\colorbox{yellow}{",haa5_max, "}"), ifelse(
    haa5_max>=60, paste0("\\colorbox{red}{",haa5_max, "}"), haa5_max
  ))

  )

#dbp kable-------------------

dbp_kable<- dbp_table %>%
  mutate(across(matches("tthm_min|tthm_max|haa5_min|haa5_max"), ~round(.x, 1))) %>%
  mutate(
    tthm_min = cell_spec(tthm_min, background = ifelse(tthm_min>80, "red", ifelse(tthm_min>60, "yellow", "white"))),
    tthm_max = cell_spec(tthm_max, background = ifelse(tthm_max>80, "red", ifelse(tthm_max>60, "yellow", "white"))),
    haa5_min = cell_spec(haa5_min, background = ifelse(haa5_min>60, "red", ifelse(haa5_min>45, "yellow", "white"))),
    haa5_max = cell_spec(haa5_max, background = ifelse(haa5_max>60, "red", ifelse(haa5_max>45, "yellow", "white")))
  ) %>%
  kbl(caption = "Monthly Minimum and Maximum Total Trihalomethanes and 5 Haloacetic Acids", col.names = c("Date", "Total THM Minimum", "Total THM Maximum", "5 HAA Minimum (ug/L)", "5 HAA Maximum (ug/L)"), escape = FALSE, align = rep('c', 5), booktabs =T, linesep = "") %>% 
  kable_styling(latex_options = c("HOLD_position"), position = "center")%>% 
  column_spec(1, "6 cm") %>% 
  column_spec(2:5, "2 cm")
#kable table summarizing a years worth of DBP data


#DBP statements--------------
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
  "HAA5 values increased compared to last month.", 
  ifelse(
    dbp_table$haa5_min[1]>dbp_table$haa5_min[2] && dbp_table$haa5_max[1]<dbp_table$haa5_max[2],
    "HAA5 Minimum increased and Total HAA5 Maximum decreased compared to last month.",
    ifelse(
      dbp_table$haa5_min[1]<dbp_table$haa5_min[2] && dbp_table$haa5_max[1]>dbp_table$haa5_max[2],
      "HAA5 Minimum decreased and HAA5 Maximum increased compared to last month.",
      ifelse(
        dbp_table$haa5_min[1]<dbp_table$haa5_min[2] && dbp_table$haa5_max[1]<dbp_table$haa5_max[2],
        "HAA5 values decreased compared to last month."
      )
    )
  ))

#MCL exceedances -------------

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


#Generating Chlorine Residual summary table------------

column_dsc1 <- c(paste("Total DRR Samples for",paste(currentmonthlabel, year)), "Num of samples < 1.00 mg/L", "Num of samples >0.5 mg/L, <1.00 mg/L", "Num of samples >0.15 mg/L, <0.50 mg/L", "Num of samples <0.15 mg/L")
column_dsc2 <- c(" ","Below PWD distribution system goal of 1.00 mg/L",
                 "Meets PfSW goal of 0.50 mg/L. Does not meet PWD distribution system goal 1.00 mg/L",
                 "Meets PADEP minimum distribution residual required (0.15 mg/L). Does not meet: PfSW goal (0.50 mg/L), PWD distribution system goal (1.00 mg/L)",
                 "Does not meet: PADEP minimum distribution residual required (0.15 mg/L), PfSW goal (0.50 mg/L), PWD distribution system goal (1.00 mg/L)")

summarytable1 <- read_LIMS(
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
counts_vector <- c(drr_count, summarytable1$count_lessthan1[1], summarytable1$countbetween1and0.5[1], summarytable1$countbetween0.15and0.5[1], summarytable1$count_lessthan0.15[1])

summarytable <- data.frame(column_dsc1, counts_vector, column_dsc2 ) %>%
  kbl(caption = "Sample Summary",col.names = c("-", "Count", "Description"), escape = FALSE, align = rep('c', 5), booktabs = T, linesep = "")%>%
  column_spec(c(3), width = "20em")%>%
  kable_styling(latex_options = c("HOLD_position", "striped"))


#-----------------------------------------------
#Updating DSOP DB and creating table------------
#-----------------------------------------------

add_dsop_monthly_counts(
  year = year,
  month = month(start_date),
  n_total = total_count, #total_count
  n_drr = drr_count, #drr_count
  n_less_point5 = numless0.5, #numless0.5
  n_less_one = numless1, #numless1
  comments =NA
)

monthly_table1 <- get_dsop_monthly_counts(
  start_year = year - 12,
  end_year = year) %>% 
  arrange(desc(year), desc(month)) %>% 
  mutate(date_temp = as.Date(paste0(year,"-", month,"-01" ))) %>% 
  filter(date_temp<= start_date)

# monthly_table1 <- read.csv("prep.csv")

#table with all months that failed(greater than 5)
great5 <- monthly_table1 %>% 
  filter(p_less_point5>5.0) %>% 
  arrange(desc(date_temp))

#table with greatest percentages in last 3 years
great3yr <- monthly_table1 %>% filter(date_temp > (start_date - years(3)))


failuretext <- ifelse(
  start_date == great5$date_temp[1],
  paste0("In the month of ",currentmonthlabel, " ", year, ", ", great5$p_less_point5[1], "% of samples were less than 0.5 mg/L and therefore PWD failed to meet the 95% goal for Partnership. Prior to ", currentmonthlabel, " ", year, ", PWD had not exceeded the Partnership goal since ", format(as.Date(great5$date_temp[2]), '%B'), " ", great5$year[2], " (", great5$p_less_point5[2], "%). Elevated temperatures may impact disinfectant residual in PWD's distribution system. The highest Partnership percentage in the previous three calendar years (", year-3, " - ", year, ") was ", max(great3yr$p_less_point5), "% in ", format(as.Date(great3yr$date_temp[which.max(great3yr$p_less_point5)]), '%B'), " ", great3yr$year[which.max(great3yr$p_less_point5)], "."
), 
ifelse(
  start_date != great5$date_temp[1],
  paste0("PWD has continuosly met the 95% Partnership goal from ", format(as.Date(great5$date_temp[1]), '%B'), " ", great5$year[1], " to ", currentmonthlabel," ", year, " (the partnership year begins June 1st and runs through May 31st).  Prior to ", format(as.Date(great5$date_temp[1]), '%B'), " ", great5$year[1], ", PWD had not exceeded the Partnership goal since ", format(as.Date(great5$date_temp[2]), '%B'), " ", great5$year[2], " (", great5$p_less_point5[2],"%). Elevated temperatures may impact disinfectant residual in PWDs distribution system. The highest Partnership percentage in the previous three calendar years (", year-3, " - ", year, ") was ", max(great3yr$p_less_point5), "% in ", format(as.Date(great3yr$date_temp[which.max(great3yr$p_less_point5)]), '%B'), " ", great3yr$year[which.max(great3yr$p_less_point5)], "."), "CHECK"

))



monthly_table <- monthly_table1%>% 
  filter(year > year(start_date)-7) %>% 
  mutate(month = month(month, label=TRUE),
         p_less_point5 = round(p_less_point5, digits = 2)) %>% 
  rename(Year = "year") %>% 
  pivot_wider(
    id_cols = Year,
    names_from = month,
    values_from = p_less_point5
  ) %>% mutate(
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
)%>% pivot_longer(cols= c('Jul', 'Jun', 'May', 'Apr', 'Mar', 'Feb', 'Jan', 'Dec', 'Nov', 'Oct', 'Sep', 'Aug'),
                                            names_to = "month",
                                            values_to = 'percent') %>% mutate(Year = as.character(Year))


monthly_table2 <- monthly_table %>% 
  mutate(across(where(is.numeric), ~ifelse(percent <= 3.99 & percent>=3.00, paste0("\\colorbox{yellow}{", percent, "}"), ifelse(
percent>=4.00 & percent<5.00, paste0("\\colorbox{orange}{", percent, "}"), ifelse(
  percent>= 5.00, paste0("\\colorbox{red}{", percent, "}"), percent
)

  ))))%>% 
  pivot_wider(
    id_cols = "Year",
    names_from = "month",
    values_from = "percent"
  ) %>%
  select("Year", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")







source(here::here('simulation_fns.R'))

gait_colnames <- names(read_csv(here::here("gait","gait_data","WBDSc3d_csv","WBDS15walkT07.csv")))
test_gait_data <- read_csv(here::here("gait","gait_data","WBDSc3d_csv","WBDS15walkT02.csv"), col_names = F)
names(test_gait_data) <- gait_colnames

test_gait_data_subsample <- test_gait_data %>%
  select(Time, L.PSISY, L.AnkleY)

test_gait_data_subsample %>%
  head(500) %>%
  ggplot(aes(x = L.PSISY, y= L.AnkleY)) +
  geom_point()


test_gait_data_subsample %>%
	ggplot() +
	geom_point(aes(x=Time, y = L.PSISY), color = "blue") +
	geom_point(aes(x=Time, y = L.AnkleY), color = "red")

noisy_samples <- test_gait_data_subsample[,c(2,3)] %>%
  head(500)
names(noisy_samples) <- c("x","y")
sampled_data <- noisy_samples %>% mutate(f_x = NA, f_y = NA)



gait_smooth <- spline_smooth_noisy_samples(sampled_data, noisy_samples, nbasis = 10,
                              return_type = "bspline", title = "")


unit_list <-  c("14", "15","23","37")
speed_list <- c("2","4","6","8")
sensor_list <- c("L.ASISY","L.PSISY","L.Iliac.CrestY","L.KneeY","L.AnkleY","L.HeelY")
sensor_pairs <- t(combn(sensor_list,2)) 
sensor_pairs_list <- paste(sensor_pairs[,1]," / ",sensor_pairs[,2])
colnames = c("Unit","Speed","Pair","x","y","Time")
# https://stackoverflow.com/questions/48833807/initialize-an-empty-tibble-with-column-names-and-0-rows
comp_plot_tibble <- colnames %>% purrr::map_dfc(~tibble::tibble(!!.x := double())) %>%
  mutate(Unit = as.factor(Unit), Speed = as.factor(Speed), Pair = as.factor(Pair))

for (unit in unit_list){
  for (speed in speed_list){
    csv_name <- paste0("WBDS",unit,"walkT0",speed,".csv")
    loaded_data <- read_csv(paste0("../gait/gait_data/WBDSc3d_csv/",csv_name), col_names = F)
    names(loaded_data) <- gait_colnames
    loaded_data <- loaded_data %>%
      head(300) %>%
      select(Time,sensor_list) %>%
      mutate(Unit = as.factor(unit), Speed = as.factor(speed))
    
    single_unit_tibble <- c("Pair","x","y") %>% purrr::map_dfc(~tibble::tibble(!!.x := double())) %>%
      mutate(Pair = as.factor(Pair))
    
    for (i in 1:nrow(sensor_pairs)){
      sensor_pair <- sensor_pairs[i,]
      sensor_pair_data <- loaded_data %>% 
        select(sensor_pair, Time) %>%
        rename(x = sensor_pair[1], y = sensor_pair[2]) %>%
        mutate(Pair = as.factor(sensor_pairs_list[i]), Unit = as.factor(unit), Speed = as.factor(speed))
        single_unit_tibble <- bind_rows(single_unit_tibble, sensor_pair_data)
    }
    comp_plot_tibble <- bind_rows(comp_plot_tibble, single_unit_tibble)
  }
}

ggplot(comp_plot_tibble,aes(x=x,y=y)) +
  geom_point(aes(color=Unit)) +
  facet_grid(Pair~Speed,scales = "free") + 
  theme(strip.text.y = element_text(size = 6)) +
  labs(title = "First 300 Samples of Various Sensor Pairs:LY")
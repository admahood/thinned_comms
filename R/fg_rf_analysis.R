# functional group RF analysis

source("R/a_tc_data_prep.R")
library(randomForest)
library(iml)
library(tidymodels)
species_list <- read_csv("data/species_list_alm_modified.csv") |> 
  unique() |> # one duplicate in species list
  dplyr::select(-notes) |>
  mutate(SpeciesCode = str_to_lower(SpeciesCode))

plot_level_metrics <- read_csv("data/plot_level_data.csv")

fg_cover_plot <- 
  cover_plot |>
    left_join(species_list) |>
    filter(!CodeType %in% c('Substrate', 'FWD', "Fuel in Air", "UNK")) |>
    group_by(new_visit_code, PlotTreatmentStatus, NativityL48, Lifespan, CodeType) |>
    summarise(cover = sum(cover)) |>
    ungroup() |>
    mutate(fg = str_c(str_sub(NativityL48, 1,1),
                      str_sub(Lifespan, 1,1),
                      str_sub(CodeType, 1,1)) |> str_to_upper()) |>
    tidyr::separate(new_visit_code, c("PlotCode", 'phase_adj'), "\\.", remove = F) |>
    left_join(plot_level_metrics) |>
  mutate_if(is.character, as.factor)

glimpse(fg_cover_plot)

fg_cover_plot |>
  filter(str_sub(fg, 1,1) != "U",
         str_sub(fg, 3,3) != "T",
         str_sub(fg, 2,2) != "U") |>
ggplot(aes(x=phase_adj, y=cover, fill = PlotTreatmentStatus)) +
  geom_boxplot(outliers = F) +
  facet_wrap(~fg, scales='free_y', ncol=4)

# functions =================

treeplot <- function(rfmod, yvar, datavar, fgvar, title = "treeplot"){
  pr <- iml::Predictor$new(model = rfmod,
                             data = datavar |> dplyr::filter(fg == fgvar) |> dplyr::select(-fg, -yvar) ,
                             y = datavar |> dplyr::filter(fg == fgvar) |> pull(yvar))
  
  plot(TreeSurrogate$new(predictor = pr,
                         maxdepth = 3)$tree, main = title)
}

get_rsqs <- function(drf){
  cell_split <- initial_split(drf)
  cell_train <- training(cell_split)
  cell_test  <- testing(cell_split)
  rf_mod <- 
    rand_forest(trees = 5000) %>% 
    set_engine("ranger") %>% 
    set_mode("regression")
  rf_fit <- 
    rf_mod %>% 
    fit(cover ~ ., data = cell_train)
  rf_training_pred <- 
    predict(rf_fit, cell_train) %>% 
    # Add the true outcome data back in
    bind_cols(cell_train %>% 
                select(cover))
  train_rsq <- rf_training_pred %>%                # training set predictions
    yardstick::rsq(truth = cover, .pred) |>
    dplyr::rename(train_rsq = 3)
  
  rf_testing_pred <- 
    predict(rf_fit, cell_test) %>% 
    bind_cols(cell_test %>% select(cover))
  test_rsq <- rf_testing_pred %>%                   # test set predictions
    yardstick::rsq(truth = cover, .pred) |>
    dplyr::rename(test_rsq = 3) |>
    left_join(train_rsq)
  return(test_rsq)
}
# random forests ========================
drf <- fg_cover_plot |> 
  dplyr::select(-c(1:3, 5:7), -plot, -trt_u_adj,
                -TreatmentUnit, -native_cover, -n, -invaded,
                -contains('exotic'), -contains('native')) |> 
  na.omit() |>
  mutate(yst = case_when(phase == '01_Pre' ~ -1,
                         phase == 'post0-1' ~ 1,
                         phase == 'post04-5' ~ 5,
                         phase == 'post10-11' ~ 10,
                         treated == 'not_treated' ~ 0)) |>
  dplyr::select(-phase)

glimpse(drf)
rf_nlf <- randomForest(cover ~ ., data = drf |> filter(fg == "NLF") |> dplyr::select(-fg))
varImpPlot(rf_nlf)
for(i in 1:3) get_rsqs(drf |> filter(fg == "NLF") |> dplyr::select(-fg)) |> print()
for(i in 1:3) get_rsqs(drf |> filter(fg == "NSF") |> dplyr::select(-fg)) |> print()
for(i in 1:3) get_rsqs(drf |> filter(fg == "NLG") |> dplyr::select(-fg)) |> print()
for(i in 1:3) get_rsqs(drf |> filter(fg == "NLT") |> dplyr::select(-fg)) |> print()

treeplot(rfmod = rf_nlf, yvar = 'cover', fgvar = 'NLF', datavar = drf, title = "Native Long-Lived Forb Cover")

rf_nsf <- randomForest(cover ~ ., data = drf |> filter(fg == "NSF") |> dplyr::select(-fg))
varImpPlot(rf_nsf)
treeplot(rfmod = rf_nsf, yvar = 'cover', fgvar = 'NSF', datavar = drf, title = "Native Short-Lived Forb Cover")

rf_nlg <- randomForest(cover ~ ., data = drf |> filter(fg == "NLG") |> dplyr::select(-fg))
varImpPlot(rf_nlg)
treeplot(rfmod = rf_nlg, yvar = 'cover', fgvar = 'NLG', datavar = drf, title = "Native Long-Lived Graminoid Cover")


rf_nlt <- randomForest(cover ~ ., data = drf |> filter(fg == "NLT") |> dplyr::select(-fg))
varImpPlot(rf_nlt)

treeplot(rfmod = rf_nlt, yvar = 'cover', fgvar = 'NLT', datavar = drf, title = "Native Long-Lived Tree Cover")


# maybe do seedling and sapling density as well

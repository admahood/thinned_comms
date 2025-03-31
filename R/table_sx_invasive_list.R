# make a species list table

read_csv('data/species_list_alm_modified.csv') |>
  filter(NativityL48 == 'exotic') |>
  dplyr::select(FinalName, Family) |>
  dplyr::arrange(Family) |>
  write_csv('out/table_sx_invasives.csv')

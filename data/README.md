## Data

This project uses the **Vinho Verde (wine quality)** dataset (commonly shared via Kaggle).

Because Kaggle datasets can have redistribution restrictions, this repo is set up so you **add the CSV locally**:

1. Download the dataset from Kaggle (wine quality, Vinho Verde).
2. Place the CSV at: `data/raw/vino_qualita.csv` (recommended file name without accents).

If your CSV has a different name, update the path in `R/data_prep.R` (`load_wine_data()`).

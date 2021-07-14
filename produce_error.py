import matplotlib.pyplot as plt
from virocon import (get_OMAE2020_Hs_Tz, get_OMAE2020_V_Hs, 
    GlobalHierarchicalModel, plot_dependence_functions, read_ec_benchmark_dataset)
	
file_path = "datasets/D.txt"
sample = read_ec_benchmark_dataset(file_path)


# Define the structure of the joint distribution model and fit it to the data.
dist_descriptions, fit_descriptions, semantics = get_OMAE2020_V_Hs()
model = GlobalHierarchicalModel(dist_descriptions)
model.fit(sample)

# The renaming to $\alpha$ causes an error.
par_rename = {"alpha": "$\alpha$", "beta": "$\beta$"}
plot_dependence_functions(model, semantics, par_rename)

plt.show()
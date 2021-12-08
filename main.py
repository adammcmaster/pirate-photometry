#from calibration import calibrate_data
from indexing import index_data
from schedule import generate_schedules
#from plotting import plot_folded, plot_targets


#index_data(reprocess=True)
index_data()

generate_schedules()

# Find ref stars for targets which don't have any yet

#    calibrate_data()
#    plot_targets()
#    plot_folded()
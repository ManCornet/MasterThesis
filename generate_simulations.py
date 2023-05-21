#------------------------------------------------------------------------------
#
#                            Sensibility Analysis
#
#                              Graduation work
#
#------------------------------------------------------------------------------
# @ Manon Cornet                                             @ Victor Mangeleer
#
# -----------------
#   Documentation
# -----------------
# Hey love !
#
# So, please don't rush and read carrefully these very few instructions (I know you :-)) :
#
# 1 - Add a parser in your julia code and then, modify line 65 of the code to add all the
#                           the new arguments your code can take !
#
# 2 - Add, in as lists, all the parameters you want to test just below this nice documentation !
#
# 3 - Add the missing for loops in the "Generation of scripts" section
#
# 4 - Once your new for loops are declared, go back to line 65 and give all the new parameters
#                       you need to write down in your shell script
#
# 5 - Run this python code in your terminal ! You will obtain in the same location as this code
#                              a new file with the extension .sh
#
# -------------------- YOU ARE ALMOST THERE BUT DO NOT SKIP THIS IMPORTANT STEP --------------------
#
# 6 - Before running your script and generating all the good results, you must give the autorization
#        to your script to be executed on your laptop ! This is simply achivee by typing in the
#             terminal (while being in THE SAME LOCATION AS YOUR .SH FILE !):
#
#                             chmod 744 sensitivity_analysis_manon.sh
#
# 7 - You can enjoy and relax while all your simulations are running in the background
#                                   using simply the command :
#
#                                ./sensitivity_analysis_manon.sh
#
#
# That's it ! <3
#
# --------------
#   Parameters
# --------------
#
#       All the parameters you want to test and their corresponding values
#                      can be initialized here in a list !
#

# Inclusion of EVs in load profiles
EVs = [False, True]
# Inclusion of HPs in load profiles
HPs = [False, True]
# Maximum PV capacity per load bus [MVA]
MAX_PV_CAPACITY_PER_NODE = [0.4, 0.0, 0.8, 1.6]
# Cost of the energy that is imported [k€/kWh]
IMP_ELECTRICITY_ENRG_COST = [0.3, 0.6, 0.9]
# Cost of the energy that is exported [k€/kWh]
EXP_ELECTRICITY_ENRG_COST = [0.1, 0.2, 0.3]
# DSO cost of the energy that is imported [k€/kWh]
IMP_ELECTRICITY_DSO_COST = [0.1, 0.2, 0.3]
# Grid connection cost [k€/kWh]
GRID_CONNECTION_COST = [80, 120, 160]


# ----------------------
#    Script Properties
# ----------------------
# Name of the shell script containing all the tests to run
script_name = "sensitivity_analysis_manon"

# -------------------------
#    Generation of script
# -------------------------
# Opening script
script_file = open(f"{script_name}.sh", "w")

# Adding good ass looking header
script_file.write("""#!/usr/bin/env bash
#------------------------------------------------------------------------------
#
#                           Sensibility Analysis
#
#                             Graduation work
#
#------------------------------------------------------------------------------
# @ Manon Cornet
\n""")

# Going through all parameters defined previously (more can be added with other for loops)
for l in nb_elec_lines:
    for ps in power_station_voltage:

        # Writting down new simulation test in command line in the script
        script_file.write(f"""julia --project --sysimage JuliaSysimage.so src/Bilevel/main.jl 
                            --EV {l} --power_stations {ps} \n""")

# Closing the script
script_file.close()

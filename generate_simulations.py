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
import copy

# Inclusion of EVs in load profiles
EV = [False, True]
# Inclusion of HPs in load profiles
HP = [False, True]
# Inclusion of storage in the simulations
storage = [False, True]
Storage_cost = [500, 300, 600]
# Network reconfiguration allowed or not
network_reconfig = [False, True]
# Bilevel or single model
bilevel = [True, False]

# Maximum PV capacity per load bus [MVA]
PV_CAPA = [0.4, 0.0, 0.8, 1.6]
PVC = [500, 300, 600]
# Cost of the energy that is imported [k€/kWh]
EIC = [0.3, 0.6, 0.9]
# Cost of the energy that is exported [k€/kWh]
EEC = [0.1, 0.2, 0.3]
# DSO cost of the energy that is imported [k€/kWh]
DSOEC = [0.1, 0.2, 0.3]
# Grid connection cost [k€/kWh]
GCC = [80, 120, 160]
# Weight I relaxed 
weight_I = [1.0e-2, 1.0e-3, 1.0e-1]

network_file = "network_graph_simu"
plot_file    = "plot_simu"

# Reference configuration
ref_config = [EV[0], HP[0], storage[0], Storage_cost[0], network_reconfig[0], bilevel[0], PV_CAPA[0], PVC[0], EIC[0], EEC[0], DSOEC[0], GCC[0], weight_I[0]]

# Used to make each cute line of our shell script
iterator = [EV, HP, storage, Storage_cost, network_reconfig, bilevel, PV_CAPA, PVC, EIC, EEC, DSOEC, GCC, weight_I]

radiality = ["spanning_tree"] #["single_commodity", "multi_commodity"]

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

for r in radiality:
        count = 0
        for i, variable in enumerate(iterator):
            # We skip storage by itself
                if i == 2:
                        continue
                for v, value in enumerate(variable):
                        network_file_name = network_file + str(count)
                        plot_file_name    = plot_file + str(count)
                        simu_name = "simu_"+r+"_"+str(count)
                        # Current configuration
                        curr_conf = copy.deepcopy(ref_config)

                        # If we are change storage value, we change storage to true to take modification
                        if i == 3:
                                curr_conf[2] = True


                        if (count > 0 and value != ref_config[i]) or count == 0:
                                # Applying modification
                                curr_conf[i] = value
                                #print(curr_conf)

                                # Writting down new simulation test in command line in the script
                                terminal_command = f"""julia --project src/main_bilevel.jl --radiality {r} --EV {str(curr_conf[0])} --EHP {str(curr_conf[1])} --storage {str(curr_conf[2])} --Storage_cost {str(curr_conf[3])} --network_reconfig {str(curr_conf[4])} --bilevel {str(curr_conf[5])} --PV_CAPA {str(curr_conf[6])} --PVC {str(curr_conf[7])} --EIC {str(curr_conf[8])} --EEC {str(curr_conf[9])} --DSOEC {str(curr_conf[10])} --GCC {str(curr_conf[11])} --weight_I {str(curr_conf[12])} --network_graph_name {network_file_name} --plot_file_name {plot_file_name} --simu_name {simu_name}\n"""

                                # Checking for bitches
                                terminal_command = terminal_command.replace("True", "true")
                                terminal_command = terminal_command.replace("False", "false")

                                script_file.write(terminal_command)
                                count +=1 


# Closing the script
script_file.close()

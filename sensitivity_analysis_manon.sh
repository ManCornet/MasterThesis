#!/usr/bin/env bash
#------------------------------------------------------------------------------
#
#                           Sensibility Analysis
#
#                             Graduation work
#
#------------------------------------------------------------------------------
# @ Manon Cornet

julia --project src/main_bilevel.jl --EV false --EHP false --storage false --Storage_cost 300 --network_reconfig false --bilevel true --PV_CAPA 0.4 --PVC 500 --EIC 0.3 --EEC 0.1 --DSOEC 0.1 --GCC 80 --weight_I 0.01
julia --project src/main_bilevel.jl --EV true --EHP false --storage false --Storage_cost 300 --network_reconfig false --bilevel true --PV_CAPA 0.4 --PVC 500 --EIC 0.3 --EEC 0.1 --DSOEC 0.1 --GCC 80 --weight_I 0.01
julia --project src/main_bilevel.jl --EV true --EHP false --storage false --Storage_cost 300 --network_reconfig false --bilevel true --PV_CAPA 0.4 --PVC 500 --EIC 0.3 --EEC 0.1 --DSOEC 0.1 --GCC 80 --weight_I 0.01
julia --project src/main_bilevel.jl --EV true --EHP true --storage false --Storage_cost 300 --network_reconfig false --bilevel true --PV_CAPA 0.4 --PVC 500 --EIC 0.3 --EEC 0.1 --DSOEC 0.1 --GCC 80 --weight_I 0.01
julia --project src/main_bilevel.jl --EV true --EHP true --storage true --Storage_cost 300 --network_reconfig false --bilevel true --PV_CAPA 0.4 --PVC 500 --EIC 0.3 --EEC 0.1 --DSOEC 0.1 --GCC 80 --weight_I 0.01
julia --project src/main_bilevel.jl --EV true --EHP true --storage true --Storage_cost 150 --network_reconfig false --bilevel true --PV_CAPA 0.4 --PVC 500 --EIC 0.3 --EEC 0.1 --DSOEC 0.1 --GCC 80 --weight_I 0.01
julia --project src/main_bilevel.jl --EV true --EHP true --storage true --Storage_cost 500 --network_reconfig false --bilevel true --PV_CAPA 0.4 --PVC 500 --EIC 0.3 --EEC 0.1 --DSOEC 0.1 --GCC 80 --weight_I 0.01
julia --project src/main_bilevel.jl --EV true --EHP true --storage true --Storage_cost 500 --network_reconfig false --bilevel true --PV_CAPA 0.4 --PVC 500 --EIC 0.3 --EEC 0.1 --DSOEC 0.1 --GCC 80 --weight_I 0.01
julia --project src/main_bilevel.jl --EV true --EHP true --storage true --Storage_cost 500 --network_reconfig true --bilevel true --PV_CAPA 0.4 --PVC 500 --EIC 0.3 --EEC 0.1 --DSOEC 0.1 --GCC 80 --weight_I 0.01
julia --project src/main_bilevel.jl --EV true --EHP true --storage true --Storage_cost 500 --network_reconfig true --bilevel true --PV_CAPA 0.4 --PVC 500 --EIC 0.3 --EEC 0.1 --DSOEC 0.1 --GCC 80 --weight_I 0.01
julia --project src/main_bilevel.jl --EV true --EHP true --storage true --Storage_cost 500 --network_reconfig true --bilevel false --PV_CAPA 0.4 --PVC 500 --EIC 0.3 --EEC 0.1 --DSOEC 0.1 --GCC 80 --weight_I 0.01
julia --project src/main_bilevel.jl --EV true --EHP true --storage true --Storage_cost 500 --network_reconfig true --bilevel false --PV_CAPA 0.4 --PVC 500 --EIC 0.3 --EEC 0.1 --DSOEC 0.1 --GCC 80 --weight_I 0.01
julia --project src/main_bilevel.jl --EV true --EHP true --storage true --Storage_cost 500 --network_reconfig true --bilevel false --PV_CAPA 0.0 --PVC 500 --EIC 0.3 --EEC 0.1 --DSOEC 0.1 --GCC 80 --weight_I 0.01
julia --project src/main_bilevel.jl --EV true --EHP true --storage true --Storage_cost 500 --network_reconfig true --bilevel false --PV_CAPA 0.8 --PVC 500 --EIC 0.3 --EEC 0.1 --DSOEC 0.1 --GCC 80 --weight_I 0.01
julia --project src/main_bilevel.jl --EV true --EHP true --storage true --Storage_cost 500 --network_reconfig true --bilevel false --PV_CAPA 1.6 --PVC 500 --EIC 0.3 --EEC 0.1 --DSOEC 0.1 --GCC 80 --weight_I 0.01
julia --project src/main_bilevel.jl --EV true --EHP true --storage true --Storage_cost 500 --network_reconfig true --bilevel false --PV_CAPA 1.6 --PVC 500 --EIC 0.3 --EEC 0.1 --DSOEC 0.1 --GCC 80 --weight_I 0.01
julia --project src/main_bilevel.jl --EV true --EHP true --storage true --Storage_cost 500 --network_reconfig true --bilevel false --PV_CAPA 1.6 --PVC 150 --EIC 0.3 --EEC 0.1 --DSOEC 0.1 --GCC 80 --weight_I 0.01
julia --project src/main_bilevel.jl --EV true --EHP true --storage true --Storage_cost 500 --network_reconfig true --bilevel false --PV_CAPA 1.6 --PVC 300 --EIC 0.3 --EEC 0.1 --DSOEC 0.1 --GCC 80 --weight_I 0.01
julia --project src/main_bilevel.jl --EV true --EHP true --storage true --Storage_cost 500 --network_reconfig true --bilevel false --PV_CAPA 1.6 --PVC 300 --EIC 0.3 --EEC 0.1 --DSOEC 0.1 --GCC 80 --weight_I 0.01
julia --project src/main_bilevel.jl --EV true --EHP true --storage true --Storage_cost 500 --network_reconfig true --bilevel false --PV_CAPA 1.6 --PVC 300 --EIC 0.6 --EEC 0.1 --DSOEC 0.1 --GCC 80 --weight_I 0.01
julia --project src/main_bilevel.jl --EV true --EHP true --storage true --Storage_cost 500 --network_reconfig true --bilevel false --PV_CAPA 1.6 --PVC 300 --EIC 0.9 --EEC 0.1 --DSOEC 0.1 --GCC 80 --weight_I 0.01
julia --project src/main_bilevel.jl --EV true --EHP true --storage true --Storage_cost 500 --network_reconfig true --bilevel false --PV_CAPA 1.6 --PVC 300 --EIC 0.9 --EEC 0.1 --DSOEC 0.1 --GCC 80 --weight_I 0.01
julia --project src/main_bilevel.jl --EV true --EHP true --storage true --Storage_cost 500 --network_reconfig true --bilevel false --PV_CAPA 1.6 --PVC 300 --EIC 0.9 --EEC 0.2 --DSOEC 0.1 --GCC 80 --weight_I 0.01
julia --project src/main_bilevel.jl --EV true --EHP true --storage true --Storage_cost 500 --network_reconfig true --bilevel false --PV_CAPA 1.6 --PVC 300 --EIC 0.9 --EEC 0.3 --DSOEC 0.1 --GCC 80 --weight_I 0.01
julia --project src/main_bilevel.jl --EV true --EHP true --storage true --Storage_cost 500 --network_reconfig true --bilevel false --PV_CAPA 1.6 --PVC 300 --EIC 0.9 --EEC 0.3 --DSOEC 0.1 --GCC 80 --weight_I 0.01
julia --project src/main_bilevel.jl --EV true --EHP true --storage true --Storage_cost 500 --network_reconfig true --bilevel false --PV_CAPA 1.6 --PVC 300 --EIC 0.9 --EEC 0.3 --DSOEC 0.2 --GCC 80 --weight_I 0.01
julia --project src/main_bilevel.jl --EV true --EHP true --storage true --Storage_cost 500 --network_reconfig true --bilevel false --PV_CAPA 1.6 --PVC 300 --EIC 0.9 --EEC 0.3 --DSOEC 0.3 --GCC 80 --weight_I 0.01
julia --project src/main_bilevel.jl --EV true --EHP true --storage true --Storage_cost 500 --network_reconfig true --bilevel false --PV_CAPA 1.6 --PVC 300 --EIC 0.9 --EEC 0.3 --DSOEC 0.3 --GCC 80 --weight_I 0.01
julia --project src/main_bilevel.jl --EV true --EHP true --storage true --Storage_cost 500 --network_reconfig true --bilevel false --PV_CAPA 1.6 --PVC 300 --EIC 0.9 --EEC 0.3 --DSOEC 0.3 --GCC 120 --weight_I 0.01
julia --project src/main_bilevel.jl --EV true --EHP true --storage true --Storage_cost 500 --network_reconfig true --bilevel false --PV_CAPA 1.6 --PVC 300 --EIC 0.9 --EEC 0.3 --DSOEC 0.3 --GCC 160 --weight_I 0.01
julia --project src/main_bilevel.jl --EV true --EHP true --storage true --Storage_cost 500 --network_reconfig true --bilevel false --PV_CAPA 1.6 --PVC 300 --EIC 0.9 --EEC 0.3 --DSOEC 0.3 --GCC 160 --weight_I 0.01
julia --project src/main_bilevel.jl --EV true --EHP true --storage true --Storage_cost 500 --network_reconfig true --bilevel false --PV_CAPA 1.6 --PVC 300 --EIC 0.9 --EEC 0.3 --DSOEC 0.3 --GCC 160 --weight_I 0.001
julia --project src/main_bilevel.jl --EV true --EHP true --storage true --Storage_cost 500 --network_reconfig true --bilevel false --PV_CAPA 1.6 --PVC 300 --EIC 0.9 --EEC 0.3 --DSOEC 0.3 --GCC 160 --weight_I 0.1
